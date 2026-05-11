#!/usr/bin/env python3
"""
03 — for each aligned read at the Srrm3 locus, do two things:

  (A) Walk CIGAR + reference position to classify the read as
        INCLUSION   — read alignment covers the cassette via M ops
        SKIP        — read alignment has an N (skip) op spanning the cassette
                      (i.e. the read uses the splice sites that define the
                       cassette, but skips the cassette body)
        UNINFORMATIVE — read does not span the cassette region at all,
                        OR is too short on one side to be sure

  (B) Slice the 16-bp 10x cell barcode from the read SEQUENCE by polyA
      anchoring. Match against the author's barcode whitelist (from obs.tsv)
      with up to BARCODE_MAX_HAMMING mismatches. Output the matched barcode
      (NOT the raw extraction) so downstream cluster lookup is exact.

Output: TSV with one row per aligned read:
    sample    read_id    classification    cb_extracted    cb_matched    cb_orientation

NB: this script does NOT do UMI dedup. For PSI, we want raw read counts
of inclusion vs skip — duplicate detection isn't required since the same
mRNA molecule can only manifest in one way (included or skipped), and
duplicates from PCR amplification preserve that ratio.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Iterator

import pandas as pd
import pysam

# ---------------------------------------------------------------------------
# Config (mirrors 00_config.sh)
# ---------------------------------------------------------------------------

BASE_DIR = os.environ.get("BASE_DIR", "/mnt/e/RiboSTAMP_SRRM3_GSE314176/data")
TGT_DIR = os.environ.get("TGT_DIR", f"{BASE_DIR}/targeted_psi")
TGT_ALIGN_DIR = Path(os.environ.get("TGT_ALIGN_DIR", f"{TGT_DIR}/aligned_locus"))
TGT_REF_DIR = Path(os.environ.get("TGT_REF_DIR", f"{TGT_DIR}/reference"))
TGT_PERREAD_DIR = Path(os.environ.get("TGT_PERREAD_DIR", f"{TGT_DIR}/per_read"))
META_DIR = Path(os.environ.get("META_DIR", f"{BASE_DIR}/metadata"))

CB_LEN = int(os.environ.get("CB_LEN", "16"))
UMI_LEN = int(os.environ.get("UMI_LEN", "12"))
POLYA_MIN = int(os.environ.get("POLYA_MIN", "8"))
BARCODE_MAX_HAMMING = int(os.environ.get("BARCODE_MAX_HAMMING", "1"))

# Map SRR ↔ mouse / GSM (from the SRA metadata)
SRR_TO_GSM = {
    "SRR36480452": "GSM9380801",  # mouse3
    "SRR36480453": "GSM9380800",  # mouse2
    "SRR36480454": "GSM9380799",  # mouse1
}

# ---------------------------------------------------------------------------
# Locus / anchor coordinates (read from 01's output)
# ---------------------------------------------------------------------------


def load_locus_meta() -> dict:
    meta_path = TGT_REF_DIR / "locus_origin.tsv"
    if not meta_path.exists():
        sys.exit(f"ERROR: {meta_path} not found. Run 01_extract_locus.sh first.")
    m = {}
    with open(meta_path) as f:
        next(f)  # header
        for line in f:
            k, v = line.rstrip("\n").split("\t", 1)
            m[k] = v
    return m


# ---------------------------------------------------------------------------
# (A) CIGAR-based inclusion/skip classifier
# ---------------------------------------------------------------------------

# pysam's CigarTuple uses these op codes:
#   0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
CONSUMES_REF = {0, 2, 3, 7, 8}
CONSUMES_QRY = {0, 1, 4, 7, 8}


def classify_read(read: pysam.AlignedSegment, anchor_start: int, anchor_end: int,
                  tol: int = 20) -> str:
    """
    Walk this read's CIGAR against the reference. Return one of:
        INCLUSION   — read has M-aligned bases covering >= 50% of the cassette
        SKIP        — read has an N op whose start/end is within `tol` bp of
                      the cassette boundaries (i.e. uses the cassette's
                      splice sites)
        UNINFORMATIVE — neither condition holds
    """
    if read.is_unmapped or read.cigartuples is None:
        return "UNINFORMATIVE"

    # Read must have alignment that flanks the cassette on BOTH sides for
    # the classification to be meaningful.
    rs = read.reference_start  # 0-based
    re_ = read.reference_end    # exclusive
    if rs is None or re_ is None:
        return "UNINFORMATIVE"
    if rs > anchor_start - 50 or re_ < anchor_end + 50:
        return "UNINFORMATIVE"  # doesn't span both sides

    # Walk CIGAR
    ref_pos = rs
    m_in_anchor = 0
    has_n_engulfing_cassette = False
    for op, length in read.cigartuples:
        if op == 3:  # N — splice/skip
            n_start = ref_pos
            n_end = ref_pos + length
            # SKIP detection: an N op that engulfs the cassette region.
            # This is the canonical splice signature for cassette skipping —
            # the read goes upstream-exon → big N spanning (intron + cassette
            # + intron) → downstream-exon. The N's boundaries don't need to
            # match the cassette's boundaries (the splice sites are AT the
            # flanking exons, not at the cassette itself).
            if n_start < anchor_start and n_end > anchor_end:
                has_n_engulfing_cassette = True
        if op in (0, 7, 8):  # M, =, X — counts as alignment to ref
            seg_start = ref_pos
            seg_end = ref_pos + length
            overlap = max(0, min(seg_end, anchor_end) - max(seg_start, anchor_start))
            m_in_anchor += overlap
        if op in CONSUMES_REF:
            ref_pos += length

    cassette_len = anchor_end - anchor_start
    inclusion_frac = m_in_anchor / cassette_len if cassette_len > 0 else 0.0

    # Decision: prefer the SKIP-via-engulfing-N signal first. If a single N
    # spans the cassette region, the read definitively skipped it — regardless
    # of any tiny M overlap from alignment wobble at junction boundaries.
    if has_n_engulfing_cassette and inclusion_frac < 0.5:
        return "SKIP"

    if inclusion_frac >= 0.5:
        return "INCLUSION"

    return "UNINFORMATIVE"


# ---------------------------------------------------------------------------
# (B) PolyA-anchored barcode extraction
# ---------------------------------------------------------------------------

POLYA_FWD = re.compile(r"A{%d,}" % POLYA_MIN)
POLYT_FWD = re.compile(r"T{%d,}" % POLYA_MIN)

REVCOMP = str.maketrans("ACGTNacgtn", "TGCANTGCAN")
def revcomp(s: str) -> str:
    return s.translate(REVCOMP)[::-1]


def extract_cb_candidates(seq: str) -> list[tuple[str, str]]:
    """
    Return list of (cb_candidate_16bp, orientation) candidates.

    10x Chromium 3' v3 read structure has UMI BETWEEN polyA and CB.

    Forward (R1-side) orientation:
       [R1 primer 22bp][CB:16][UMI:12][polyT(N)][cDNA rc][TSO rc]
    Reverse-complement (cDNA-side) orientation:
       [TSO][cDNA][polyA(N)][UMI rc:12][CB rc:16][R1 primer rc]

    So:
      - When we find a polyA stretch (forward read of cDNA), the CB(rc) is
        at offset UMI_LEN AFTER polyA: seq[polyA_end + UMI : polyA_end + UMI + CB]
        Then RC to get the CB in author's orientation.
      - When we find a polyT stretch (reverse read of cDNA), the CB is
        at offset UMI_LEN BEFORE polyT (going further 5'):
        seq[polyT_start - UMI - CB : polyT_start - UMI]
        Already in author's orientation, no RC needed.

    UMI is nominally 12 bp but PacBio sequencing introduces small length
    variations (indels). We try a small range of offsets ±2 around the
    nominal UMI_LEN to capture this.
    """
    candidates = []

    # Try a few offsets around the nominal UMI length to absorb indel slop.
    umi_offsets = [UMI_LEN, UMI_LEN - 1, UMI_LEN + 1, UMI_LEN - 2, UMI_LEN + 2]

    # polyA-anchored: CB(rc) sits AFTER polyA, separated by UMI rc.
    # Read sequence: [...][polyA][UMI rc][CB rc][...]
    # We extract the 16 bp at polyA_end + offset, then RC.
    for m in POLYA_FWD.finditer(seq):
        pa_end = m.end()
        for off in umi_offsets:
            cb_start = pa_end + off
            cb_end = cb_start + CB_LEN
            if cb_end > len(seq):
                continue
            cb_rc = seq[cb_start : cb_end]
            if "N" in cb_rc:
                continue
            candidates.append((revcomp(cb_rc), "polyA"))

    # polyT-anchored: CB sits BEFORE polyT, separated by UMI.
    # Read sequence: [...][CB][UMI][polyT][...]
    # We extract the 16 bp at polyT_start - offset - CB_LEN, no RC.
    for m in POLYT_FWD.finditer(seq):
        pt_start = m.start()
        for off in umi_offsets:
            cb_end = pt_start - off
            cb_start = cb_end - CB_LEN
            if cb_start < 0:
                continue
            cb = seq[cb_start : cb_end]
            if "N" in cb:
                continue
            candidates.append((cb, "polyT"))

    return candidates


# ---------------------------------------------------------------------------
# Whitelist match with hamming tolerance
# ---------------------------------------------------------------------------


def hamming(a: str, b: str) -> int:
    if len(a) != len(b):
        return max(len(a), len(b))
    return sum(1 for x, y in zip(a, b) if x != y)


class WhitelistMatcher:
    """O(1) exact match + O(N) one-hamming fallback."""

    def __init__(self, whitelist: set[str], max_hamming: int = 1):
        self.exact = whitelist
        self.max_hamming = max_hamming
        # Pre-compute a list for the linear fallback (only used when exact misses)
        self._wl_list = list(whitelist) if max_hamming >= 1 else []

    def match(self, cb: str) -> str | None:
        if cb in self.exact:
            return cb
        if self.max_hamming == 0:
            return None
        # Linear fallback
        for w in self._wl_list:
            if hamming(cb, w) <= self.max_hamming:
                return w
        return None


def load_author_barcodes_for_mouse(gsm: str) -> set[str]:
    """Load the bare 16-mer cell barcodes from the author obs.tsv for one mouse."""
    obs_paths = list(META_DIR.glob(f"{gsm}_longread_normed_counts_*__obs.tsv"))
    if not obs_paths:
        sys.exit(f"ERROR: no obs.tsv for {gsm} under {META_DIR}")
    df = pd.read_csv(obs_paths[0], sep="\t")
    # Strip the "_sampleN" suffix that the long-read obs uses
    bc = df["barcode"].astype(str).str.split("_", n=1).str[0]
    return set(bc.tolist())


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def process_one(srr: str, locus_meta: dict) -> Path:
    bam_path = TGT_ALIGN_DIR / f"{srr}.locus.bam"
    if not bam_path.exists():
        sys.exit(f"ERROR: {bam_path} not found. Run 02_align_to_locus.sh first.")

    gsm = SRR_TO_GSM.get(srr)
    if gsm is None:
        sys.exit(f"ERROR: unknown SRR {srr} (not in SRR_TO_GSM map)")
    whitelist = load_author_barcodes_for_mouse(gsm)
    matcher = WhitelistMatcher(whitelist, max_hamming=BARCODE_MAX_HAMMING)
    print(f"[{srr}] loaded {len(whitelist)} author barcodes for {gsm}")

    # Anchor coordinates in the locus FASTA's coordinate system
    anchor_start_in_locus = int(locus_meta["anchor_start_in_locus"])
    anchor_end_in_locus = int(locus_meta["anchor_end_in_locus"])

    # The locus FASTA uses a single seqname like "chr5:135698575-136098652"
    # — pysam handles that fine.

    out_path = TGT_PERREAD_DIR / f"{srr}.per_read.tsv"
    n_total = 0
    n_inclusion = 0
    n_skip = 0
    n_unin = 0
    n_cb_matched = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam, \
         open(out_path, "w") as out:
        out.write("sample\tread_id\tclassification\tcb_extracted\tcb_matched\tcb_orientation\n")
        for read in bam.fetch(until_eof=True):
            n_total += 1
            cls = classify_read(read, anchor_start_in_locus, anchor_end_in_locus)

            cb_extracted = ""
            cb_matched = ""
            cb_ori = ""
            seq = read.query_sequence
            if seq is not None:
                # Read is reverse-complemented if it aligned to negative strand —
                # restore original orientation for CB extraction.
                if read.is_reverse:
                    seq = revcomp(seq)
                for cand, ori in extract_cb_candidates(seq):
                    matched = matcher.match(cand)
                    if matched is not None:
                        cb_extracted = cand
                        cb_matched = matched
                        cb_ori = ori
                        n_cb_matched += 1
                        break

            out.write(f"{srr}\t{read.query_name}\t{cls}\t{cb_extracted}\t"
                      f"{cb_matched}\t{cb_ori}\n")

            if cls == "INCLUSION":   n_inclusion += 1
            elif cls == "SKIP":      n_skip += 1
            else:                    n_unin += 1

    print(f"[{srr}] {n_total} reads | INCLUSION={n_inclusion} SKIP={n_skip} "
          f"UNINFORMATIVE={n_unin} | CB matched={n_cb_matched} "
          f"({100*n_cb_matched/n_total:.1f}%)")
    if n_inclusion + n_skip > 0:
        print(f"[{srr}] bulk PSI = {n_inclusion / (n_inclusion + n_skip):.3f}")
    return out_path


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--srr", nargs="*", default=list(SRR_TO_GSM.keys()),
                   help="SRR accessions to process (default: all 3)")
    args = p.parse_args()

    locus_meta = load_locus_meta()
    print(f"Anchor (locus FASTA coords): {locus_meta['anchor_start_in_locus']}-"
          f"{locus_meta['anchor_end_in_locus']}")
    print(f"Anchor (mm10 coords): {locus_meta['anchor_start_mm10']}-"
          f"{locus_meta['anchor_end_mm10']}  strand {locus_meta['anchor_strand']}")
    print()

    out_paths = []
    for srr in args.srr:
        out_paths.append(process_one(srr, locus_meta))

    # Pool all per-read TSVs into a single file for downstream
    pooled = TGT_PERREAD_DIR / "all_samples.per_read.tsv"
    dfs = [pd.read_csv(p, sep="\t") for p in out_paths]
    pd.concat(dfs, ignore_index=True).to_csv(pooled, sep="\t", index=False)
    print(f"\nWrote pooled per-read TSV: {pooled}  ({sum(len(d) for d in dfs)} rows)")


if __name__ == "__main__":
    main()
