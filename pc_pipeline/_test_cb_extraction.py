#!/usr/bin/env python3
"""Diagnostic: how often does the polyA-anchored CB extractor generate
ANY candidate from a sample of reads, regardless of whitelist match?"""
import re
import sys

import pandas as pd
import pysam

BAM = "/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/aligned_locus/SRR36480452.locus.bam"
WHITELIST = "/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/metadata/GSM9380801_longread_normed_counts_transcript_adata_mouse3__obs.tsv"

POLYA_FWD = re.compile(r"A{8,}")
POLYT_FWD = re.compile(r"T{8,}")
CB_LEN = 16

REVCOMP = str.maketrans("ACGTNacgtn", "TGCANTGCAN")
def revcomp(s): return s.translate(REVCOMP)[::-1]


def extract_candidates(seq):
    cands = []
    for m in POLYA_FWD.finditer(seq):
        if m.end() + CB_LEN > len(seq): continue
        cb = seq[m.end():m.end() + CB_LEN]
        if "N" in cb: continue
        cands.append((cb, "fwd"))
    for m in POLYT_FWD.finditer(seq):
        if m.start() - CB_LEN < 0: continue
        cb_rc = seq[m.start() - CB_LEN:m.start()]
        if "N" in cb_rc: continue
        cands.append((revcomp(cb_rc), "rc"))
    return cands


def main():
    wl_df = pd.read_csv(WHITELIST, sep="\t")
    whitelist = set(wl_df["barcode"].astype(str).str.split("_", n=1).str[0])
    print(f"Whitelist: {len(whitelist)} barcodes")

    n_total = 0
    n_with_polyA = 0
    n_with_polyT = 0
    n_with_any_candidate = 0
    n_match = 0
    polyA_lengths = []
    polyT_lengths = []
    sample_seqs = []

    with pysam.AlignmentFile(BAM, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            n_total += 1
            seq = read.query_sequence
            if seq is None: continue
            if read.is_reverse:
                seq = revcomp(seq)

            polyA_matches = list(POLYA_FWD.finditer(seq))
            polyT_matches = list(POLYT_FWD.finditer(seq))
            if polyA_matches:
                n_with_polyA += 1
                polyA_lengths.append(max(m.end() - m.start() for m in polyA_matches))
            if polyT_matches:
                n_with_polyT += 1
                polyT_lengths.append(max(m.end() - m.start() for m in polyT_matches))

            cands = extract_candidates(seq)
            if cands:
                n_with_any_candidate += 1
                # Try to match
                for cand, ori in cands:
                    if cand in whitelist:
                        n_match += 1
                        break

            if n_total <= 5:
                sample_seqs.append((read.query_name, seq[:60], seq[-60:], len(seq), cands[:3]))

            if n_total >= 5000: break

    print(f"\nFrom first {n_total} reads:")
    print(f"  reads with >=1 polyA stretch: {n_with_polyA} ({100*n_with_polyA/n_total:.1f}%)")
    print(f"  reads with >=1 polyT stretch: {n_with_polyT} ({100*n_with_polyT/n_total:.1f}%)")
    print(f"  reads with >=1 CB candidate:  {n_with_any_candidate} ({100*n_with_any_candidate/n_total:.1f}%)")
    print(f"  reads with EXACT match:        {n_match} ({100*n_match/n_total:.1f}%)")

    if polyA_lengths:
        polyA_lengths.sort()
        print(f"\n  polyA length stats: min={polyA_lengths[0]} median={polyA_lengths[len(polyA_lengths)//2]} max={polyA_lengths[-1]}")
    if polyT_lengths:
        polyT_lengths.sort()
        print(f"  polyT length stats: min={polyT_lengths[0]} median={polyT_lengths[len(polyT_lengths)//2]} max={polyT_lengths[-1]}")

    print("\nSample reads:")
    for name, head, tail, length, cands in sample_seqs:
        print(f"  {name[:50]}  len={length}")
        print(f"    head: {head}")
        print(f"    tail: {tail}")
        for c, o in cands:
            print(f"    cand({o}): {c}  in_whitelist={c in whitelist}")


if __name__ == "__main__":
    main()
