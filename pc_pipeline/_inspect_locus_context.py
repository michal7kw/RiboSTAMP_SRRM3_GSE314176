#!/usr/bin/env python3
"""Manually inspect the locus FASTA around the cassette to understand whether
the 'failures' from _verify_liftover.py are real liftOver bugs or script
artifacts (wrong RC convention, off-by-one boundaries, etc.).
"""
from pathlib import Path
import pyfaidx

LOCUS_FA = Path("/mnt/e/RiboSTAMP_SRRM3_GSE314176/data/targeted_psi/reference/srrm3_locus_mm10.fa")

LAB_SEQ = "GAACCTCCTGCTGCGTGGCCTGTGGCCAGCCTGGGGGCTGTGTGGACCGGGCTCCACAGGCAGCCCAAGTGATACTCAT"

REVCOMP = str.maketrans("ACGT", "TGCA")
def revcomp(s): return s.translate(REVCOMP)[::-1]

fa = pyfaidx.Fasta(str(LOCUS_FA))
seqname = list(fa.keys())[0]
seq = str(fa[seqname][:]).upper()
print(f"Locus FASTA: {seqname}, length {len(seq):,} bp\n")

print(f"Lab published cassette sequence (79 bp):")
print(f"  {LAB_SEQ}")
print(f"\nLab sequence reverse-complement (would be mRNA strand IF lab printed + strand):")
print(f"  {revcomp(LAB_SEQ)}")

print("\n" + "=" * 80)
print("DIRECT (NO RC) COMPARISON: locus + strand vs lab sequence")
print("=" * 80)

# Try our recorded boundaries
ours_30000_30078 = seq[30000:30078]
print(f"\nOur locus[30000:30078] (78 bp, + strand):")
print(f"  {ours_30000_30078}")

# Compare to lab seq
matches_full = sum(1 for a, b in zip(LAB_SEQ, ours_30000_30078) if a == b)
matches_drop_first = sum(1 for a, b in zip(LAB_SEQ[1:], ours_30000_30078) if a == b)
matches_drop_last = sum(1 for a, b in zip(LAB_SEQ[:-1], ours_30000_30078) if a == b)
print(f"\nDirect overlap with lab seq:")
print(f"  vs lab[0:78]:   {matches_full}/78")
print(f"  vs lab[1:79]:   {matches_drop_first}/78  (drop lab's first base)")
print(f"  vs lab[0:78]:   {matches_drop_last}/78  (drop lab's last base — same as full)")

# Also try shifted boundaries
print(f"\nTrying various 79-bp windows in our locus:")
for offset in range(-3, 4):
    start = 30000 + offset - (1 if offset < 0 else 0)
    end = start + 79
    win = seq[start:end]
    matches = sum(1 for a, b in zip(LAB_SEQ, win) if a == b)
    print(f"  locus[{start}:{end}] (offset {offset:+d}, 79 bp):  matches lab = {matches}/79")
    if matches == 79:
        print(f"    *** EXACT 79-bp match starting at locus position {start} ***")

print("\n" + "=" * 80)
print("WIDER CONTEXT — 30 bp on each side of recorded cassette boundaries")
print("=" * 80)
ctx_left = seq[30000-30:30000]
ctx_right = seq[30078:30078+30]
print(f"\n30 bp upstream of position 30000 on + strand:")
print(f"  {ctx_left}")
print(f"30 bp downstream of position 30078 on + strand:")
print(f"  {ctx_right}")

print(f"\nFull cassette region [29990:30090] with markers:")
context = seq[30000-10:30078+10]
marker = " " * 10 + "[" + " " * 76 + "]"
print(f"  {context}")
print(f"  {marker}")

print("\n" + "=" * 80)
print("SPLICE SITE CHECK — try multiple boundary interpretations")
print("=" * 80)

# Recall: − strand cassette. For canonical GT/AG splicing on the mRNA:
#  - intron BEFORE cassette in mRNA direction: ends in AG
#  - intron AFTER cassette in mRNA direction:  starts with GT
# Translating to + strand (cassette is on −, so mRNA = RC of + strand at cassette region):
#  - "intron BEFORE cassette in mRNA direction" = intron AFTER cassette in + direction
#    (because reading + L→R is reading mRNA R→L for − strand genes)
#  - "intron AFTER cassette in mRNA direction" = intron BEFORE cassette in + direction
# So on the + strand:
#  - intron AFTER cassette (in + direction) — its mRNA-3' end "AG" appears as "CT" reading + L→R
#    NO WAIT. The intron AFTER cassette on + strand IS itself; we read it L→R on +. Its mRNA
#    direction reading is RIGHT→LEFT for − strand. So the intron's mRNA-3' end (which is "AG"
#    on mRNA) is at the LEFT end of the intron on +. Read on + L→R, the leftmost 2 bp of the
#    intron (immediately after the cassette on +) are the RC of "AG" reversed = RC of mRNA-3'
#    end = "CT"... wait no.
#
# Just do it empirically:
#   - take 4 bp on each side of cassette
#   - take RC of the cassette region + flanks
#   - check the RC for canonical AG-cassette-GT
#   - i.e. (RC of right_flank) + (RC of cassette) + (RC of left_flank) is the mRNA-strand seq

left_flank4 = seq[30000-4:30000]
cassette_plus = seq[30000:30078]
right_flank4 = seq[30078:30078+4]

mrna_left = revcomp(right_flank4)   # mRNA-direction LEFT of cassette = RC of + strand RIGHT
mrna_cas  = revcomp(cassette_plus)
mrna_right = revcomp(left_flank4)   # mRNA-direction RIGHT of cassette = RC of + strand LEFT

print(f"\n+ strand context (our boundaries 30000:30078):")
print(f"  ...{left_flank4}|{cassette_plus[:6]}...{cassette_plus[-6:]}|{right_flank4}...")
print(f"\nmRNA strand (RC) reading 5' → 3':")
print(f"  ...{mrna_left}|{mrna_cas[:6]}...{mrna_cas[-6:]}|{mrna_right}...")
print(f"\nFor canonical splicing on mRNA strand, we expect:")
print(f"  {' '*4}AG | (cassette) | GT")
print(f"\nObserved:")
print(f"  Last 2 bp BEFORE cassette on mRNA: '{mrna_left[-2:]}'  ", "✓ canonical AG" if mrna_left[-2:] == "AG" else "✗ not AG")
print(f"  First 2 bp AFTER cassette on mRNA: '{mrna_right[:2]}'  ", "✓ canonical GT" if mrna_right[:2] == "GT" else "✗ not GT")

# Try shifted boundaries to find the canonical splice sites
print(f"\nTrying 1-bp boundary shifts to find canonical AG..GT:")
for d_start in range(-3, 4):
    for d_end in range(-3, 4):
        s0 = 30000 + d_start
        e0 = 30078 + d_end
        if e0 - s0 < 50:
            continue
        lf = seq[s0-4:s0]
        rf = seq[e0:e0+4]
        ml = revcomp(rf)
        mr = revcomp(lf)
        if ml[-2:] == "AG" and mr[:2] == "GT":
            cas_len = e0 - s0
            print(f"  Cassette [{s0}:{e0}] (length {cas_len}):  "
                  f"upstream...{ml[-4:]}|cassette|{mr[:4]}...downstream  ✓ AG..GT")
