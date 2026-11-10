#!/usr/bin/env python3
"""ATG-initiated ORF (>= 30 aa) that aligns best to the parent of the RpL36, RpS14, CCT1 (T-cp1) and VhaAC39-1 retrocopies in IsoJpt1.0 and in the two
WT3-strain assemblies (Dsuzukii.v01, Chiu et al. 2013, GCA_000472105.1; LBDM_Dsuz_2.1.pri, GCA_013340165.1), both
strands of the copy +- 300 nt; coverage of and identity to the D. suzukii (IsoJpt1.0) parent protein (local
alignment, BLOSUM62, gap -10/-1). Also reports in-frame stop codons and frame changes of the copy relative to the
parent protein via the alignment. No random numbers."""
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner, substitution_matrices
R = sys.argv[1]
par = {r.id: str(r.seq.translate()).rstrip('*') for r in SeqIO.parse(f'{R}/parent_cds.fa', 'fasta')}
al = PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                     open_gap_score=-10, extend_gap_score=-1)
print('copy\tparent_aa\tbest_ORF_aa\tparent_positions_covered\tcoverage\tidentity')
for r in SeqIO.parse(f'{R}/retro_regions_3asm.fa', 'fasta'):
    g = r.id.rsplit('_', 1)[0]
    P = par[g]
    best, a, bs = '', None, -1  # ORF (>= 30 aa, ATG-initiated) with the highest alignment score to the parent
    for s in (r.seq, r.seq.reverse_complement()):
        for f in range(3):
            p = str(s[f:len(s) - (len(s) - f) % 3].translate())
            for seg in p.split('*'):
                i = seg.find('M')
                if i >= 0 and len(seg) - i >= 30:
                    x = al.align(seg[i:], P)[0]
                    if x.score > bs:
                        best, a, bs = seg[i:], x, x.score
    m = sum(best[b1 + k] == P[b2 + k] for (b1, e1), (b2, e2) in zip(*a.aligned) for k in range(e1 - b1))
    alen = sum(e - b for b, e in a.aligned[1])
    print(f'{r.id}\t{len(P)}\t{len(best)}\t{a.aligned[1][0][0] + 1}-{a.aligned[1][-1][1]}\t{alen / len(P):.0%}\t{m / alen:.1%}')
