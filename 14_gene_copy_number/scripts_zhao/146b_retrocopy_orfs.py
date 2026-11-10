#!/usr/bin/env python3
"""Longest ATG-initiated ORF in each D. suzukii retrocopy (IsoJpt1.0 region +- 300 nt) and its coverage of the
parent protein (D. suzukii parent CDS translated; local protein alignment, BLOSUM62). No random numbers."""
import sys
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
R = sys.argv[1]
par = {r.id: str(r.seq.translate()).rstrip('*') for r in SeqIO.parse(f'{R}/parent_cds.fa', 'fasta')}
al = PairwiseAligner(mode='local', substitution_matrix=substitution_matrices.load('BLOSUM62'), open_gap_score=-10, extend_gap_score=-1)
strand = {'RpL36': '+', 'RpS14': '-', 'CCT1': '-', 'VhaAC39-1': '-'}
pname = {'RpS14': 'RpS14', 'RpL36': 'RpL36', 'CCT1': 'CCT1', 'VhaAC39-1': 'VhaAC39-1'}
for r in SeqIO.parse(f'{R}/retro_regions.fa', 'fasta'):
    s = r.seq if strand[r.id] == '+' else r.seq.reverse_complement()
    best = ''
    for f in range(3):
        p = str(s[f:len(s) - (len(s) - f) % 3].translate())
        for seg in p.split('*'):
            i = seg.find('M')
            if i >= 0 and len(seg) - i > len(best):
                best = seg[i:]
    P = par[pname[r.id]]
    a = al.align(best, P)[0]
    ident = sum(x == y for x, y in zip(*[str(z) for z in a])) if False else None
    cov = sum(e - b for b, e in a.aligned[1]) / len(P)
    m = sum(best[b1 + k] == P[b2 + k] for (b1, e1), (b2, e2) in zip(*a.aligned) for k in range(e1 - b1))
    alen = sum(e - b for b, e in a.aligned[1])
    print(f'{r.id}\tparent {len(P)} aa\tlongest ORF {len(best)} aa\tcovers {cov:.0%} of parent (positions {a.aligned[1][0][0]+1}-{a.aligned[1][-1][1]})\tidentity {m/alen:.1%}')
