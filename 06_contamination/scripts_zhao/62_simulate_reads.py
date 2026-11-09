#!/usr/bin/env python3
"""Simulate paired-end reads from one genome, to measure species misassignment rates.

Reads: 150 bp pairs, insert size Normal(350, 50) (truncated at 200), positions uniform over
the genome (all sequences, weighted by length). Each read base differs from the reference
with probability 0.015 (stands in for within-species polymorphism) plus 0.003 (sequencing
error). Sequences with N in the fragment are resampled.
Usage: 62_simulate_reads.py GENOME.fna N_PAIRS SEED OUT_PREFIX
"""
import gzip
import sys

import numpy as np

fna, npairs, seed, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
rng = np.random.default_rng(seed)
seqs, name, buf = [], None, []
with open(fna) as f:
    for line in f:
        if line.startswith('>'):
            if buf:
                seqs.append(''.join(buf).upper())
            buf = []
        else:
            buf.append(line.strip())
seqs.append(''.join(buf).upper())
seqs = [s for s in seqs if len(s) > 1000]
w = np.array([len(s) for s in seqs], dtype=float)
w /= w.sum()
comp = str.maketrans('ACGTN', 'TGCAN')
B = np.array(list('ACGT'))


def mutate(s):
    a = np.array(list(s))
    hit = rng.random(len(a)) < 0.018
    a[hit] = B[rng.integers(0, 4, hit.sum())]
    return ''.join(a)


q = 'I' * 150
with gzip.open(f'{out}_1.fq.gz', 'wt') as f1, gzip.open(f'{out}_2.fq.gz', 'wt') as f2:
    n = 0
    while n < npairs:
        s = seqs[rng.choice(len(seqs), p=w)]
        ins = max(200, int(rng.normal(350, 50)))
        if len(s) <= ins:
            continue
        st = int(rng.integers(0, len(s) - ins))
        frag = s[st:st + ins]
        if 'N' in frag:
            continue
        if rng.random() < 0.5:
            frag = frag.translate(comp)[::-1]
        r1, r2 = mutate(frag[:150]), mutate(frag[-150:].translate(comp)[::-1])
        f1.write(f'@sim{n}/1\n{r1}\n+\n{q}\n')
        f2.write(f'@sim{n}/2\n{r2}\n+\n{q}\n')
        n += 1
print(out, npairs)
