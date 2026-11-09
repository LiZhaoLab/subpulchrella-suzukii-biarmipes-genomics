#!/usr/bin/env python3
"""Depth histograms and thresholds per sample and chromosome class.

Usable depth = A+C+G+T reads after MAPQ >= 20, base quality >= 20, proper pairs,
duplicates removed (per-site counts from 01_pileup_counts.py). Among covered sites
(depth > 0): minimum = max(10, 5th percentile); maximum = 98th and 99th percentiles.
Writes results/depth_thresholds.tsv.
"""
import glob

import numpy as np
import pandas as pd

A = '/ru-auth/local/home/asato/scratch/population_polymorphism/reanalysis_20260821/config/regions'
QC = '/ru-auth/local/home/lzhao/Data_scratch/zhao/suzukii_QC'
OUT = '/ru-auth/local/home/lzhao/vscode_projects/suzukii/results/depth_thresholds.tsv'
SAMPLES = {'Dsub_L5_7': 'Dsub', 'JP-Kan': 'Dsuz', 'JP-Sap': 'Dsuz', 'JP-Sap2': 'Dsuz', 'JP-Tok': 'Dsuz'}


def pct(hist, q):
    c = np.cumsum(hist)
    return int(np.searchsorted(c, q / 100 * c[-1]))


rows = []
for sample, sp in SAMPLES.items():
    xchrom = open(f'{A}/{sp}.X.bed').read().split()[0]
    hist = {'autosome': np.zeros(65536, dtype=np.int64), 'X': np.zeros(65536, dtype=np.int64)}
    files = sorted(glob.glob(f'{QC}/counts/{sample}/*.npz'))
    for f in files:
        d = np.load(f)
        depth = d['counts'].sum(axis=1, dtype=np.int64)
        cls = 'X' if str(d['chrom']) == xchrom else 'autosome'
        hist[cls] += np.bincount(np.minimum(depth, 65535), minlength=65536)
    for cls, h in hist.items():
        cov = h.copy()
        cov[0] = 0
        rows.append(dict(sample=sample, cls=cls, chunks=len(files), total_sites=int(h.sum()),
                         covered_sites=int(cov.sum()),
                         mean_depth_covered=float((np.arange(65536) * cov).sum() / cov.sum()),
                         median=pct(cov, 50), p05=pct(cov, 5), min_depth=max(10, pct(cov, 5)),
                         p98=pct(cov, 98), p99=pct(cov, 99)))
        np.save(f'{QC}/depth_hist_{sample}_{cls}.npy', h)
t = pd.DataFrame(rows)
t.to_csv(OUT, sep='\t', index=False)
print(t.to_string())
