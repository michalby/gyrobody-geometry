#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Michał Bykowski
"""Gyrobody detection probability p in 1 µm sections (Supplementary Fig. 9)."""

import argparse
import csv
from collections import defaultdict

import numpy as np


def mean_apparent_k(kappa, position, rng, n):
    """Mean gyrobody/chloroplast width ratio seen in sections that hit the
    gyrobody (lengths in units of the chloroplast semi-axis)."""
    s = rng.uniform(-1, 1, n)                      # section offset
    if position == "central":
        d = np.zeros(n)
    else:                                          # random gyrobody centre
        r = 1 - kappa
        d = rng.uniform(-r, r, n)
        keep = rng.uniform(0, 1, n) < 1 - (d / r) ** 2
        s, d = s[keep], d[keep]
    u = s - d                                      # offset from gyrobody centre
    hit = np.abs(u) < kappa
    s, u = s[hit], u[hit]
    return np.mean(kappa * np.sqrt(1 - (u / kappa) ** 2) / np.sqrt(1 - s ** 2))


def calibrate_kappa(mean_k, position, rng, n=200000):
    """Bisection for kappa such that the simulated mean apparent k equals mean_k."""
    lo, hi = mean_k, 0.999
    for _ in range(42):
        mid = (lo + hi) / 2
        if mean_apparent_k(mid, position, rng, n) < mean_k:
            lo = mid
        else:
            hi = mid
    return min(1, (lo + hi) / 2)


def detection_probability(major_chloro, major_gyro, position, t, rng):
    k = np.clip(major_gyro / major_chloro, 0, 1)
    a = major_chloro / 2
    factor = calibrate_kappa(k.mean(), position, rng) / k.mean()
    kappa_i = np.minimum(1, factor * k)
    return np.mean((kappa_i * a + t / 2) / (a + t / 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tem-axes", required=True)
    ap.add_argument("--lm-counts")
    ap.add_argument("--position", choices=["central", "random"], default="central")
    ap.add_argument("--thickness-nm", type=float, default=1000.0)
    ap.add_argument("--replicates", type=int, default=10)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    axes = defaultdict(list)
    for r in csv.DictReader(open(args.tem_axes, encoding="utf-8-sig")):
        axes[r["genotype"]].append((float(r["chloroplast_major_nm"]),
                                    float(r["gyrobody_major_nm"])))
    counts = defaultdict(lambda: [0.0, 0.0])
    if args.lm_counts:
        for r in csv.DictReader(open(args.lm_counts, encoding="utf-8-sig")):
            counts[r["genotype"]][0] += float(r["chloroplasts_with_gyrobody"])
            counts[r["genotype"]][1] += float(r["chloroplasts_total"])

    print(f"position: {args.position}, t = {args.thickness_nm:g} nm")
    print("genotype\tn_TEM\tmean_k\tp\tp_MC_SD\tobserved")
    for g, rows in axes.items():
        mc, mg = np.array(rows).T
        p = [detection_probability(mc, mg, args.position, args.thickness_nm,
                                   np.random.default_rng(args.seed + i))
             for i in range(args.replicates)]
        w, n = counts[g]
        obs = f"{w / n:.3f}" if n else "-"
        print(f"{g}\t{len(mc)}\t{np.mean(np.clip(mg / mc, 0, 1)):.3f}\t{np.mean(p):.3f}\t"
              f"{np.std(p):.4f}\t{obs}")


if __name__ == "__main__":
    main()
