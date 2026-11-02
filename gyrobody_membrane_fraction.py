#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Michał Bykowski
"""Gyrobody share of the stroma-region membrane per chloroplast (Supplementary Fig. 12e)."""

import argparse
import csv
import math
from collections import defaultdict


def ellipsoid_volume(l1, l2):
    """(pi/6) l1 l2 l3 with l3 = the longer measured axis."""
    return math.pi / 6 * l1 * l2 * max(l1, l2)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chloroplasts", required=True)
    ap.add_argument("--gyrobodies", required=True)
    ap.add_argument("--stroma-traces", required=True)
    ap.add_argument("--unit-cell-nm", type=float, default=519.0)
    ap.add_argument("--unit-cell-area-um2", type=float, default=1.62)
    args = ap.parse_args()
    v_uc = (args.unit_cell_nm / 1e3) ** 3                      # um^3

    v_gb, a_gb = defaultdict(float), defaultdict(float)        # um^3, um^2
    for r in csv.DictReader(open(args.gyrobodies, encoding="utf-8-sig")):
        c = r["chloroplast_id"]
        v_gb[c] += ellipsoid_volume(float(r["gyrobody_major_nm"]),
                                    float(r["gyrobody_minor_nm"])) / 1e9
        a_gb[c] += float(r["gyrobody_area_nm2"]) / 1e6
    L = defaultdict(float)                                     # um
    for r in csv.DictReader(open(args.stroma_traces, encoding="utf-8-sig")):
        L[r["chloroplast_id"]] += float(r["length_nm"]) / 1e3

    print("chloroplast_id\tgenotype\tV_cp_um3\tV_gb_um3\tS_st_um2\tS_gb_um2\tf_gb_pct")
    for r in csv.DictReader(open(args.chloroplasts, encoding="utf-8-sig")):
        c = r["chloroplast_id"]
        v_cp = ellipsoid_volume(float(r["chloroplast_major_nm"]),
                                float(r["chloroplast_minor_nm"])) / 1e9
        a_cp = float(r["chloroplast_area_nm2"]) / 1e6
        s_st = 2 * L[c] * (v_cp - v_gb[c]) / (a_cp - a_gb[c])
        s_gb = args.unit_cell_area_um2 * v_gb[c] / v_uc
        f_gb = 100 * s_gb / (s_gb + s_st)
        print(f"{c}\t{r['genotype']}\t{v_cp:.3f}\t{v_gb[c]:.3f}\t"
              f"{s_st:.2f}\t{s_gb:.2f}\t{f_gb:.3f}")


if __name__ == "__main__":
    main()
