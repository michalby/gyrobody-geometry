#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Michał Bykowski
"""Gyroid unit-cell membrane area, gyrobody and granum membrane areas,
and grana equivalent to one gyrobody (Methods; Fig. 2c)."""

import argparse
import csv
import math
from collections import defaultdict


def mean(values):
    values = [float(v) for v in values if v.strip()]
    return sum(values) / len(values)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lattice-nm", type=float, nargs="+", default=[424.0, 551.0])
    ap.add_argument("--lumen-nm", type=float, default=30.0)
    ap.add_argument("--gyrobodies")
    ap.add_argument("--grana")
    ap.add_argument("--unit-cell-nm", type=float, default=519.0)
    ap.add_argument("--unit-cell-area-um2", type=float, default=1.62)
    args = ap.parse_args()

    # A_unit = 2 (3.019 a^2 + 2 pi chi d^2), chi = -8, two membranes
    d = args.lumen_nm / 2
    for a in args.lattice_nm:
        a_unit = 2 * (3.019 * a ** 2 + 2 * math.pi * -8 * d ** 2) / 1e6
        print(f"A_unit(a = {a:g} nm) = {a_unit:.2f} um^2")

    if args.gyrobodies:
        v = defaultdict(float)
        for r in csv.DictReader(open(args.gyrobodies, encoding="utf-8-sig")):
            l1, l2 = float(r["gyrobody_major_nm"]), float(r["gyrobody_minor_nm"])
            v[r["chloroplast_id"]] += math.pi / 6 * l1 * l2 * max(l1, l2) / 1e9
        v_gb = sum(v.values()) / len(v)
        n_uc = v_gb / (args.unit_cell_nm / 1e3) ** 3
        s_gb = n_uc * args.unit_cell_area_um2
        print(f"mean V_gb = {v_gb:.3f} um^3 ({len(v)} chloroplasts), "
              f"unit cells = {n_uc:.1f}, S_gb = {s_gb:.1f} um^2")

    if args.grana:
        rows = list(csv.DictReader(open(args.grana, encoding="utf-8-sig")))
        dia = mean(r["granum_diameter_nm"] for r in rows)
        n_th = mean(r["n_thylakoids"] for r in rows)
        s_gs = 2 * n_th * math.pi * (dia / 2) ** 2 / 1e6
        print(f"mean granum diameter = {dia:.1f} nm, thylakoids = {n_th:.2f}, "
              f"S_gs = {s_gs:.2f} um^2")
        if args.gyrobodies:
            print(f"grana per gyrobody = S_gb / S_gs = {s_gb / s_gs:.2f}")


if __name__ == "__main__":
    main()
