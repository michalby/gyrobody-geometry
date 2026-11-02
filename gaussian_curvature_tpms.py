#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright (C) 2026 Michał Bykowski
"""Gaussian curvature of nodal gyroid and diamond surfaces (Fig. 1g)."""

import argparse
import os

import numpy as np
from skimage import measure


def gyroid(x, y, z):
    sx, sy, sz, cx, cy, cz = np.sin(x), np.sin(y), np.sin(z), np.cos(x), np.cos(y), np.cos(z)
    f = sx * cy + sy * cz + sz * cx
    grad = np.array([cx * cy - sz * sx, -sx * sy + cy * cz, -sy * sz + cz * cx])
    hess = np.array([[-sx * cy - sz * cx, -cx * sy, -cz * sx],
                     [-cx * sy, -sx * cy - sy * cz, -cy * sz],
                     [-cz * sx, -cy * sz, -sy * cz - sz * cx]])
    return f, grad, hess


def diamond(x, y, z):
    sx, sy, sz, cx, cy, cz = np.sin(x), np.sin(y), np.sin(z), np.cos(x), np.cos(y), np.cos(z)
    f = sx * sy * sz + sx * cy * cz + cx * sy * cz + cx * cy * sz
    grad = np.array([cx * sy * sz + cx * cy * cz - sx * sy * cz - sx * cy * sz,
                     sx * cy * sz - sx * sy * cz + cx * cy * cz - cx * sy * sz,
                     sx * sy * cz - sx * cy * sz - cx * sy * sz + cx * cy * cz])
    fxy = cx * cy * sz - cx * sy * cz - sx * cy * cz + sx * sy * sz
    fxz = cx * sy * cz - cx * cy * sz + sx * sy * sz - sx * cy * cz
    fyz = sx * cy * cz + sx * sy * sz - cx * cy * sz - cx * sy * cz
    hess = np.array([[-f, fxy, fxz], [fxy, -f, fyz], [fxz, fyz, -f]])
    return f, grad, hess


def gaussian_curvature(surface, x, y, z):
    """Bordered-Hessian formula: K = -det([[H, g], [g^T, 0]]) / |g|^4."""
    _, g, h = surface(x, y, z)
    m = np.zeros((4, 4))
    m[:3, :3], m[:3, 3], m[3, :3] = h, g, g
    return -np.linalg.det(m) / np.dot(g, g) ** 2


def curvature_extremes(surface, t, n=80000, seed=0):
    """Min and max K over random points moved onto F - t = 0 by Newton steps."""
    ks = []
    for p in np.random.default_rng(seed).random((n, 3)) * 2 * np.pi:
        for _ in range(7):
            f, g, _ = surface(*p)
            p = p - (f - t) * g / np.dot(g, g)
        p = np.mod(p, 2 * np.pi)
        if abs(surface(*p)[0] - t) < 1e-6:
            ks.append(gaussian_curvature(surface, *p))
    return min(ks), max(ks)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--structure", choices=["gyroid", "diamond"], required=True)
    ap.add_argument("--unit-cell-nm", type=float, required=True)
    ap.add_argument("--t-level", type=float, default=0.0)
    ap.add_argument("--out-dir", default="curvature_output")
    args = ap.parse_args()
    surface = {"gyroid": gyroid, "diamond": diamond}[args.structure]
    a, t = args.unit_cell_nm, args.t_level
    scale = (2 * np.pi / a) ** 2                  # radians^-2 -> nm^-2

    for name, k in zip(("K min", "K max"), curvature_extremes(surface, t)):
        print(f"{name} = {k * scale:.6g} nm^-2, R = {1 / np.sqrt(abs(k * scale)):.2f} nm")

    # marching-cubes mesh of one unit cell with K at each vertex
    grid = np.linspace(0, 2 * np.pi, 96, endpoint=False)
    vol = surface(*np.meshgrid(grid, grid, grid, indexing="ij"))[0] - t
    verts, faces, _, _ = measure.marching_cubes(vol, 0.0, spacing=(a / 96,) * 3)
    k_nm = np.array([gaussian_curvature(surface, *(v * 2 * np.pi / a)) for v in verts]) * scale

    os.makedirs(args.out_dir, exist_ok=True)
    path = os.path.join(args.out_dir, f"{args.structure}_a{a:g}nm_t{t:g}.ply")
    with open(path, "w") as fh:
        fh.write(f"ply\nformat ascii 1.0\nelement vertex {len(verts)}\n")
        for prop in ("x", "y", "z", "K_nm", "radius_nm"):
            fh.write(f"property float {prop}\n")
        fh.write(f"element face {len(faces)}\nproperty list uchar int vertex_indices\nend_header\n")
        for (x, y, z), k in zip(verts, k_nm):
            fh.write(f"{x} {y} {z} {k} {1 / np.sqrt(max(abs(k), 1e-30))}\n")
        for f in faces:
            fh.write(f"3 {f[0]} {f[1]} {f[2]}\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
