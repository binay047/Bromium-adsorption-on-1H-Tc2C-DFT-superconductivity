#!/usr/bin/env python3
"""
export_fs_dat.py

Export tiled 2D Fermi-surface contours (kz=0 slice) from a QE BXSF file,
plus the Gamma-M-K-Gamma high-symmetry path, as plain xmgrace-ready .dat
files. This mirrors the tiling/placement logic used in the matplotlib
preview that already matched the paper's look -- just writing .dat files
instead of a PNG.

Usage:
    python3 export_fs_dat.py aiida_fs.bxsf --tile 3 --outdir fs_dat

Writes into --outdir:
    band16_fs.dat, band17_fs.dat, ...   (one file per band crossing E_F,
                                          each contour piece as its own
                                          xmgrace set, separated by '&')
    gmk_path.dat                        (Gamma-M-K-Gamma, in the same
                                          [0,1]x[0,1] tiled coordinate
                                          system as the FS contours)
"""

import argparse
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_bxsf(path):
    with open(path, "r") as f:
        text = f.read()

    m = re.search(r"Fermi Energy:\s*([-\d.Ee+]+)", text, re.IGNORECASE)
    if m is None:
        m = re.search(r"FERMI\s*ENERGY\s*[:=]?\s*([-\d.Ee+]+)", text, re.IGNORECASE)
    ef = float(m.group(1)) if m else None

    lines = text.splitlines()

    def find_line(keyword, start=0):
        for i in range(start, len(lines)):
            if keyword in lines[i]:
                return i
        raise ValueError(f"Could not find '{keyword}' in BXSF file")

    idx = find_line("BANDGRID_3D_BANDS")
    n1, n2, n3 = map(int, lines[idx + 2].split()[:3])
    origin = np.array(list(map(float, lines[idx + 3].split()[:3])))
    span1 = np.array(list(map(float, lines[idx + 4].split()[:3])))
    span2 = np.array(list(map(float, lines[idx + 5].split()[:3])))
    span3 = np.array(list(map(float, lines[idx + 6].split()[:3])))

    bands = {}
    i = idx + 7
    while i < len(lines):
        if "BAND:" in lines[i]:
            band_idx = int(lines[i].split()[-1])
            i += 1
            vals = []
            while i < len(lines) and "BAND:" not in lines[i] and "END_BANDGRID_3D" not in lines[i]:
                vals.extend(float(x) for x in lines[i].split())
                i += 1
            arr = np.array(vals[: n1 * n2 * n3]).reshape(n1, n2, n3)
            bands[band_idx] = arr
        elif "END_BANDGRID_3D" in lines[i]:
            break
        else:
            i += 1

    return dict(ef=ef, origin=origin, span1=span1, span2=span2, span3=span3,
                n1=n1, n2=n2, n3=n3, bands=bands)


def bands_crossing_ef(bands, ef, tol=1e-6):
    return sorted(i for i, arr in bands.items()
                  if arr.min() - tol <= ef <= arr.max() + tol)


def gamma_fractional_position(data):
    origin_xy = data["origin"][:2]
    M = np.column_stack([data["span1"][:2], data["span2"][:2]])
    return np.linalg.solve(M, -origin_xy)


def high_symmetry_fractional(data):
    g = gamma_fractional_position(data)
    return {
        "Gamma": g,
        "M": g + np.array([0.5, 0.0]),
        "K": g + np.array([1 / 3, 1 / 3]),
    }


def tile_grid(E, n_tile):
    return np.tile(E, (n_tile, n_tile))


def write_xmgrace_multiset(path, blocks, header_comment="", include_legend=False):
    with open(path, "w") as f:
        if header_comment:
            f.write(f"# {header_comment}\n")
        for i, (label, arr) in enumerate(blocks):
            if include_legend:
                f.write(f"@    s{i} legend \"{label}\"\n")
            else:
                f.write(f"# {label}\n")
            for row in arr:
                f.write(" ".join(f"{v:.6f}" for v in row) + "\n")
            f.write("&\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bxsf")
    p.add_argument("--tile", type=int, default=3)
    p.add_argument("--kz-index", type=int, default=0)
    p.add_argument("--outdir", default="fs_dat")
    p.add_argument("--with-legend", action="store_true",
                    help="embed @sN legend commands in the .dat files "
                         "(off by default, since xmgrace auto-displays "
                         "them as a legend box on load)")
    args = p.parse_args()

    data = parse_bxsf(args.bxsf)
    ef = data["ef"]
    n1, n2 = data["n1"], data["n2"]
    N = args.tile

    print(f"Fermi energy: {ef:.4f} eV")
    crossing = bands_crossing_ef(data["bands"], ef)
    print(f"Bands crossing E_F: {crossing}")

    os.makedirs(args.outdir, exist_ok=True)

    kx_tiled = np.linspace(0, N, N * n1, endpoint=False) / N
    ky_tiled = np.linspace(0, N, N * n2, endpoint=False) / N
    KX, KY = np.meshgrid(kx_tiled, ky_tiled, indexing="ij")

    for band_idx in crossing:
        E = data["bands"][band_idx][:, :, args.kz_index] - ef
        E_tiled = tile_grid(E, N)

        cs = plt.contour(KX, KY, E_tiled, levels=[0])
        segs = cs.allsegs[0]
        plt.close()

        outfile = os.path.join(args.outdir, f"band{band_idx}_fs.dat")
        blocks = [(f"band{band_idx}_piece{k}", seg) for k, seg in enumerate(segs)]
        if not blocks:
            print(f"  band {band_idx}: no contour at E_F (skipped)")
            continue
        write_xmgrace_multiset(
            outfile, blocks,
            header_comment=f"Fermi surface, band index {band_idx}, "
                            f"tiled {N}x{N}, E_F={ef:.4f} eV",
            include_legend=args.with_legend
        )
        print(f"  band {band_idx}: {len(segs)} piece(s) -> {outfile}")

    pts = high_symmetry_fractional(data)
    path_order = ["Gamma", "M", "K", "Gamma"]
    path_xy = np.array([pts[p] for p in path_order]) / N

    path_file = os.path.join(args.outdir, "gmk_path.dat")
    write_xmgrace_multiset(
        path_file, [("Gamma-M-K-Gamma", path_xy)],
        header_comment="High-symmetry path, same [0,1]x[0,1] tiled coordinate "
                        "system as the FS contours",
        include_legend=args.with_legend
    )
    print(f"  path -> {path_file}")

    points_file = os.path.join(args.outdir, "high_sym_points.dat")
    with open(points_file, "w") as f:
        f.write("# label   kx           ky\n")
        for name in ["Gamma", "M", "K"]:
            x, y = pts[name] / N
            f.write(f"{name:6s} {x:.6f} {y:.6f}\n")
    print(f"  labelled points -> {points_file}")


if __name__ == "__main__":
    main()
