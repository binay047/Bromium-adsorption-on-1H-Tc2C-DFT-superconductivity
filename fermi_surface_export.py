#!/usr/bin/env python3
"""
fermi_surface_export.py

One script, three outputs, from a single QE BXSF file:
  1. fermi_surface.png   -- matplotlib preview (hexagonal BZ, velocity-
                             colored contours, Gamma-M-K-Gamma path)
  2. fermi_surface.dat   -- plain multi-column xmgrace-readable data
                             (velocity-binned sets, colors NOT applied --
                             for manual/other use)
  3. fermi_surface.agr   -- ready-to-open Grace PROJECT file with colors
                             already assigned per set. Just run:
                                 xmgrace fermi_surface.agr
                             and the velocity coloring shows immediately,
                             no manual color assignment needed.

------------------------------------------------------------------
USAGE

    python3 fermi_surface_export.py aiida_fs.bxsf \
        --a1 2.9426199818 0.0 \
        --a2 -1.4713099909 2.5483836579 \
        --alat-bohr 5.560746 \
        --outdir fs_out

If --a1/--a2/--alat-bohr are omitted, velocities are RELATIVE only.
------------------------------------------------------------------
"""

import argparse
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

HBAR_EVS = 6.582119569e-16
ANG_TO_M = 1e-10
N_VBINS = 12


# ---------------------------------------------------------------- parsing
def parse_bxsf(path):
    with open(path, "r") as f:
        text = f.read()

    m = re.search(r"Fermi Energy:\s*([-\d.Ee+]+)", text, re.IGNORECASE)
    ef = float(m.group(1)) if m else None
    if ef is None:
        raise ValueError("Could not find Fermi energy in BXSF file.")

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


def physical_scale_factor(a1_ang, a2_ang, alat_bohr, span1_xy):
    A = np.array([a1_ang, a2_ang])
    B_ang = 2 * np.pi * np.linalg.inv(A).T
    return np.linalg.norm(B_ang[0]) / np.linalg.norm(span1_xy)


def high_symmetry_cart(data):
    g = gamma_fractional_position(data)
    span1_xy, span2_xy = data["span1"][:2], data["span2"][:2]
    origin_xy = data["origin"][:2]
    gamma_cart = origin_xy + g[0] * span1_xy + g[1] * span2_xy
    return {
        "Gamma": gamma_cart,
        "M": gamma_cart + 0.5 * span1_xy,
        "K": gamma_cart + (1 / 3) * span1_xy + (1 / 3) * span2_xy,
    }, gamma_cart


def ws_hexagon(data, gamma_cart):
    from scipy.spatial import Voronoi
    span1_xy, span2_xy = data["span1"][:2], data["span2"][:2]
    lattice_pts = [gamma_cart]
    for n1i in range(-2, 3):
        for n2i in range(-2, 3):
            if n1i == 0 and n2i == 0:
                continue
            lattice_pts.append(gamma_cart + n1i * span1_xy + n2i * span2_xy)
    vor = Voronoi(np.array(lattice_pts))
    region = vor.regions[vor.point_region[0]]
    if -1 in region or len(region) == 0:
        return None
    poly = vor.vertices[region]
    return np.vstack([poly, poly[0]])


# -------------------------------------------------------- core computation
def compute_all(data, args):
    ef = data["ef"]
    n1, n2 = data["n1"], data["n2"]
    origin, span1, span2 = data["origin"], data["span1"], data["span2"]
    span1_xy, span2_xy = span1[:2], span2[:2]

    crossing = bands_crossing_ef(data["bands"], ef)

    have_calibration = args.a1 is not None and args.a2 is not None and args.alat_bohr is not None
    scale = physical_scale_factor(args.a1, args.a2, args.alat_bohr, span1_xy) if have_calibration else None

    f1 = np.linspace(0, 1, n1, endpoint=False)
    f2 = np.linspace(0, 1, n2, endpoint=False)
    df1, df2 = 1.0 / n1, 1.0 / n2
    S = np.array([span1_xy, span2_xy])
    Sinv = np.linalg.inv(S)

    pts, gamma_cart = high_symmetry_cart(data)

    TILE = 3
    F1s = f1[None, :] + np.arange(-1, TILE - 1)[:, None]
    F2s = f2[None, :] + np.arange(-1, TILE - 1)[:, None]
    F1_full, F2_full = F1s.reshape(-1), F2s.reshape(-1)
    F1G, F2G = np.meshgrid(F1_full, F2_full, indexing="ij")
    KX = origin[0] + F1G * span1_xy[0] + F2G * span2_xy[0]
    KY = origin[1] + F1G * span1_xy[1] + F2G * span2_xy[1]

    Gs = []
    for n1i in (-1, 0, 1):
        for n2i in (-1, 0, 1):
            if n1i == 0 and n2i == 0:
                continue
            Gs.append(n1i * span1_xy + n2i * span2_xy)
    dvec_x, dvec_y = KX - gamma_cart[0], KY - gamma_cart[1]
    d0 = np.hypot(dvec_x, dvec_y)
    ws_mask = np.ones_like(d0, dtype=bool)
    for Gv in Gs:
        dG = np.hypot(dvec_x - Gv[0], dvec_y - Gv[1])
        ws_mask &= d0 <= dG + 1e-9

    band_seg_v = {}  # band_idx -> list of (seg_xy (N,2), seg_v (N,))
    for band_idx in crossing:
        E_single = data["bands"][band_idx][:, :, args.kz_index] - ef

        dEdf1 = (np.roll(E_single, -1, axis=0) - np.roll(E_single, 1, axis=0)) / (2 * df1)
        dEdf2 = (np.roll(E_single, -1, axis=1) - np.roll(E_single, 1, axis=1)) / (2 * df2)
        dEdkx = Sinv[0, 0] * dEdf1 + Sinv[0, 1] * dEdf2
        dEdky = Sinv[1, 0] * dEdf1 + Sinv[1, 1] * dEdf2

        if have_calibration:
            v_x = (dEdkx / scale / HBAR_EVS) * ANG_TO_M
            v_y = (dEdky / scale / HBAR_EVS) * ANG_TO_M
            vmag = np.hypot(v_x, v_y) / 1e6
        else:
            vmag = np.hypot(dEdkx, dEdky)

        E_tiled = np.tile(E_single, (TILE, TILE))
        E_masked = np.where(ws_mask, E_tiled, np.nan)

        cs_fig, cs_ax = plt.subplots()
        cs = cs_ax.contour(KX, KY, E_masked, levels=[0])
        segs = cs.allsegs[0]
        plt.close(cs_fig)

        seg_list = []
        for seg in segs:
            if len(seg) < 2:
                continue
            f_local = Sinv @ (seg - origin[:2]).T
            i_idx = np.mod(np.round(f_local[0] * n1).astype(int), n1)
            j_idx = np.mod(np.round(f_local[1] * n2).astype(int), n2)
            seg_v = vmag[i_idx, j_idx]
            seg_list.append((seg, seg_v))
        band_seg_v[band_idx] = seg_list

    hexagon = ws_hexagon(data, gamma_cart)

    return dict(crossing=crossing, band_seg_v=band_seg_v, pts=pts,
                gamma_cart=gamma_cart, hexagon=hexagon,
                have_calibration=have_calibration)


# -------------------------------------------------------------- PNG output
def write_png(path, results, vmax):
    fig, ax = plt.subplots(figsize=(6, 6))

    for band_idx, seg_list in results["band_seg_v"].items():
        for seg, seg_v in seg_list:
            points = seg.reshape(-1, 1, 2)
            segments = np.concatenate([points[:-1], points[1:]], axis=1)
            lc = LineCollection(segments, cmap="turbo", norm=plt.Normalize(0, vmax))
            lc.set_array((seg_v[:-1] + seg_v[1:]) / 2)
            lc.set_linewidth(1.6)
            ax.add_collection(lc)

    if results["hexagon"] is not None:
        h = results["hexagon"]
        ax.plot(h[:, 0], h[:, 1], color="0.3", lw=1.0, zorder=1)

    path_order = ["Gamma", "M", "K", "Gamma"]
    path_xy = np.array([results["pts"][p] for p in path_order])
    ax.plot(path_xy[:, 0], path_xy[:, 1], color="green", lw=1.0, zorder=2)

    label_map = {"Gamma": r"$\Gamma$", "M": "M", "K": "K"}
    gc = results["gamma_cart"]
    for name, txt in label_map.items():
        px, py = results["pts"][name]
        # push label slightly outward from Gamma so it doesn't sit on the line
        dx, dy = px - gc[0], py - gc[1]
        norm = np.hypot(dx, dy) or 1.0
        offset = 0.04 * max(np.linalg.norm(results["pts"]["M"] - gc), 1e-9) / (norm / max(norm, 1e-9))
        lx, ly = px + 0.15 * dx / norm, py + 0.15 * dy / norm
        ax.annotate(txt, (px, py), xytext=(lx, ly), color="green",
                    fontsize=13, fontweight="bold", ha="center", va="center")

    sm = plt.cm.ScalarMappable(cmap="turbo", norm=plt.Normalize(0, vmax))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", pad=0.08, shrink=0.8)
    cbar.set_label(r"$v$ ($10^6$ m/s)" if results["have_calibration"] else "relative $|v_F|$")

    gc = results["gamma_cart"]
    span_est = np.linalg.norm(results["pts"]["M"] - gc) * 2.5
    ax.set_xlim(gc[0] - span_est, gc[0] + span_est)
    ax.set_ylim(gc[1] - span_est, gc[1] + span_est)
    ax.set_aspect("equal")
    ax.axis("off")

    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close(fig)


# -------------------------------------------------------------- DAT output
def band_avg_velocities(results):
    """Average |v_F| per band, weighted by number of contour points --
    used to assign ONE representative color per band (clean look),
    instead of speckled per-segment coloring."""
    avg = {}
    for band_idx, seg_list in results["band_seg_v"].items():
        all_v = np.concatenate([sv for _, sv in seg_list]) if seg_list else np.array([0.0])
        avg[band_idx] = float(np.mean(all_v))
    return avg


def write_dat(outdir, results):
    """One continuous polyline per contour piece (no choppy sub-segments),
    one file per band. Each PIECE within a band's file is its own xmgrace
    set (so disconnected loops don't get joined by a spurious line), but
    all pieces of a band share the band's single representative color --
    see color_bands.txt for which color to assign."""
    for band_idx, seg_list in results["band_seg_v"].items():
        outfile = os.path.join(outdir, f"band{band_idx}_fs.dat")
        with open(outfile, "w") as f:
            f.write(f"# Fermi surface, band {band_idx} -- one set per "
                     f"disconnected contour piece, assign ONE color to all "
                     f"sets in this file (see color_bands.txt)\n")
            if not seg_list:
                f.write("# (no contour found for this band)\n")
                continue
            for k, (seg, seg_v) in enumerate(seg_list):
                f.write(f"# band{band_idx}_piece{k}\n")
                for row in seg:
                    f.write(f"{row[0]:.6f} {row[1]:.6f}\n")
                f.write("&\n")

    path_order = ["Gamma", "M", "K", "Gamma"]
    path_xy = np.array([results["pts"][p] for p in path_order])
    with open(os.path.join(outdir, "gmk_path.dat"), "w") as f:
        f.write("# Gamma-M-K-Gamma path\n")
        for row in path_xy:
            f.write(f"{row[0]:.6f} {row[1]:.6f}\n")
        f.write("&\n")


    if results["hexagon"] is not None:
        with open(os.path.join(outdir, "ws_hexagon.dat"), "w") as f:
            f.write("# Wigner-Seitz hexagon boundary\n")
            for row in results["hexagon"]:
                f.write(f"{row[0]:.6f} {row[1]:.6f}\n")
            f.write("&\n")

    with open(os.path.join(outdir, "high_sym_points.dat"), "w") as f:
        f.write("# label   kx           ky\n")
        for name in ["Gamma", "M", "K"]:
            x, y = results["pts"][name]
            f.write(f"{name:6s} {x:.6f} {y:.6f}\n")


# -------------------------------------------------------------- AGR output
def write_agr(path, results, colors_by_band):
    """Write a ready-to-open Grace project file: ONE solid color per band
    (continuous polylines, no choppy sub-segments), plus the WS hexagon,
    the Gamma-M-K-Gamma path, and text labels for Gamma/M/K."""
    sets = []       # (color_key, xy array)
    color_start = 2  # grace reserves 0 (white), 1 (black)

    band_ids = list(results["band_seg_v"].keys())
    for bi, band_idx in enumerate(band_ids):
        for seg, seg_v in results["band_seg_v"][band_idx]:
            sets.append((bi, seg))  # bi = index into colors_by_band

    if results["hexagon"] is not None:
        sets.append((None, results["hexagon"]))  # black

    path_order = ["Gamma", "M", "K", "Gamma"]
    path_xy = np.array([results["pts"][p] for p in path_order])
    sets.append(("green", path_xy))

    with open(path, "w") as f:
        f.write("@version 50122\n")
        f.write("@page size 800, 800\n")
        f.write("@default linewidth 1.8\n")

        for i, rgb in enumerate(colors_by_band):
            r, g, b_ = [int(round(c * 255)) for c in rgb[:3]]
            f.write(f'@map color {color_start + i} to ({r}, {g}, {b_}), "band{i}color"\n')
        green_idx = color_start + len(colors_by_band)
        f.write(f'@map color {green_idx} to (0, 150, 0), "pathgreen"\n')

        f.write("@with g0\n")
        f.write("@g0 on\n")
        f.write("@g0 hidden false\n")
        f.write("@g0 type XY\n")

        for idx, (color_key, arr) in enumerate(sets):
            if color_key == "green":
                f.write(f"@    s{idx} line color {green_idx}\n")
                f.write(f"@    s{idx} line linewidth 1.2\n")
            elif color_key is None:
                f.write(f"@    s{idx} line color 1\n")
                f.write(f"@    s{idx} line linewidth 1.0\n")
            else:
                f.write(f"@    s{idx} line color {color_start + color_key}\n")
                f.write(f"@    s{idx} line linewidth 1.8\n")
            f.write(f"@    s{idx} symbol 0\n")
            f.write(f"@    s{idx} line type 1\n")

        # Gamma/M/K text labels as Grace string objects
        for name, txt in [("Gamma", "\\xG\\f{}"), ("M", "M"), ("K", "K")]:
            x, y = results["pts"][name]
            f.write(f'@with string\n@    string on\n@    string loctype world\n')
            f.write(f"@    string {x:.6f}, {y:.6f}\n")
            f.write(f'@    string color 3\n@    string char size 1.2\n')
            f.write(f'@    string def "{txt}"\n')

        for idx, (color_key, arr) in enumerate(sets):
            f.write(f"@target g0.s{idx}\n@type xy\n")
            for row in arr:
                f.write(f"{row[0]:.6f} {row[1]:.6f}\n")
            f.write("&\n")


# --------------------------------------------------------------------- main
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bxsf")
    p.add_argument("--kz-index", type=int, default=0)
    p.add_argument("--a1", nargs=2, type=float, default=None, metavar=("A1X", "A1Y"))
    p.add_argument("--a2", nargs=2, type=float, default=None, metavar=("A2X", "A2Y"))
    p.add_argument("--alat-bohr", type=float, default=None)
    p.add_argument("--vmax", type=float, default=None)
    p.add_argument("--outdir", default="fs_out")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    data = parse_bxsf(args.bxsf)
    print(f"Fermi energy: {data['ef']:.4f} eV")

    results = compute_all(data, args)
    print(f"Bands crossing E_F: {results['crossing']}")

    all_v = np.concatenate([sv for seg_list in results["band_seg_v"].values()
                             for _, sv in seg_list]) if results["band_seg_v"] else np.array([0, 1])
    vmax = args.vmax if args.vmax is not None else np.percentile(all_v, 98)

    avg_v = band_avg_velocities(results)
    band_ids = list(results["band_seg_v"].keys())
    cmap = matplotlib.colormaps.get_cmap("turbo")
    colors_by_band = [cmap(min(avg_v[b] / vmax, 1.0)) for b in band_ids]

    png_path = os.path.join(args.outdir, "fermi_surface.png")
    write_png(png_path, results, vmax)
    print(f"Wrote: {png_path}")

    write_dat(args.outdir, results)
    print(f"Wrote .dat files in: {args.outdir}")

    color_file = os.path.join(args.outdir, "color_bands.txt")
    with open(color_file, "w") as f:
        f.write("# band   avg_v   R      G      B\n")
        for b, rgb in zip(band_ids, colors_by_band):
            r, g, bl, _ = rgb
            f.write(f"{b:5d}  {avg_v[b]:.4f}  {r:.3f}  {g:.3f}  {bl:.3f}\n")
    print(f"Wrote: {color_file}  (assign these colors to each band<N>_fs.dat's sets)")

    agr_path = os.path.join(args.outdir, "fermi_surface.agr")
    write_agr(agr_path, results, colors_by_band)
    print(f"Wrote: {agr_path}  (open directly with: xmgrace {agr_path})")


if __name__ == "__main__":
    main()
