#!/usr/bin/env python3
"""01 Walker circulation edges + weighted center: compute and plot."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.paths import EXP_IDS, RESULT_FIG, RESULT_MAT, WALKER_DIR, ensure_result_dirs
from lib.preprocess import lon_to_360
from lib.style import COLORS, apply_journal_style, dual_axis_errorbar, precession_xticks, save_figure


def find_crossing(lon, psi, win, kind: str, mode: str = "first"):
    lon = np.asarray(lon, float)
    psi = np.asarray(psi, float)
    mask = (lon >= win[0]) & (lon <= win[1])
    lon_w = lon[mask]
    psi_w = psi[mask]
    n = lon_w.size
    if n < 2:
        return np.nan
    idx_list = range(n - 1) if mode == "first" else range(n - 2, -1, -1)
    lon0 = np.nan
    for t in idx_list:
        y1, y2 = psi_w[t], psi_w[t + 1]
        if kind == "neg2pos" and y1 <= 0 and y2 > 0:
            return float(lon_w[t] - y1 * (lon_w[t + 1] - lon_w[t]) / (y2 - y1))
        if kind == "pos2neg" and y1 > 0 and y2 <= 0:
            lon0 = float(lon_w[t] - y1 * (lon_w[t + 1] - lon_w[t]) / (y2 - y1))
            if mode == "first":
                return lon0
    return lon0


def weighted_center_lon(lon, psi, win):
    mask = (lon >= win[0]) & (lon <= win[1]) & (psi > 0)
    if not np.any(mask):
        return np.nan
    w = psi[mask]
    return float(np.sum(lon[mask] * w) / np.sum(w))


def main():
    ensure_result_dirs()
    plev_range = (200.0, 500.0)
    west_window = (100.0, 140.0)
    east_window = (220.0, 275.0)
    center_window = (100.0, 260.0)

    lon = sio.loadmat(WALKER_DIR / "lon.mat")["lon"].ravel().astype(float)
    plev = sio.loadmat(WALKER_DIR / "plev.mat")["plev"].ravel().astype(float)
    psi24 = sio.loadmat(WALKER_DIR / "psi24.mat")["psi24"].astype(float)  # lon x lev x exp

    # Ocean/atmosphere note: lon.mat is -180..180; convert to 0..360 then sort
    lon360 = lon_to_360(lon)
    order = np.argsort(lon360)
    lon360 = lon360[order]
    psi24 = psi24[order, :, :]

    lev_idx = np.where((plev >= plev_range[0]) & (plev <= plev_range[1]))[0]
    if lev_idx.size == 0:
        raise RuntimeError("No levels in plev range")

    n_exp = len(EXP_IDS)
    west = np.full(n_exp, np.nan)
    east = np.full(n_exp, np.nan)
    center = np.full(n_exp, np.nan)

    for k in range(n_exp):
        psi_z = np.nanmean(psi24[:, lev_idx, k], axis=1)
        west[k] = find_crossing(lon360, psi_z, west_window, "neg2pos", "first")
        east[k] = find_crossing(lon360, psi_z, east_window, "pos2neg", "last")
        center[k] = weighted_center_lon(lon360, psi_z, center_window)

    exp = np.array(EXP_IDS, dtype=float)
    mat_path = RESULT_MAT / "walker_edges_precession_results.mat"
    sio.savemat(
        mat_path,
        {
            "expIDs": exp.reshape(-1, 1),
            "westEdgeLon": west.reshape(-1, 1),
            "eastEdgeLon": east.reshape(-1, 1),
            "walk_center_lon": center.reshape(-1, 1),
            "plevRange": np.array(plev_range),
            "westWindow": np.array(west_window),
            "eastWindow": np.array(east_window),
            "centerWindow": np.array(center_window),
            "lon360": lon360.reshape(-1, 1),
        },
    )
    print(f"Saved {mat_path}")

    # Figure: three edges vs precession
    apply_journal_style()
    fig, ax = plt.subplots(figsize=(20 / 2.54, 11 / 2.54), constrained_layout=True)
    ax.plot(exp, west, "o-", color=COLORS["west"], lw=2, ms=7, mfc=COLORS["west"], mec="w", label="West edge (ascending)")
    ax.plot(exp, east, "s-", color=COLORS["east"], lw=2, ms=7, mfc=COLORS["east"], mec="w", label="East edge (descending)")
    ax.plot(exp, center, "^-", color=COLORS["center"], lw=2, ms=7, mfc=COLORS["center"], mec="w", label="Weighted center")
    ax.set_xlabel(r"Precession phase ($^\circ$)", fontweight="bold")
    ax.set_ylabel(r"Longitude ($^\circ$, 0–360)", fontweight="bold")
    precession_xticks(ax, EXP_IDS)
    ax.set_ylim(100, 280)
    ax.legend(frameon=False, loc="best")
    save_figure(fig, "01_walker_edges_vs_precession", RESULT_FIG)

    # Dual axis with Nino3
    from lib.paths import NINO3_DIR

    nino = sio.loadmat(NINO3_DIR / "nino3_enso_chi2_results.mat")
    enso = nino["enso_amp"].ravel()
    ci = nino["chi2_ci_std"]
    err_lo = enso - ci[:, 0]
    err_hi = ci[:, 1] - enso
    fig, _, _ = dual_axis_errorbar(
        EXP_IDS,
        enso,
        err_lo,
        err_hi,
        center,
        np.zeros_like(center),
        np.zeros_like(center),
        color_left=COLORS["nino3"],
        color_right=COLORS["walker"],
        ylabel_left="Nino3 amplitude (K)",
        ylabel_right=r"Walker center longitude ($^\circ$E)",
        label_left="Nino3 amplitude",
        label_right="Walker center longitude",
        mean_left=float(np.nanmean(enso)),
        mean_right=float(np.nanmean(center)),
        right_has_err=False,
    )
    save_figure(fig, "01_nino3_walker_dual_axis", RESULT_FIG)
    print("Walker done.")


if __name__ == "__main__":
    main()
