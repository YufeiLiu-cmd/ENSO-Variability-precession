#!/usr/bin/env python3
"""08 WWV proxy: area-mean equatorial Z20 over 140E–150W (Scheme 1 + dual-axis).

WWV proxy = lon-mean of precomputed 5°S–5°N Z20 in 140°E–150°W.
Same box as USTRW / stratification. Not the full basin Meinen–McPhaden WWV
(120°E–80°W), but a consistent regional warm-water / thermocline metric.

Plot style matches 01–05 (06_combined_nino3.plot_dual).
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from netCDF4 import Dataset

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.paths import (
    CI_LEVEL,
    EXP_IDS,
    NINO3_DIR,
    RESULT_FIG,
    RESULT_MAT,
    Z20_EQ_NC,
    ensure_result_dirs,
)
from lib.preprocess import scheme1_mean_and_ci
from lib.style import COLORS

FIGSIZE_DUAL = (12.5 / 2.54, 4.0 / 2.54)
LAT_NOTE = (-5.0, 5.0)  # already built into z20_eq file
LON_RANGE = (140.0, 210.0)  # 140E–150W
N_TIME = 1200


def apply_style():
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans", "sans-serif"],
            "font.size": 9,
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "pdf.fonttype": 42,
            "savefig.dpi": 300,
        }
    )


def _phase_ticklabels(every: int = 1):
    labs = []
    for i, v in enumerate(EXP_IDS):
        labs.append(f"{v:03d}\u00b0" if i % every == 0 else "")
    return labs


def precession_xticks(ax, show_labels: bool = True, every: int = 2):
    ax.set_xlim(-5, 350)
    ax.set_xticks(EXP_IDS)
    if show_labels:
        ax.set_xticklabels(_phase_ticklabels(every), fontsize=7, rotation=0, ha="center")
    else:
        ax.set_xticklabels([])


def plot_dual(enso, n_lo, n_hi, d):
    """Dual-axis panel styled identically to 01–05."""
    apply_style()
    fig, ax = plt.subplots(figsize=FIGSIZE_DUAL)
    fig.subplots_adjust(left=0.14, right=0.86, bottom=0.18, top=0.92)
    col_n, col_r = COLORS["nino3"], d["color"]
    fs_lab, fs_tick, fs_leg = 8, 6.5, 6.5

    ax.set_ylabel("Nino3 amp. (K)", color=tuple(col_n), fontsize=fs_lab, fontweight="bold")
    ax.tick_params(axis="y", colors=tuple(col_n), labelsize=fs_tick, direction="in", length=3)
    h1 = ax.errorbar(
        EXP_IDS, enso, yerr=np.vstack([n_lo, n_hi]),
        fmt="-o", color=col_n, lw=1.3, ms=3.5, mfc=col_n, mec="w", mew=0.35,
        capsize=2.0, elinewidth=0.9, label="Nino3 amplitude", zorder=3,
    )
    h1m = ax.axhline(np.nanmean(enso), ls="--", color=col_n, lw=0.9, alpha=0.85)
    y0 = np.nanmin(enso - n_lo)
    y1 = np.nanmax(enso + n_hi)
    pad = max((y1 - y0) * 0.16, 0.05)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.spines["left"].set_color(tuple(col_n))
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["left"].set_bounds(y0 - pad, y1 + pad)

    ax2 = ax.twinx()
    ax2.set_ylabel(d["ylabel"], color=tuple(col_r), fontsize=fs_lab, fontweight="bold")
    ax2.tick_params(axis="y", colors=tuple(col_r), labelsize=fs_tick, direction="in", length=3)
    h2 = ax2.errorbar(
        EXP_IDS, d["y"], yerr=np.vstack([d["elo"], d["ehi"]]),
        fmt="-s", color=col_r, lw=1.3, ms=3.5, mfc=col_r, mec="w", mew=0.35,
        capsize=2.0, elinewidth=0.9, label=d["label"], zorder=3,
    )
    h2m = ax2.axhline(np.nanmean(d["y"]), ls="--", color=col_r, lw=0.9, alpha=0.85)
    r0 = np.nanmin(d["y"] - d["elo"])
    r1 = np.nanmax(d["y"] + d["ehi"])
    rpad = max((r1 - r0) * 0.16, 1e-6)
    ax2.set_ylim(r0 - rpad, r1 + rpad)
    ax2.spines["right"].set_color(tuple(col_r))
    ax2.spines["right"].set_linewidth(1.2)
    ax2.spines["right"].set_bounds(r0 - rpad, r1 + rpad)
    ax2.spines["left"].set_visible(False)

    ax.set_xlabel("")
    precession_xticks(ax, show_labels=True, every=2)
    ax.tick_params(axis="x", labelsize=fs_tick, direction="in", length=3)
    for spine in ("top", "bottom"):
        ax.spines[spine].set_linewidth(1.0)
        ax2.spines[spine].set_linewidth(1.0)

    handles = [h1, h1m, h2, h2m]
    labels = [
        "Nino3 amplitude",
        f"mean = {np.nanmean(enso):.2f} K",
        d["label"],
        f"mean = {np.nanmean(d['y']):.3g}",
    ]
    ax.legend(
        handles, labels,
        loc="upper left", frameon=False, fontsize=fs_leg,
        handlelength=1.2, labelspacing=0.2, borderaxespad=0.2,
        handletextpad=0.35, ncol=1,
    )

    out = RESULT_FIG / f"{d['stem']}.png"
    fig.set_size_inches(*FIGSIZE_DUAL)
    fig.savefig(out, dpi=300, facecolor="white")
    plt.close(fig)
    print("saved", out.name)


def load_wwv_timeseries() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return exp_ids, wwv_ts (time×exp), lon_used."""
    if not Z20_EQ_NC.is_file():
        raise FileNotFoundError(Z20_EQ_NC)
    with Dataset(Z20_EQ_NC, "r") as ds:
        lon = np.asarray(ds.variables["lon"][:], float)
        exp_id = np.asarray(ds.variables["exp_id"][:], float)
        z20 = np.asarray(ds.variables["z20"][:, :N_TIME, :], float)  # lon, time, exp

    lon_mask = (lon >= LON_RANGE[0]) & (lon <= LON_RANGE[1])
    if not np.any(lon_mask):
        raise RuntimeError(f"No lon in {LON_RANGE}")
    # Already 5S–5N mean; lon-mean over the WWV box
    wwv = np.nanmean(z20[lon_mask, :, :], axis=0)  # time × exp
    return exp_id, wwv, lon[lon_mask]


def main():
    ensure_result_dirs()
    exp_id, wwv_ts, lon_used = load_wwv_timeseries()
    # Align columns to EXP_IDS order
    n_exp = len(EXP_IDS)
    ts_all = np.full((N_TIME, n_exp), np.nan)
    means = np.full(n_exp, np.nan)
    ci = np.full((n_exp, 2), np.nan)
    n_samp = np.full(n_exp, np.nan)

    for k, eid in enumerate(EXP_IDS):
        j = int(np.where(np.isclose(exp_id, eid))[0][0])
        ts = wwv_ts[:, j]
        ts_all[:, k] = ts
        mu, lo, hi, diag = scheme1_mean_and_ci(ts, CI_LEVEL)
        means[k] = mu
        ci[k] = (lo, hi)
        n_samp[k] = diag["n"]
        print(
            f"EXP {eid:03d} WWV(Z20)={mu:.2f} m CI=[{lo:.2f},{hi:.2f}] n={diag['n']}",
            flush=True,
        )

    mat_path = RESULT_MAT / "wwv_precession_results.mat"
    sio.savemat(
        mat_path,
        {
            "expIDs": np.array(EXP_IDS, float).reshape(-1, 1),
            "wwv_ts": ts_all,
            "wwv_mean": means.reshape(-1, 1),
            "wwv_ci": ci,
            "nSamples": n_samp.reshape(-1, 1),
            "ciLevel": CI_LEVEL,
            "latRange": np.array(LAT_NOTE),
            "lonRange": np.array(LON_RANGE),
            "lonUsed": lon_used,
            "z20_file": str(Z20_EQ_NC),
            "preprocess": (
                "WWV proxy = lon-mean of 5S-5N Z20 over 140E-150W; "
                "Scheme1 mean (3mo runmean) + annual-mean t-CI (n~100)"
            ),
        },
    )
    print(f"Saved {mat_path}")

    nino = sio.loadmat(NINO3_DIR / "nino3_enso_chi2_results.mat")
    enso = nino["enso_amp"].ravel().astype(float)
    nci = nino["chi2_ci_std"]
    n_lo = enso - nci[:, 0]
    n_hi = nci[:, 1] - enso
    err_lo = means - ci[:, 0]
    err_hi = ci[:, 1] - means

    plot_dual(
        enso, n_lo, n_hi,
        dict(
            key="wwv",
            stem="08_nino3_wwv_dual_axis",
            y=means,
            elo=err_lo,
            ehi=err_hi,
            color=COLORS["wwv"],
            ylabel="WWV (m)",
            label="WWV",
        ),
    )
    print("WWV done.")
    print(f"r(Niño3, WWV) = {np.corrcoef(enso, means)[0, 1]:+.3f}")


if __name__ == "__main__":
    main()
