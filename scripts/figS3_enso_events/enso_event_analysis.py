#!/usr/bin/env python3
"""
ENSO event statistics for 24 precession-phase experiments (v3).

Pipeline
--------
1. Read XXXtsw.nc (000 … 345, step 15°)
2. Niño3 area-mean SST (5°S–5°N, 210–270°E)
3. Remove seasonal cycle → SSTA; remove linear trend
4. 3-month centered running mean → Niño3 index I(t)
5. Events on I(t):
      El Niño : I >  +0.5 K for ≥ 5 months
      La Niña : I <  −0.5 K for ≥ 5 months
6. Figures (English, Nature-style); all CIs are 95%
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
from scipy import stats
from scipy.io import savemat

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
DATA_DIR = Path(
    __import__("os").environ.get(
        "ENSO_DATA_ROOT",
        str(HERE.parents[2] / "data"),
    )
) / "tsw"
EXP_IDS = np.arange(0, 346, 15, dtype=int)
N_TIME = 1200
N_YEARS = N_TIME / 12.0

LAT_MIN, LAT_MAX = -5.0, 5.0
LON_MIN, LON_MAX = 210.0, 270.0

THRESH = 0.5
MIN_DUR = 5
SMOOTH_WIN = 3
CI_LEVEL = 0.95

OUT_MAT = HERE / "enso_event_stats.mat"
OUT_CSV = HERE / "enso_event_stats.csv"
OUT_NINO3 = HERE / "nino3_index_24exp.npy"
OUT_NINO3_RAW = HERE / "nino3_index_raw_24exp.npy"

COL_EN = "#C0392B"
COL_LN = "#2471A3"
COL_LINE = "#2C3E50"
COL_GRID = "#E8E8E8"

CRITERIA_TEXT = (
    f"Criteria: 3-mo running-mean Niño3\n"
    f"|SSTA| > {THRESH} K for ≥ {MIN_DUR} months"
)


def log(msg: str) -> None:
    print(msg, flush=True)


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 13,
            "axes.labelsize": 15,
            "axes.titlesize": 15,
            "xtick.labelsize": 11,
            "ytick.labelsize": 12,
            "legend.fontsize": 11.5,
            "axes.linewidth": 1.35,
            "xtick.major.width": 1.1,
            "ytick.major.width": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


# ---------------------------------------------------------------------------
# Data processing
# ---------------------------------------------------------------------------
def read_nino3_sst(fpath: Path) -> np.ndarray:
    with Dataset(fpath, "r") as nc:
        lon = np.asarray(nc.variables["lon"][:], dtype=np.float64)
        lat = np.asarray(nc.variables["lat"][:], dtype=np.float64)
        tsw = np.asarray(nc.variables["tsw"][:N_TIME], dtype=np.float64)
    lat_idx = np.where((lat >= LAT_MIN) & (lat <= LAT_MAX))[0]
    lon_idx = np.where((lon >= LON_MIN) & (lon <= LON_MAX))[0]
    sub = tsw[:, lat_idx][:, :, lon_idx]
    w_lat = np.cos(np.deg2rad(lat[lat_idx]))
    w_lat /= w_lat.sum()
    return np.einsum("tjl,j->tl", sub, w_lat).mean(axis=1)


def remove_seasonal_cycle(sst: np.ndarray) -> np.ndarray:
    months = np.arange(sst.size) % 12
    clim = np.array([sst[months == m].mean() for m in range(12)])
    return sst - clim[months]


def detrend_linear(x: np.ndarray) -> np.ndarray:
    t = np.arange(x.size, dtype=np.float64)
    return x - np.polyval(np.polyfit(t, x, 1), t)


def running_mean(x: np.ndarray, win: int = SMOOTH_WIN) -> np.ndarray:
    pad = win // 2
    xp = np.pad(x, (pad, pad), mode="edge")
    return np.convolve(xp, np.ones(win) / win, mode="valid")[: x.size]


def compute_nino3_index(sst: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    raw = detrend_linear(remove_seasonal_cycle(sst))
    return raw, running_mean(raw, SMOOTH_WIN)


def find_events(index: np.ndarray, sign: int, thresh: float, min_dur: int):
    mask = (sign * index) > thresh
    events = []
    i, n = 0, mask.size
    while i < n:
        if not mask[i]:
            i += 1
            continue
        j = i
        while j < n and mask[j]:
            j += 1
        dur = j - i
        if dur >= min_dur:
            seg = index[i:j]
            peak_rel = int(np.argmax(seg) if sign > 0 else np.argmin(seg))
            events.append(
                {
                    "start": i,
                    "end": j - 1,
                    "duration": dur,
                    "peak": float(seg[peak_rel]),
                    "peak_month": i + peak_rel,
                }
            )
        i = j
    return events


def mean_ci_95(arr: np.ndarray) -> tuple[float, float, float]:
    """Return (mean, err_low, err_high) for 95% t-based CI; err = distance from mean."""
    arr = np.asarray(arr, dtype=float)
    n = arr.size
    if n == 0:
        return np.nan, np.nan, np.nan
    m = float(arr.mean())
    if n == 1:
        return m, 0.0, 0.0
    sem = float(arr.std(ddof=1) / np.sqrt(n))
    tcrit = float(stats.t.ppf(0.5 + CI_LEVEL / 2.0, n - 1))
    half = tcrit * sem
    return m, half, half


def chi2_std_ci_95(sigma: float, n: int) -> tuple[float, float]:
    """
    95% CI for population std from sample std (chi-square).
    Returns (ci_low, ci_high) absolute bounds.
    """
    if n < 2 or not np.isfinite(sigma) or sigma <= 0:
        return np.nan, np.nan
    nu = n - 1
    # P(χ²_{1-α/2} < ν s²/σ² < χ²_{α/2})  → invert for σ
    chi_lo = float(stats.chi2.ppf(0.5 - CI_LEVEL / 2.0, nu))  # 0.025
    chi_hi = float(stats.chi2.ppf(0.5 + CI_LEVEL / 2.0, nu))  # 0.975
    # σ_low = s * sqrt(ν / χ²_hi),  σ_high = s * sqrt(ν / χ²_lo)
    return sigma * np.sqrt(nu / chi_hi), sigma * np.sqrt(nu / chi_lo)


def summarize_events(events: list[dict]) -> dict:
    empty = {
        "count": 0,
        "mean_duration": np.nan,
        "ci_dur_lo": np.nan,
        "ci_dur_hi": np.nan,
        "mean_peak": np.nan,
        "ci_peak_lo": np.nan,
        "ci_peak_hi": np.nan,
        "mean_abs_peak": np.nan,
        "ci_abs_peak_lo": np.nan,
        "ci_abs_peak_hi": np.nan,
    }
    if not events:
        return empty
    durs = np.array([e["duration"] for e in events], dtype=float)
    peaks = np.array([e["peak"] for e in events], dtype=float)
    abs_peaks = np.abs(peaks)
    md, dlo, dhi = mean_ci_95(durs)
    mp, plo, phi = mean_ci_95(peaks)
    ma, alo, ahi = mean_ci_95(abs_peaks)
    return {
        "count": len(events),
        "mean_duration": md,
        "ci_dur_lo": dlo,
        "ci_dur_hi": dhi,
        "mean_peak": mp,
        "ci_peak_lo": plo,
        "ci_peak_hi": phi,
        "mean_abs_peak": ma,
        "ci_abs_peak_lo": alo,
        "ci_abs_peak_hi": ahi,
    }


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------
def run_analysis():
    n_exp = EXP_IDS.size
    nino3_raw = np.full((N_TIME, n_exp), np.nan)
    nino3_all = np.full((N_TIME, n_exp), np.nan)

    n_elnino = np.zeros(n_exp, dtype=int)
    n_lanina = np.zeros(n_exp, dtype=int)
    mean_dur_en = np.full(n_exp, np.nan)
    mean_dur_ln = np.full(n_exp, np.nan)
    ci_dur_en = np.full((n_exp, 2), np.nan)
    ci_dur_ln = np.full((n_exp, 2), np.nan)
    mean_peak_en = np.full(n_exp, np.nan)
    mean_abs_peak_ln = np.full(n_exp, np.nan)
    ci_peak_en = np.full((n_exp, 2), np.nan)
    ci_abs_peak_ln = np.full((n_exp, 2), np.nan)
    std_nino3 = np.full(n_exp, np.nan)
    chi2_ci_std = np.full((n_exp, 2), np.nan)

    all_en_events: list[list[dict]] = []
    all_ln_events: list[list[dict]] = []

    log("=" * 64)
    log("ENSO event analysis v3 | 24 precession experiments")
    log(f"  Smooth: {SMOOTH_WIN}-mo mean | thresh +/- {THRESH} K | min {MIN_DUR} mo")
    log(f"  All confidence intervals: {int(CI_LEVEL * 100)}%")
    log("=" * 64)

    t_all = time.time()
    for k, exp_id in enumerate(EXP_IDS):
        fpath = DATA_DIR / f"{exp_id:03d}tsw.nc"
        if not fpath.exists():
            raise FileNotFoundError(fpath)
        t0 = time.time()
        log(f"[{k + 1:02d}/{n_exp}] EXP {exp_id:03d}")
        sst = read_nino3_sst(fpath)
        raw, smooth = compute_nino3_index(sst)
        nino3_raw[:, k] = raw
        nino3_all[:, k] = smooth

        sigma = float(np.std(smooth, ddof=1))
        std_nino3[k] = sigma
        lo, hi = chi2_std_ci_95(sigma, N_TIME)
        chi2_ci_std[k] = [lo, hi]

        en = find_events(smooth, +1, THRESH, MIN_DUR)
        ln = find_events(smooth, -1, THRESH, MIN_DUR)
        all_en_events.append(en)
        all_ln_events.append(ln)
        s_en = summarize_events(en)
        s_ln = summarize_events(ln)

        n_elnino[k] = s_en["count"]
        n_lanina[k] = s_ln["count"]
        mean_dur_en[k] = s_en["mean_duration"]
        mean_dur_ln[k] = s_ln["mean_duration"]
        ci_dur_en[k] = [s_en["ci_dur_lo"], s_en["ci_dur_hi"]]
        ci_dur_ln[k] = [s_ln["ci_dur_lo"], s_ln["ci_dur_hi"]]
        mean_peak_en[k] = s_en["mean_peak"]
        mean_abs_peak_ln[k] = s_ln["mean_abs_peak"]
        ci_peak_en[k] = [s_en["ci_peak_lo"], s_en["ci_peak_hi"]]
        ci_abs_peak_ln[k] = [s_ln["ci_abs_peak_lo"], s_ln["ci_abs_peak_hi"]]

        log(
            f"         EN={s_en['count']:2d}  LN={s_ln['count']:2d}  "
            f"σ={sigma:.3f} [{lo:.3f}, {hi:.3f}]  ({time.time() - t0:.1f} s)"
        )

    log(f"All done in {time.time() - t_all:.1f} s")

    np.save(OUT_NINO3, nino3_all)
    np.save(OUT_NINO3_RAW, nino3_raw)

    savemat(
        OUT_MAT,
        {
            "expIDs": EXP_IDS.astype(np.uint16),
            "nino3_index": nino3_all,
            "nino3_index_raw": nino3_raw,
            "n_elnino": n_elnino,
            "n_lanina": n_lanina,
            "mean_dur_elnino": mean_dur_en,
            "mean_dur_lanina": mean_dur_ln,
            "ci95_dur_elnino": ci_dur_en,
            "ci95_dur_lanina": ci_dur_ln,
            "mean_peak_elnino": mean_peak_en,
            "mean_abs_peak_lanina": mean_abs_peak_ln,
            "ci95_peak_elnino": ci_peak_en,
            "ci95_abs_peak_lanina": ci_abs_peak_ln,
            "std_nino3": std_nino3,
            "chi2_ci95_std": chi2_ci_std,
            "ciLevel": CI_LEVEL,
            "thresh": THRESH,
            "min_duration": MIN_DUR,
            "smooth_window": SMOOTH_WIN,
            "n_years": N_YEARS,
            "n_time": N_TIME,
            "latMin": LAT_MIN,
            "latMax": LAT_MAX,
            "lonMin": LON_MIN,
            "lonMax": LON_MAX,
        },
    )

    with OUT_CSV.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "exp_id",
                "n_elnino",
                "n_lanina",
                "mean_dur_elnino",
                "ci95_dur_en_lo",
                "ci95_dur_en_hi",
                "mean_dur_lanina",
                "ci95_dur_ln_lo",
                "ci95_dur_ln_hi",
                "mean_peak_elnino",
                "ci95_peak_en_lo",
                "ci95_peak_en_hi",
                "mean_abs_peak_lanina",
                "ci95_abs_peak_ln_lo",
                "ci95_abs_peak_ln_hi",
                "std_nino3",
                "chi2_ci95_std_lo",
                "chi2_ci95_std_hi",
            ]
        )
        for k, exp_id in enumerate(EXP_IDS):
            w.writerow(
                [
                    f"{exp_id:03d}",
                    n_elnino[k],
                    n_lanina[k],
                    f"{mean_dur_en[k]:.3f}",
                    f"{ci_dur_en[k, 0]:.3f}",
                    f"{ci_dur_en[k, 1]:.3f}",
                    f"{mean_dur_ln[k]:.3f}",
                    f"{ci_dur_ln[k, 0]:.3f}",
                    f"{ci_dur_ln[k, 1]:.3f}",
                    f"{mean_peak_en[k]:.3f}",
                    f"{ci_peak_en[k, 0]:.3f}",
                    f"{ci_peak_en[k, 1]:.3f}",
                    f"{mean_abs_peak_ln[k]:.3f}",
                    f"{ci_abs_peak_ln[k, 0]:.3f}",
                    f"{ci_abs_peak_ln[k, 1]:.3f}",
                    f"{std_nino3[k]:.4f}",
                    f"{chi2_ci_std[k, 0]:.4f}",
                    f"{chi2_ci_std[k, 1]:.4f}",
                ]
            )

    log(f"Saved: {OUT_MAT.name}, {OUT_CSV.name}")

    return {
        "expIDs": EXP_IDS,
        "nino3_all": nino3_all,
        "n_elnino": n_elnino,
        "n_lanina": n_lanina,
        "mean_dur_en": mean_dur_en,
        "mean_dur_ln": mean_dur_ln,
        "ci_dur_en": ci_dur_en,
        "ci_dur_ln": ci_dur_ln,
        "mean_peak_en": mean_peak_en,
        "mean_abs_peak_ln": mean_abs_peak_ln,
        "ci_peak_en": ci_peak_en,
        "ci_abs_peak_ln": ci_abs_peak_ln,
        "std_nino3": std_nino3,
        "chi2_ci_std": chi2_ci_std,
        "all_en_events": all_en_events,
        "all_ln_events": all_ln_events,
    }


# ---------------------------------------------------------------------------
# Figure helpers
# ---------------------------------------------------------------------------
def _save(fig: plt.Figure, stem: str) -> None:
    for ext in ("png", "pdf"):
        fig.savefig(HERE / f"{stem}.{ext}", dpi=300, facecolor="white")
    plt.close(fig)
    log(f"  → {stem}.png / .pdf")


def _add_criteria(ax, x=0.98, y=0.98, ha="right", va="top", fontsize=9.0):
    ax.text(
        x,
        y,
        CRITERIA_TEXT,
        transform=ax.transAxes,
        fontsize=fontsize,
        ha=ha,
        va=va,
        ma="right",
        color="#333333",
        bbox=dict(
            boxstyle="round,pad=0.32",
            facecolor="white",
            edgecolor="#CCCCCC",
            alpha=0.92,
            linewidth=0.8,
        ),
        zorder=10,
    )


def _panel_label(ax, lab: str) -> None:
    ax.text(
        0.02,
        0.98,
        lab,
        transform=ax.transAxes,
        fontsize=17,
        fontweight="bold",
        va="top",
        ha="left",
        zorder=20,
        bbox=dict(boxstyle="square,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85),
    )


def _style_ax(ax) -> None:
    ax.yaxis.grid(True, color=COL_GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_linewidth(1.35)


def _xticks_all(ax, exp, rotation=0):
    ax.set_xticks(exp)
    ax.set_xticklabels([f"{int(e):03d}" for e in exp], rotation=rotation, ha="center")


def _xticks_30(ax, exp):
    """Label every 30° (000, 030, …, 330); keep all tick marks optional via major only."""
    majors = exp[exp % 30 == 0]
    ax.set_xticks(majors)
    ax.set_xticklabels([f"{int(e):03d}" for e in majors], rotation=0, ha="center")


def _errbar_kwargs(color):
    return dict(
        color=color,
        lw=1.9,
        ms=6.5,
        mfc=color,
        mec="white",
        mew=0.9,
        capsize=3.0,
        elinewidth=1.15,
        ecolor=color,
        zorder=3,
    )


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def _draw_counts_panel(ax, r, *, every_30=False, ylim=40, show_criteria=True, numbers_fs=7.5):
    exp = r["expIDs"]
    n = len(exp)
    x = np.arange(n)
    w = 0.38
    b1 = ax.bar(x - w / 2, r["n_elnino"], w, color=COL_EN, edgecolor="white", lw=0.4, label="El Niño", zorder=3)
    b2 = ax.bar(x + w / 2, r["n_lanina"], w, color=COL_LN, edgecolor="white", lw=0.4, label="La Niña", zorder=3)
    for bar, val in zip(b1, r["n_elnino"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.3,
            f"{int(val)}",
            ha="center",
            va="bottom",
            fontsize=numbers_fs,
            color=COL_EN,
            fontweight="bold",
        )
    for bar, val in zip(b2, r["n_lanina"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.3,
            f"{int(val)}",
            ha="center",
            va="bottom",
            fontsize=numbers_fs,
            color=COL_LN,
            fontweight="bold",
        )
    if every_30:
        majors = np.where(exp % 30 == 0)[0]
        ax.set_xticks(majors)
        ax.set_xticklabels([f"{int(exp[i]):03d}" for i in majors], rotation=0)
    else:
        ax.set_xticks(x)
        ax.set_xticklabels([f"{int(e):03d}" for e in exp], rotation=0)
    ax.set_xlabel("Precession phase (°)")
    ax.set_ylabel("Number of events")
    ax.set_xlim(-0.45, n - 0.55)
    ax.set_ylim(0, ylim)
    _style_ax(ax)
    ax.legend(frameon=False, loc="upper left", fontsize=10, ncol=2, bbox_to_anchor=(0.08, 1.0))
    if show_criteria:
        _add_criteria(ax, x=0.98, y=0.98, ha="right", va="top", fontsize=8.0)


def _draw_sigma_panel(ax, r) -> None:
    exp = r["expIDs"]
    sig = r["std_nino3"]
    ci = r["chi2_ci_std"]
    yerr = np.vstack([sig - ci[:, 0], ci[:, 1] - sig])
    ax.errorbar(
        exp,
        sig,
        yerr=yerr,
        fmt="-D",
        color=COL_LINE,
        lw=1.9,
        ms=6.5,
        mfc=COL_LINE,
        mec="white",
        mew=0.9,
        capsize=3.0,
        elinewidth=1.15,
        ecolor=COL_LINE,
        label=r"Niño3 $\sigma$ (3-mo mean)",
        zorder=3,
    )
    _xticks_30(ax, exp)
    ax.set_xlabel("Precession phase (°)")
    ax.set_ylabel("Niño3 standard deviation (K)")
    ax.set_xlim(-2, 347)
    _style_ax(ax)
    ax.legend(frameon=False, loc="upper left", fontsize=10, bbox_to_anchor=(0.22, 1.0))
    ax.text(
        0.98,
        0.02,
        r"95% CI ($\chi^2$)",
        transform=ax.transAxes,
        fontsize=8.5,
        ha="right",
        va="bottom",
        color="#555",
    )


def _draw_timeseries_panel(ax, r, k: int, *, show_note: bool = False) -> None:
    years = np.arange(N_TIME) / 12.0
    idx = r["nino3_all"][:, k]
    exp = int(r["expIDs"][k])
    ax.fill_between(years, 0, idx, where=idx > 0, color=COL_EN, alpha=0.22, lw=0)
    ax.fill_between(years, 0, idx, where=idx < 0, color=COL_LN, alpha=0.22, lw=0)
    ax.plot(years, idx, color=COL_LINE, lw=0.55, zorder=3)
    ax.axhline(THRESH, color=COL_EN, ls="--", lw=1.0, alpha=0.85)
    ax.axhline(-THRESH, color=COL_LN, ls="--", lw=1.0, alpha=0.85)
    ax.axhline(0, color="k", lw=0.55, alpha=0.45)
    for e in r["all_en_events"][k]:
        ax.axvspan(e["start"] / 12.0, (e["end"] + 1) / 12.0, color=COL_EN, alpha=0.12, lw=0)
    for e in r["all_ln_events"][k]:
        ax.axvspan(e["start"] / 12.0, (e["end"] + 1) / 12.0, color=COL_LN, alpha=0.12, lw=0)
    ax.set_ylabel("Niño3 (°C)")
    ax.set_xlim(0, N_YEARS)
    ylim = max(2.5, np.nanpercentile(np.abs(idx), 99.5) * 1.15)
    ax.set_ylim(-ylim, ylim)
    _style_ax(ax)
    ax.text(
        0.99,
        0.92,
        f"Precession {exp:03d}°   El Niño = {r['n_elnino'][k]}   La Niña = {r['n_lanina'][k]}",
        transform=ax.transAxes,
        fontsize=12,
        ha="right",
        va="top",
    )
    if show_note:
        ax.text(
            0.99,
            0.08,
            "3-month running-mean Niño3",
            transform=ax.transAxes,
            fontsize=9.5,
            ha="right",
            va="bottom",
            color="#555555",
        )


def fig_combined(r: dict) -> None:
    """
    Single combined figure:
      a  event counts (Fig.4a)
      b  Niño3 σ with χ² 95% CI (Fig.4d)
      c  090° Niño3 timeseries (from Fig.6)
      d  270° Niño3 timeseries (from Fig.6)
    """
    # experiment indices: 090° → 6, 270° → 18
    k_090 = int(np.where(r["expIDs"] == 90)[0][0])
    k_270 = int(np.where(r["expIDs"] == 270)[0][0])

    fig = plt.figure(figsize=(13.6, 11.2))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1.0, 1.0], hspace=0.32, wspace=0.28)

    ax_a = fig.add_subplot(gs[0, 0])
    _draw_counts_panel(ax_a, r, every_30=True, ylim=40, show_criteria=True, numbers_fs=7)
    _panel_label(ax_a, "a")

    ax_b = fig.add_subplot(gs[0, 1])
    _draw_sigma_panel(ax_b, r)
    _panel_label(ax_b, "b")

    ax_c = fig.add_subplot(gs[1, :])
    _draw_timeseries_panel(ax_c, r, k_090, show_note=True)
    _panel_label(ax_c, "c")

    ax_d = fig.add_subplot(gs[2, :], sharex=ax_c)
    _draw_timeseries_panel(ax_d, r, k_270)
    _panel_label(ax_d, "d")
    ax_d.set_xlabel("Model year")

    _save(fig, "fig_enso_combined")


def make_all_figures(r: dict) -> None:
    setup_style()
    log("Making combined figure ...")
    fig_combined(r)
    log("Figure saved.")


def main():
    results = run_analysis()
    make_all_figures(results)
    log("Done.")


if __name__ == "__main__":
    main()
