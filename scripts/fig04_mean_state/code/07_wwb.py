#!/usr/bin/env python3
"""07 Westerly wind-burst (WWB) statistics from ustrw: compute + plot (Scheme-1 CI).

Definition (monthly τx proxy of WWBs):
  Region: 5°S–5°N, 140°E–150°W (same equatorial Pacific box as USTRW/MLD/strat)
  Anomaly: remove seasonal cycle + linear detrend
  WWB month: τx' > THRESHOLD (Pa)
  Event: contiguous WWB months (duration ≥ 1 month)
  Plotted statistic: mean annual event frequency (yr⁻¹)
  95% CI: Student-t on annual event counts (n = 100 years)
           — same uncertainty philosophy as Scheme 1 for mean-state vars
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from netCDF4 import Dataset
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lib.paths import (
    CI_LEVEL,
    EXP_IDS,
    NINO3_DIR,
    RESULT_FIG,
    RESULT_MAT,
    USTRW_DIR,
    ensure_result_dirs,
)
from lib.preprocess import linear_detrend, remove_seasonal_cycle
from lib.style import COLORS

# Match 01–05 paper panels from 06_combined_nino3.plot_dual
FIGSIZE_DUAL = (12.5 / 2.54, 4.0 / 2.54)

LAT_RANGE = (-5.0, 5.0)
LON_RANGE = (140.0, 210.0)  # 140E–150W; match USTRW / MLD / strat box
THRESHOLD = 0.015  # Pa; monthly τx' threshold
N_TIME = 1200
N_YEAR = N_TIME // 12


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
        if i % every == 0:
            labs.append(f"{v:03d}\u00b0")
        else:
            labs.append("")
    return labs


def precession_xticks(ax, show_labels: bool = True, every: int = 2):
    ax.set_xlim(-5, 350)
    ax.set_xticks(EXP_IDS)
    if show_labels:
        ax.set_xticklabels(_phase_ticklabels(every), fontsize=7, rotation=0, ha="center")
        for lab in ax.get_xticklabels():
            lab.set_rotation(0)
    else:
        ax.set_xticklabels([])


def plot_dual(enso, n_lo, n_hi, d):
    """Dual-axis panel styled identically to 01–05 (06_combined_nino3.plot_dual)."""
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


def area_mean_ustrw_fast(fpath: Path) -> np.ndarray:
    with Dataset(fpath, "r") as ds:
        lon = np.asarray(ds.variables["lon"][:], dtype=float)
        lat = np.asarray(ds.variables["lat"][:], dtype=float)
        lat_idx = np.where((lat >= LAT_RANGE[0]) & (lat <= LAT_RANGE[1]))[0]
        lon_idx = np.where((lon >= LON_RANGE[0]) & (lon <= LON_RANGE[1]))[0]
        la0, la1 = int(lat_idx.min()), int(lat_idx.max()) + 1
        lo0, lo1 = int(lon_idx.min()), int(lon_idx.max()) + 1
        u = np.asarray(ds.variables["ustrw"][0:N_TIME, la0:la1, lo0:lo1], dtype=float)
        lat_sub = lat[la0:la1]
        lon_sub = lon[lo0:lo1]
        lat_mask = (lat_sub >= LAT_RANGE[0]) & (lat_sub <= LAT_RANGE[1])
        lon_mask = (lon_sub >= LON_RANGE[0]) & (lon_sub <= LON_RANGE[1])
        u = u[:, lat_mask][:, :, lon_mask]
        wlat = np.cos(np.deg2rad(lat_sub[lat_mask]))
        w2d = np.ones((wlat.size, int(lon_mask.sum())), dtype=float) * wlat[:, None]
        w2d /= w2d.sum()
        ts = np.nansum(u * w2d[None, :, :], axis=(1, 2))
    return ts.astype(float)


def identify_events(mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return start/end indices (end exclusive) of True runs."""
    m = np.asarray(mask, dtype=bool).ravel()
    d = np.diff(np.concatenate([[0], m.astype(np.int8), [0]]))
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0]
    return starts, ends


def annual_event_counts(starts: np.ndarray, n_year: int = N_YEAR) -> np.ndarray:
    if starts.size == 0:
        return np.zeros(n_year, dtype=float)
    yr = (starts // 12).astype(int)
    yr = yr[(yr >= 0) & (yr < n_year)]
    return np.bincount(yr, minlength=n_year).astype(float)[:n_year]


def annual_month_counts(mask: np.ndarray, n_year: int = N_YEAR) -> np.ndarray:
    m = np.asarray(mask, dtype=float).ravel()[: n_year * 12]
    return m.reshape(n_year, 12).sum(axis=1)


def mean_and_ci_from_annual(annual: np.ndarray, ci_level: float = CI_LEVEL):
    valid = np.asarray(annual, dtype=float)
    valid = valid[np.isfinite(valid)]
    n = valid.size
    if n < 3:
        raise ValueError(f"Not enough annual samples for CI (n={n})")
    mu = float(np.mean(valid))
    sigma = float(np.std(valid, ddof=1))
    sem = sigma / np.sqrt(n)
    alpha = 1.0 - ci_level
    tcrit = float(stats.t.ppf(1.0 - alpha / 2.0, df=n - 1))
    half = tcrit * sem
    return mu, mu - half, mu + half, {
        "n": n,
        "sigma_annual": sigma,
        "sem": sem,
        "tcrit": tcrit,
        "half_width": half,
        "annual": valid,
    }


def wwb_for_series(ts: np.ndarray, threshold: float = THRESHOLD) -> dict:
    anom = linear_detrend(remove_seasonal_cycle(ts))
    mask = anom > threshold
    starts, ends = identify_events(mask)
    ann_ev = annual_event_counts(starts)
    ann_mo = annual_month_counts(mask)
    freq, lo, hi, diag = mean_and_ci_from_annual(ann_ev)

    peaks = np.array(
        [float(np.nanmax(anom[s:e])) for s, e in zip(starts, ends)], dtype=float
    )
    if peaks.size:
        # annual mean peak intensity (NaN years with no event -> omit in CI)
        peak_year = (starts // 12).astype(int)
        ann_peak = np.full(N_YEAR, np.nan)
        for y in range(N_YEAR):
            sel = peaks[peak_year == y]
            if sel.size:
                ann_peak[y] = float(np.nanmean(sel))
        inten_valid = ann_peak[np.isfinite(ann_peak)]
        if inten_valid.size >= 3:
            inten_mu, inten_lo, inten_hi, inten_diag = mean_and_ci_from_annual(inten_valid)
        else:
            inten_mu = float(np.nanmean(peaks))
            inten_lo = inten_hi = np.nan
            inten_diag = {"n": inten_valid.size}
    else:
        inten_mu = inten_lo = inten_hi = np.nan
        inten_diag = {"n": 0}
        ann_peak = np.full(N_YEAR, np.nan)

    mo_mu, mo_lo, mo_hi, mo_diag = mean_and_ci_from_annual(ann_mo)

    return {
        "anom": anom,
        "mask": mask.astype(float),
        "n_events": int(starts.size),
        "n_months": int(mask.sum()),
        "freq_mean": freq,
        "freq_ci": (lo, hi),
        "freq_annual": ann_ev,
        "freq_diag": diag,
        "month_mean": mo_mu,
        "month_ci": (mo_lo, mo_hi),
        "month_annual": ann_mo,
        "month_diag": mo_diag,
        "inten_mean": inten_mu,
        "inten_ci": (inten_lo, inten_hi),
        "inten_annual": ann_peak,
        "inten_diag": inten_diag,
        "starts": starts,
        "ends": ends,
        "peaks": peaks,
    }


def main():
    ensure_result_dirs()
    if not USTRW_DIR.is_dir():
        raise FileNotFoundError(
            f"USTRW_DIR missing: {USTRW_DIR} (expected symlink to 24 × *.nc)"
        )

    n_exp = len(EXP_IDS)
    ts_all = np.full((N_TIME, n_exp), np.nan)
    anom_all = np.full((N_TIME, n_exp), np.nan)
    freq = np.full(n_exp, np.nan)
    freq_ci = np.full((n_exp, 2), np.nan)
    month_rate = np.full(n_exp, np.nan)
    month_ci = np.full((n_exp, 2), np.nan)
    inten = np.full(n_exp, np.nan)
    inten_ci = np.full((n_exp, 2), np.nan)
    n_events = np.full(n_exp, np.nan)
    n_samp = np.full(n_exp, np.nan)
    freq_annual = np.full((N_YEAR, n_exp), np.nan)

    for k, eid in enumerate(EXP_IDS):
        fpath = USTRW_DIR / f"{eid:03d}.nc"
        print(f"Reading {fpath.name} ...", flush=True)
        ts = area_mean_ustrw_fast(fpath)
        ts_all[:, k] = ts
        out = wwb_for_series(ts, THRESHOLD)
        anom_all[:, k] = out["anom"]
        freq[k] = out["freq_mean"]
        freq_ci[k] = out["freq_ci"]
        month_rate[k] = out["month_mean"]
        month_ci[k] = out["month_ci"]
        inten[k] = out["inten_mean"]
        inten_ci[k] = out["inten_ci"]
        n_events[k] = out["n_events"]
        n_samp[k] = out["freq_diag"]["n"]
        freq_annual[:, k] = out["freq_annual"]
        print(
            f"  EXP {eid:03d} freq={freq[k]:.3f}/yr "
            f"CI=[{freq_ci[k, 0]:.3f},{freq_ci[k, 1]:.3f}] "
            f"events={int(n_events[k])} months={int(out['n_months'])} "
            f"inten={inten[k]:.4f} Pa",
            flush=True,
        )

    exp = np.array(EXP_IDS, float)
    mat_path = RESULT_MAT / "wwb_precession_results.mat"
    sio.savemat(
        mat_path,
        {
            "expIDs": exp.reshape(-1, 1),
            "ustrw_ts": ts_all,
            "ustrw_anom": anom_all,
            "wwb_freq": freq.reshape(-1, 1),
            "wwb_freq_ci": freq_ci,
            "wwb_freq_annual": freq_annual,
            "wwb_month_rate": month_rate.reshape(-1, 1),
            "wwb_month_ci": month_ci,
            "wwb_inten": inten.reshape(-1, 1),
            "wwb_inten_ci": inten_ci,
            "nEvents": n_events.reshape(-1, 1),
            "nSamples": n_samp.reshape(-1, 1),
            "ciLevel": CI_LEVEL,
            "threshold": THRESHOLD,
            "latRange": np.array(LAT_RANGE),
            "lonRange": np.array(LON_RANGE),
            "preprocess": (
                "WWB: deseason+detrend τx; month if τx'>threshold; "
                "event=contiguous months; mean annual frequency; "
                "95% CI = Student-t on annual event counts (n~100)"
            ),
        },
    )
    print(f"Saved {mat_path}")

    nino = sio.loadmat(NINO3_DIR / "nino3_enso_chi2_results.mat")
    enso = nino["enso_amp"].ravel().astype(float)
    nci = nino["chi2_ci_std"]
    n_lo = enso - nci[:, 0]
    n_hi = nci[:, 1] - enso

    # Only the two dual-axis panels, styled like 01–05
    err_lo = freq - freq_ci[:, 0]
    err_hi = freq_ci[:, 1] - freq
    plot_dual(
        enso, n_lo, n_hi,
        dict(
            key="wwb",
            stem="07_nino3_wwb_dual_axis",
            y=freq,
            elo=err_lo,
            ehi=err_hi,
            color=COLORS["wwb"],
            ylabel=r"WWB freq. (yr$^{-1}$)",
            label="WWB frequency",
        ),
    )

    i_lo = inten - inten_ci[:, 0]
    i_hi = inten_ci[:, 1] - inten
    plot_dual(
        enso, n_lo, n_hi,
        dict(
            key="wwb_inten",
            stem="07_nino3_wwb_intensity_dual_axis",
            y=inten,
            elo=i_lo,
            ehi=i_hi,
            color=COLORS["wwb"],
            ylabel="WWB intensity (Pa)",
            label="WWB intensity",
        ),
    )

    # Remove obsolete wide-style singles if present
    for obsolete in (
        "07_wwb_vs_precession.png",
        "07_wwb_intensity_vs_precession.png",
    ):
        p = RESULT_FIG / obsolete
        if p.is_file():
            p.unlink()
            print("removed", obsolete)

    print("WWB done.")
    print(f"Figures in: {RESULT_FIG}")


if __name__ == "__main__":
    main()
