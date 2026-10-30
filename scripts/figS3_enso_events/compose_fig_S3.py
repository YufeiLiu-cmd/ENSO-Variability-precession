#!/usr/bin/env python3
"""
Supplementary Figure S3 — unified style from analysis data.

  (a) ENSO event counts
  (b) Niño3 σ with 95% CI
  (c) Power spectral density (computed from the same Niño3 series)
  (d) Niño3 time series at 090°  (full width)
  (e) Niño3 time series at 270°  (full width)
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

HERE = Path(__file__).resolve().parent
EVENT_DIR = HERE

EXP_IDS = np.arange(0, 346, 15, dtype=int)
N_TIME = 1200
N_YEARS = N_TIME / 12.0
THRESH = 0.5
MIN_DUR = 5
FS = 12.0  # samples per year (monthly)

COL_EN = "#C0392B"
COL_LN = "#2471A3"
COL_LINE = "#2C3E50"
COL_GRID = "#E8E8E8"

CRITERIA_TEXT = (
    f"Criteria: 3-mo running-mean Niño3\n"
    f"|SSTA| > {THRESH} K for ≥ {MIN_DUR} months"
)


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
            "font.size": 12,
            "axes.labelsize": 13,
            "axes.titlesize": 13,
            "xtick.labelsize": 10,
            "ytick.labelsize": 11,
            "legend.fontsize": 10,
            "axes.linewidth": 1.2,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def find_events(series: np.ndarray, sign: int, thresh: float, min_dur: int) -> list[dict]:
    mask = series > thresh if sign > 0 else series < -thresh
    events: list[dict] = []
    n = len(series)
    i = 0
    while i < n:
        if not mask[i]:
            i += 1
            continue
        j = i
        while j < n and mask[j]:
            j += 1
        if j - i >= min_dur:
            seg = series[i:j]
            peak = float(np.max(seg) if sign > 0 else np.min(seg))
            events.append({"start": i, "end": j - 1, "duration": j - i, "peak": peak})
        i = j
    return events


def load_event_inputs() -> dict:
    nino3 = np.load(EVENT_DIR / "nino3_index_24exp.npy")
    assert nino3.shape == (N_TIME, 24)

    n_elnino = np.zeros(24, dtype=int)
    n_lanina = np.zeros(24, dtype=int)
    std_nino3 = np.full(24, np.nan)
    chi2_ci = np.full((24, 2), np.nan)

    with (EVENT_DIR / "enso_event_stats.csv").open() as f:
        rows = list(csv.DictReader(f))
    for k, row in enumerate(rows):
        n_elnino[k] = int(row["n_elnino"])
        n_lanina[k] = int(row["n_lanina"])
        std_nino3[k] = float(row["std_nino3"])
        chi2_ci[k, 0] = float(row["chi2_ci95_std_lo"])
        chi2_ci[k, 1] = float(row["chi2_ci95_std_hi"])

    all_en, all_ln = [], []
    for k in range(24):
        all_en.append(find_events(nino3[:, k], +1, THRESH, MIN_DUR))
        all_ln.append(find_events(nino3[:, k], -1, THRESH, MIN_DUR))

    return {
        "expIDs": EXP_IDS,
        "nino3_all": nino3,
        "n_elnino": n_elnino,
        "n_lanina": n_lanina,
        "std_nino3": std_nino3,
        "chi2_ci_std": chi2_ci,
        "all_en_events": all_en,
        "all_ln_events": all_ln,
    }


def compute_psd(nino3: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Welch PSD for each experiment; same settings as original powerspectrum.m."""
    n_exp = nino3.shape[1]
    win = signal.windows.blackman(360)
    pxx_list = []
    f_ref = None
    for k in range(n_exp):
        f, pxx = signal.welch(
            nino3[:, k],
            fs=FS,
            window=win,
            noverlap=180,
            nfft=None,
            scaling="density",
        )
        if f_ref is None:
            f_ref = f
        pxx_list.append(pxx)
    return f_ref, np.vstack(pxx_list)


def panel_label(ax, lab: str) -> None:
    ax.text(
        0.015,
        0.97,
        lab,
        transform=ax.transAxes,
        fontsize=14,
        fontweight="bold",
        va="top",
        ha="left",
        zorder=20,
        color="k",
    )


def style_ax(ax) -> None:
    ax.yaxis.grid(True, color=COL_GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for sp in ax.spines.values():
        sp.set_linewidth(1.2)


def draw_counts(ax, r: dict) -> None:
    exp = r["expIDs"]
    n = len(exp)
    x = np.arange(n)
    w = 0.38
    b1 = ax.bar(
        x - w / 2, r["n_elnino"], w, color=COL_EN, edgecolor="white", lw=0.4, label="El Niño", zorder=3
    )
    b2 = ax.bar(
        x + w / 2, r["n_lanina"], w, color=COL_LN, edgecolor="white", lw=0.4, label="La Niña", zorder=3
    )
    for bar, val in zip(b1, r["n_elnino"]):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.3,
            f"{int(val)}",
            ha="center",
            va="bottom",
            fontsize=7,
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
            fontsize=7,
            color=COL_LN,
            fontweight="bold",
        )
    majors = np.where(exp % 30 == 0)[0]
    ax.set_xticks(majors)
    ax.set_xticklabels([f"{int(exp[i]):03d}" for i in majors])
    ax.set_ylabel("Number of events")
    ax.set_xlim(-0.45, n - 0.55)
    ax.set_ylim(0, 40)
    style_ax(ax)
    ax.legend(frameon=False, loc="upper left", fontsize=9.5, ncol=2, bbox_to_anchor=(0.06, 1.0))
    ax.text(
        0.98,
        0.98,
        CRITERIA_TEXT,
        transform=ax.transAxes,
        fontsize=7.5,
        ha="right",
        va="top",
        color="#333333",
        bbox=dict(
            boxstyle="round,pad=0.28",
            facecolor="white",
            edgecolor="#CCCCCC",
            alpha=0.92,
            linewidth=0.8,
        ),
        zorder=10,
    )


def draw_sigma(ax, r: dict) -> None:
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
        lw=1.8,
        ms=5.8,
        mfc=COL_LINE,
        mec="white",
        mew=0.9,
        capsize=3.0,
        elinewidth=1.1,
        ecolor=COL_LINE,
        label=r"Niño3 $\sigma$ (3-mo mean)",
        zorder=3,
    )
    majors = exp[exp % 30 == 0]
    ax.set_xticks(majors)
    ax.set_xticklabels([f"{int(e):03d}" for e in majors])
    ax.set_ylabel("Niño3 standard deviation (K)")
    ax.set_xlim(-2, 347)
    style_ax(ax)
    ax.legend(frameon=False, loc="upper left", fontsize=9.5, bbox_to_anchor=(0.16, 1.0))
    ax.text(
        0.98,
        0.02,
        r"95% CI ($\chi^2$)",
        transform=ax.transAxes,
        fontsize=8,
        ha="right",
        va="bottom",
        color="#555",
    )


def draw_timeseries(ax, r: dict, k: int, *, show_note: bool = False) -> None:
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
    ylim = max(2.5, float(np.nanpercentile(np.abs(idx), 99.5) * 1.15))
    ax.set_ylim(-ylim, ylim)
    style_ax(ax)
    ax.text(
        0.99,
        0.90,
        f"Precession {exp:03d}°   El Niño = {r['n_elnino'][k]}   La Niña = {r['n_lanina'][k]}",
        transform=ax.transAxes,
        fontsize=11,
        ha="right",
        va="top",
    )
    if show_note:
        ax.text(
            0.99,
            0.08,
            "3-month running-mean Niño3",
            transform=ax.transAxes,
            fontsize=9,
            ha="right",
            va="bottom",
            color="#555555",
        )


def draw_psd(ax, fig, f: np.ndarray, pxx: np.ndarray) -> None:
    """PSD filled by dense contourf (smooth colorbar interpolation), A2 style."""
    fmax = 0.6
    mask = (f >= 0) & (f <= fmax)
    f_plot = f[mask]
    pxx_plot = np.maximum(pxx[:, mask], 1e-12)

    phases = np.arange(0, 361, 15, dtype=float)
    z = np.vstack([pxx_plot, pxx_plot[0:1, :]])
    z_log = np.log(z)

    # Dense levels → smooth interpolated fill keyed to the colorbar
    levels = np.linspace(-1.6, 1.6, 321)
    cf = ax.contourf(
        f_plot,
        phases,
        z_log,
        levels=levels,
        cmap="RdPu",
        extend="both",
        antialiased=False,
    )

    cs = ax.contour(
        f_plot,
        phases,
        z,
        levels=np.arange(0.5, 8.0, 0.5),
        colors="k",
        linewidths=0.35,
    )
    ax.clabel(cs, fontsize=8, fmt="%1.1f", inline=True, levels=np.arange(1.0, 8.0, 1.0))

    ax.set_xlim(0.0, fmax)
    ax.set_ylim(-2.0, 362.0)
    ax.set_xticks(np.arange(0.0, 0.7, 0.1))
    yticks = np.arange(0, 361, 30)
    ax.set_yticks(yticks)
    ax.set_yticklabels([f"{i}°" for i in yticks])
    ax.set_ylabel("Precession phase")
    ax.tick_params(axis="both", length=4)

    for y0, h in ((75, 30), (255, 30)):
        ax.add_patch(
            plt.Rectangle(
                (0.0, y0),
                fmax,
                h,
                fill=False,
                edgecolor="#00A8A8",
                linewidth=1.8,
                linestyle="--",
                zorder=5,
            )
        )

    for sp in ax.spines.values():
        sp.set_linewidth(1.2)

    cbar = fig.colorbar(cf, ax=ax, fraction=0.035, pad=0.01)
    cbar.set_label(r"log(PSD) [$^\circ$C$^2$ year]")
    cbar.set_ticks(np.arange(-1.5, 1.51, 0.5))
    cbar.ax.tick_params(labelsize=9)


def compose() -> None:
    setup_style()
    r = load_event_inputs()
    f, pxx = compute_psd(r["nino3_all"])
    np.savez(HERE / "psd_from_nino3.npz", f=f, pxx=pxx, expIDs=EXP_IDS)

    k_090 = int(np.where(r["expIDs"] == 90)[0][0])
    k_270 = int(np.where(r["expIDs"] == 270)[0][0])

    fig = plt.figure(figsize=(13.0, 13.6))
    gs = fig.add_gridspec(
        3,
        1,
        height_ratios=[1.0, 1.12, 1.72],
        hspace=0.08,
        left=0.07,
        right=0.97,
        top=0.988,
        bottom=0.045,
    )

    # (a)|(b) with tighter horizontal gap
    gs_ab = gs[0].subgridspec(1, 2, wspace=0.23)
    ax_a = fig.add_subplot(gs_ab[0, 0])
    draw_counts(ax_a, r)
    panel_label(ax_a, "(a)")

    ax_b = fig.add_subplot(gs_ab[0, 1])
    draw_sigma(ax_b, r)
    panel_label(ax_b, "(b)")

    ax_c = fig.add_subplot(gs[1])
    draw_psd(ax_c, fig, f, pxx)
    panel_label(ax_c, "(c)")

    gs_de = gs[2].subgridspec(2, 1, hspace=0.05)
    ax_d = fig.add_subplot(gs_de[0, 0])
    draw_timeseries(ax_d, r, k_090, show_note=True)
    panel_label(ax_d, "(d)")

    ax_e = fig.add_subplot(gs_de[1, 0], sharex=ax_d)
    draw_timeseries(ax_e, r, k_270)
    panel_label(ax_e, "(e)")

    out_png = HERE / "Fig_S3_enso_events_psd.png"
    out_pdf = HERE / "Fig_S3_enso_events_psd.pdf"
    fig.savefig(out_png, dpi=300, facecolor="white")
    fig.savefig(out_pdf, dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {out_png}")
    print(f"Saved: {out_pdf}")


if __name__ == "__main__":
    compose()
