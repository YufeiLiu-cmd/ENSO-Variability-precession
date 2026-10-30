"""Unified journal-style plotting helpers (match 01_walker_edges / 06_combined)."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from lib.paths import EXP_IDS

# Nino3 always the same red; feedbacks use physically motivated colors
COLORS = {
    "nino3": np.array([0.78, 0.16, 0.16]),      # ENSO / warm SST
    "wwb": np.array([0.55, 0.20, 0.35]),        # westerly wind-burst intensity
    "wwv": np.array([0.20, 0.35, 0.55]),        # WWV proxy / Z20
    "mu": np.array([0.12, 0.35, 0.58]),         # atmospheric wind feedback (blue)
    "beta": np.array([0.10, 0.45, 0.55]),       # ocean dynamic response (teal)
    "alpha": np.array([0.10, 0.45, 0.38]),      # thermocline feedback (ocean green)
    "bj": np.array([0.55, 0.25, 0.55]),         # combined loop gain (purple)
}

FIGSIZE_DUAL = (12.5 / 2.54, 4.0 / 2.54)
FIGSIZE_COMBINED = (12.5 / 2.54, 16.0 / 2.54)
FIGSIZE_COMBINED_6 = (12.5 / 2.54, 19.0 / 2.54)  # 6 stacked panels
FIGSIZE_COMBINED_7 = (12.5 / 2.54, 22.0 / 2.54)  # 7 stacked panels
FIGSIZE_SINGLE = (20 / 2.54, 11 / 2.54)


def apply_journal_style() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans", "sans-serif"],
            "font.size": 9,
            "axes.linewidth": 1.1,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "axes.grid": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def _phase_ticklabels(every: int = 2):
    labs = []
    for i, v in enumerate(EXP_IDS):
        labs.append(f"{v:03d}\u00b0" if i % every == 0 else "")
    return labs


def precession_xticks(ax, show_labels: bool = True, every: int = 2) -> None:
    ax.set_xlim(-5, 350)
    ax.set_xticks(EXP_IDS)
    if show_labels:
        ax.set_xticklabels(_phase_ticklabels(every), fontsize=7, rotation=0, ha="center")
    else:
        ax.set_xticklabels([])


def save_figure(fig, stem: str, out_dir: Path) -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"{stem}.png"
    fig.savefig(png, dpi=300, facecolor="white")
    plt.close(fig)
    return png


def plot_dual_nino3(
    exp_ids,
    enso,
    n_lo,
    n_hi,
    y,
    elo,
    ehi,
    *,
    color_right,
    ylabel_right: str,
    label_right: str,
    stem: str,
    out_dir: Path,
    mean_fmt_right: str = ".3g",
    legend_ncol: int = 1,
):
    """Fig.1–3 style: left Nino3 (red), right feedback parameter."""
    apply_journal_style()
    fig, ax = plt.subplots(figsize=FIGSIZE_DUAL)
    fig.subplots_adjust(left=0.14, right=0.86, bottom=0.18, top=0.92)
    col_n = COLORS["nino3"]
    col_r = np.asarray(color_right, float)
    fs_lab, fs_tick, fs_leg = 8, 6.5, 6.5

    ax.set_ylabel("Nino3 amp. (K)", color=tuple(col_n), fontsize=fs_lab, fontweight="bold")
    ax.tick_params(axis="y", colors=tuple(col_n), labelsize=fs_tick, direction="in", length=3)
    h1 = ax.errorbar(
        exp_ids, enso, yerr=np.vstack([n_lo, n_hi]),
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

    ax2 = ax.twinx()
    ax2.set_ylabel(ylabel_right, color=tuple(col_r), fontsize=fs_lab, fontweight="bold")
    ax2.tick_params(axis="y", colors=tuple(col_r), labelsize=fs_tick, direction="in", length=3)
    h2 = ax2.errorbar(
        exp_ids, y, yerr=np.vstack([elo, ehi]),
        fmt="-s", color=col_r, lw=1.3, ms=3.5, mfc=col_r, mec="w", mew=0.35,
        capsize=2.0, elinewidth=0.9, label=label_right, zorder=3,
    )
    h2m = ax2.axhline(np.nanmean(y), ls="--", color=col_r, lw=0.9, alpha=0.85)
    r0 = np.nanmin(y - elo)
    r1 = np.nanmax(y + ehi)
    rpad = max((r1 - r0) * 0.16, 1e-6)
    ax2.set_ylim(r0 - rpad, r1 + rpad)
    ax2.spines["right"].set_color(tuple(col_r))
    ax2.spines["right"].set_linewidth(1.2)
    ax2.spines["left"].set_visible(False)

    ax.set_xlabel("")
    precession_xticks(ax, show_labels=True, every=2)
    ax.tick_params(axis="x", labelsize=fs_tick, direction="in", length=3)

    if legend_ncol >= 2:
        # 2-col: series | means — avoids covering curves at upper-left
        handles = [h1, h2, h1m, h2m]
        labels = [
            "Nino3 amplitude",
            label_right,
            f"mean = {np.nanmean(enso):.2f} K",
            f"mean = {np.nanmean(y):{mean_fmt_right}}",
        ]
        leg_kw = dict(
            loc="upper left",
            ncol=2,
            fontsize=fs_leg,
            handlelength=1.2,
            labelspacing=0.18,
            columnspacing=1.2,
            borderaxespad=0.2,
            handletextpad=0.3,
            frameon=True,
            fancybox=False,
            edgecolor="none",
            facecolor="white",
            framealpha=0.9,
        )
    else:
        handles = [h1, h1m, h2, h2m]
        labels = [
            "Nino3 amplitude",
            f"mean = {np.nanmean(enso):.2f} K",
            label_right,
            f"mean = {np.nanmean(y):{mean_fmt_right}}",
        ]
        leg_kw = dict(
            loc="upper left",
            frameon=False,
            fontsize=fs_leg,
            handlelength=1.2,
            labelspacing=0.2,
            borderaxespad=0.2,
            handletextpad=0.35,
        )
    ax.legend(handles, labels, **leg_kw)

    fig.set_size_inches(*FIGSIZE_DUAL)
    return save_figure(fig, stem, out_dir)


def plot_single_errorbar(
    exp_ids, y, elo, ehi, *, color, ylabel, label, stem, out_dir, mean_fmt=".4g",
):
    apply_journal_style()
    fig, ax = plt.subplots(figsize=FIGSIZE_SINGLE, constrained_layout=True)
    c = np.asarray(color, float)
    h = ax.errorbar(
        exp_ids, y, yerr=np.vstack([elo, ehi]),
        fmt="-o", color=c, lw=2.0, ms=7, mfc=c, mec="w", mew=0.8,
        capsize=4, elinewidth=1.2, label=label,
    )
    hm = ax.axhline(np.nanmean(y), ls="--", color=c, lw=1.6)
    ax.set_xlabel(r"Precession phase ($^\circ$)", fontweight="bold")
    ax.set_ylabel(ylabel, fontweight="bold")
    precession_xticks(ax, every=1)
    ax.set_xticklabels([f"{v:03d}" for v in exp_ids], fontsize=10, rotation=45, ha="right")
    y0 = np.nanmin(y - elo)
    y1 = np.nanmax(y + ehi)
    pad = max((y1 - y0) * 0.12, 1e-6)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.legend(
        [h, hm], [label, f"mean = {np.nanmean(y):{mean_fmt}}"],
        loc="best", frameon=False,
    )
    return save_figure(fig, stem, out_dir)


def load_nino3(mat_path):
    import scipy.io as sio

    nino = sio.loadmat(mat_path)
    enso = nino["enso_amp"].ravel().astype(float)
    ci = nino["chi2_ci_std"]
    n_lo = enso - ci[:, 0]
    n_hi = ci[:, 1] - enso
    return enso, n_lo, n_hi
