"""Unified journal-style plotting helpers."""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# Nino3 always the same red; other parameters use physically motivated colors
COLORS = {
    "nino3": np.array([0.78, 0.16, 0.16]),          # warm SST / ENSO
    "walker": np.array([0.15, 0.40, 0.65]),          # atmospheric circulation blue
    "ustrw": np.array([0.10, 0.45, 0.55]),           # wind-stress teal
    "ipwp": np.array([0.85, 0.45, 0.10]),            # warm-pool orange
    "mld": np.array([0.10, 0.45, 0.38]),             # ocean mixed-layer green
    "strat": np.array([0.45, 0.25, 0.60]),           # stratification purple
    "wwb": np.array([0.55, 0.20, 0.35]),             # westerly wind burst maroon
    "wwv": np.array([0.20, 0.35, 0.55]),             # warm-water volume / Z20 blue
    "west": np.array([0.12, 0.47, 0.71]),
    "east": np.array([0.85, 0.33, 0.10]),
    "center": np.array([0.17, 0.63, 0.17]),
}


# Wide, short panels so precession phases are not cramped
FIGSIZE_DUAL = (40 / 2.54, 8.0 / 2.54)  # cm: longer & shorter
FIGSIZE_COMBINED = (40 / 2.54, 34 / 2.54)


def apply_journal_style() -> None:
    mpl.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "font.family": "sans-serif",
            "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans", "sans-serif"],
            "font.size": 14,
            "axes.labelsize": 17,
            "axes.titlesize": 16,
            "xtick.labelsize": 12,
            "ytick.labelsize": 13,
            "legend.fontsize": 11,
            "axes.linewidth": 1.8,
            "xtick.major.width": 1.4,
            "ytick.major.width": 1.4,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "axes.grid": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def precession_xticks(ax, exp_ids=None) -> None:
    if exp_ids is None:
        exp_ids = list(range(0, 360, 15))
    ax.set_xlim(-8, 353)
    ax.set_xticks(exp_ids)
    ax.set_xticklabels([f"{v:03d}" for v in exp_ids], rotation=45, ha="right")


def dual_axis_errorbar(
    exp_ids,
    y_left,
    err_left_low,
    err_left_high,
    y_right,
    err_right_low,
    err_right_high,
    *,
    color_left,
    color_right,
    ylabel_left: str,
    ylabel_right: str,
    label_left: str,
    label_right: str,
    mean_left: float | None = None,
    mean_right: float | None = None,
    xlabel: str = r"Precession phase ($^\circ$)",
    figsize=FIGSIZE_DUAL,
    right_has_err: bool = True,
):
    """Dual-y journal panel with Nino3-style left axis."""
    apply_journal_style()
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

    ax.yaxis.set_tick_params(colors=tuple(color_left))
    h1 = ax.errorbar(
        exp_ids,
        y_left,
        yerr=np.vstack([err_left_low, err_left_high]),
        fmt="-o",
        color=color_left,
        lw=2.0,
        ms=7,
        mfc=color_left,
        mec="white",
        mew=0.8,
        capsize=4,
        label=label_left,
        zorder=3,
    )
    h1m = None
    if mean_left is not None:
        h1m = ax.axhline(mean_left, ls="--", color=color_left, lw=1.6, label=f"mean = {mean_left:.2f}")
    ax.set_ylabel(ylabel_left, fontweight="bold", color=tuple(color_left))
    y0 = np.nanmin(np.asarray(y_left) - np.asarray(err_left_low))
    y1 = np.nanmax(np.asarray(y_left) + np.asarray(err_left_high))
    pad = max((y1 - y0) * 0.12, 0.05)
    ax.set_ylim(y0 - pad, y1 + pad)

    ax2 = ax.twinx()
    ax2.yaxis.set_tick_params(colors=tuple(color_right))
    if right_has_err:
        h2 = ax2.errorbar(
            exp_ids,
            y_right,
            yerr=np.vstack([err_right_low, err_right_high]),
            fmt="-s",
            color=color_right,
            lw=2.0,
            ms=7,
            mfc=color_right,
            mec="white",
            mew=0.8,
            capsize=4,
            label=label_right,
            zorder=3,
        )
    else:
        h2 = ax2.plot(
            exp_ids,
            y_right,
            "-s",
            color=color_right,
            lw=2.0,
            ms=7,
            mfc=color_right,
            mec="white",
            mew=0.8,
            label=label_right,
            zorder=3,
        )[0]
        err_right_low = np.zeros_like(np.asarray(y_right), dtype=float)
        err_right_high = np.zeros_like(np.asarray(y_right), dtype=float)

    h2m = None
    if mean_right is not None:
        h2m = ax2.axhline(
            mean_right,
            ls="--",
            color=color_right,
            lw=1.6,
            label=f"mean = {mean_right:.3g}",
        )
    ax2.set_ylabel(ylabel_right, fontweight="bold", color=tuple(color_right))
    r0 = np.nanmin(np.asarray(y_right) - np.asarray(err_right_low))
    r1 = np.nanmax(np.asarray(y_right) + np.asarray(err_right_high))
    rpad = max((r1 - r0) * 0.12, 1e-6)
    ax2.set_ylim(r0 - rpad, r1 + rpad)

    ax.set_xlabel(xlabel, fontweight="bold")
    precession_xticks(ax, exp_ids)
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)
    for spine in ax2.spines.values():
        spine.set_linewidth(1.8)

    handles = [h1]
    labels = [label_left]
    if h1m is not None:
        handles.append(h1m)
        labels.append(h1m.get_label())
    handles.append(h2)
    labels.append(label_right)
    if h2m is not None:
        handles.append(h2m)
        labels.append(h2m.get_label())
    ax.legend(handles, labels, loc="best", frameon=False, fontsize=11)
    return fig, ax, ax2


def save_figure(fig, stem: str, out_dir: Path) -> Path:
    """Save PNG only (no PDF)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f"{stem}.png"
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return png


def single_errorbar(
    exp_ids,
    y,
    err_low,
    err_high,
    *,
    color,
    ylabel: str,
    label: str,
    mean_val: float | None = None,
    xlabel: str = r"Precession phase ($^\circ$)",
    figsize=FIGSIZE_DUAL,
):
    apply_journal_style()
    fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)
    h = ax.errorbar(
        exp_ids,
        y,
        yerr=np.vstack([err_low, err_high]),
        fmt="-o",
        color=color,
        lw=2.0,
        ms=7,
        mfc=color,
        mec="white",
        mew=0.8,
        capsize=4,
        label=label,
    )
    handles = [h]
    labels = [label]
    if mean_val is not None:
        hm = ax.axhline(mean_val, ls="--", color=color, lw=1.6)
        handles.append(hm)
        labels.append(f"mean = {mean_val:.4g}")
    ax.set_xlabel(xlabel, fontweight="bold")
    ax.set_ylabel(ylabel, fontweight="bold")
    precession_xticks(ax, exp_ids)
    y0 = np.nanmin(np.asarray(y) - np.asarray(err_low))
    y1 = np.nanmax(np.asarray(y) + np.asarray(err_high))
    pad = max((y1 - y0) * 0.12, 1e-6)
    ax.set_ylim(y0 - pad, y1 + pad)
    ax.legend(handles, labels, loc="best", frameon=False)
    return fig, ax
