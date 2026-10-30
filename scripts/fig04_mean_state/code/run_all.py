#!/usr/bin/env python3
"""Run the Figure 4 mean-state diagnostics, then the published two-column figure.

The two-column figure also reads Bjerknes-feedback results written by
scripts/fig04_bjerknes/code/run_all.py. Run that driver first.
"""
from __future__ import annotations

import runpy
import sys
from pathlib import Path

SCRIPTS = [
    "01_walker_edges.py",
    "02_ustrw.py",
    "03_ipwp_plot.py",
    "04_mld.py",
    "05_stratification.py",
    "07_wwb.py",
    "08_wwv.py",
    "09_combined_meanstate_feedbacks_2col.py",
]


def main():
    code_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(code_dir))
    for name in SCRIPTS:
        print("\n" + "=" * 60)
        print(f"Running {name}")
        print("=" * 60, flush=True)
        runpy.run_path(str(code_dir / name), run_name="__main__")
    print("\nAll scripts finished.")


if __name__ == "__main__":
    main()
