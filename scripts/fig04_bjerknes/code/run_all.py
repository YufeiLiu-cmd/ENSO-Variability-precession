#!/usr/bin/env python3
"""Run the Bjerknes-feedback calculations used by Figure 4i–l."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

CODE = Path(__file__).resolve().parent
ROOT = CODE.parent

SCRIPTS = [
    "01_wind_feedback.py",
    "02_thermocline_feedback.py",
    "03_ocean_response_beta.py",
    "04_bjerknes_index.py",
]


def main():
    env = dict(os.environ)
    env["MPLCONFIGDIR"] = str(ROOT / "result" / ".mplconfig")
    env["MPLBACKEND"] = "Agg"
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

    for name in SCRIPTS:
        print("\n" + "=" * 60)
        print(f"Running {name}")
        print("=" * 60, flush=True)
        r = subprocess.run([sys.executable, str(CODE / name)], cwd=str(CODE), env=env)
        if r.returncode != 0:
            raise SystemExit(f"{name} failed with code {r.returncode}")
    print("\nAll done.")


if __name__ == "__main__":
    main()
