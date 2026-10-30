"""Path configuration for the Figure 4 mean-state diagnostics.

Model output is not stored in this repository. Set the environment variable
ENSO_DATA_ROOT to the directory that holds the separately downloaded files.
"""
from __future__ import annotations

import os
from pathlib import Path

# .../scripts/fig04_mean_state/code/lib/paths.py
REPO = Path(__file__).resolve().parents[4]
DATA_ROOT = Path(os.environ.get("ENSO_DATA_ROOT", REPO / "data"))

ROOT = Path(__file__).resolve().parents[2]
DATA = DATA_ROOT / "fig04_mean_state"
CODE = ROOT / "code"
RESULT = ROOT / "result"
RESULT_MAT = RESULT / "mat"
RESULT_FIG = RESULT / "figures"
RESULT_TAB = RESULT / "tables"

WALKER_DIR = DATA / "walker"
USTRW_DIR = DATA / "ustrw" / "nc"
THETAO_DIR = DATA / "thetao" / "nc"
IPWP_DIR = DATA / "ipwp"
# Shared with scripts/fig04_nino3_amplitude/compute_nino3_enso.m
NINO3_DIR = DATA_ROOT / "tsw"
Z20_EQ_NC = DATA_ROOT / "fig04_bjerknes" / "prep" / "z20_eq_lon_time_exp.nc"

EXP_IDS = list(range(0, 360, 15))  # 0, 15, ..., 345
N_EXP = len(EXP_IDS)
CI_LEVEL = 0.95


def ensure_result_dirs() -> None:
    for p in (RESULT_MAT, RESULT_FIG, RESULT_TAB):
        p.mkdir(parents=True, exist_ok=True)
