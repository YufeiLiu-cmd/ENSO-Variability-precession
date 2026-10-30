"""Path configuration for the Bjerknes-feedback diagnostics in Figure 4.

Model output is not stored in this repository. Set ENSO_DATA_ROOT to the
directory that holds the separately downloaded files.
"""
from __future__ import annotations

import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[4]
DATA_ROOT = Path(os.environ.get("ENSO_DATA_ROOT", REPO / "data"))

ROOT = Path(__file__).resolve().parents[2]
DATA = DATA_ROOT / "fig04_bjerknes"
CODE = ROOT / "code"
RESULT = ROOT / "result"
RESULT_MAT = RESULT / "mat"
RESULT_FIG = RESULT / "figures"
RESULT_TAB = RESULT / "tables"

TSW_DIR = DATA / "tsw" / "nc"
USTRW_DIR = DATA / "ustrw" / "nc"
THETAO_DIR = DATA / "thetao" / "nc"
PREP_DIR = DATA / "prep"
NINO3_DIR = DATA_ROOT / "tsw"

Z20_NC = PREP_DIR / "z20_eq_lon_time_exp.nc"
SST_EQ_NC = PREP_DIR / "sst_eq_lon_time_exp.nc"
NINO3_MAT = NINO3_DIR / "nino3_enso_chi2_results.mat"

EXP_IDS = list(range(0, 360, 15))  # 0, 15, ..., 345
N_EXP = len(EXP_IDS)
N_TIME = 1200
CI_LEVEL = 0.95


def ensure_result_dirs() -> None:
    for p in (RESULT_MAT, RESULT_FIG, RESULT_TAB, PREP_DIR):
        p.mkdir(parents=True, exist_ok=True)
