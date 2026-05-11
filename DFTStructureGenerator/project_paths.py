from __future__ import annotations

import os
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "Data"
CSV_DIR = DATA_DIR / "csvs"
TS_DATA_DIR = DATA_DIR / "TS"
DESCRIPTOR_DIR = DATA_DIR / "descriptor"
FIGURE_DIR = REPO_ROOT / "Figure"

DEFAULT_RAW_CALC_ROOT = Path("E:/work/B_Cl_Nu")
RAW_CALC_ROOT = Path(
    os.environ.get("BORYLXAT_RAW_CALC_ROOT", DEFAULT_RAW_CALC_ROOT)
).expanduser()


def repo_path(*parts: str) -> Path:
    return REPO_ROOT.joinpath(*parts)


def raw_calc_path(*parts: str) -> Path:
    return RAW_CALC_ROOT.joinpath(*parts)


def raw_calc_file(*parts: str) -> str:
    return str(raw_calc_path(*parts))
