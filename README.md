# BorylXAT-DB: a transition-state database for boryl-radical-mediated C-Cl atom transfer

[中文说明](docs/zh-CN/README.md) | [Database Structure](Database_Structure.md) | [数据库结构说明](docs/zh-CN/Database_Structure.md)

This repository contains the workflow, released data interfaces, and analysis notebooks for **BorylXAT-DB**, a quantum-chemistry dataset for Lewis-base-coordinated boryl-radical-mediated C-Cl atom transfer.

Public code repository: <https://github.com/jackie-illiilli/BorylXAT-DB>  
Zenodo archive for large files and released databases: <https://doi.org/10.5281/zenodo.20134535>

## What This Repository Covers

The project combines:

- reaction-space enumeration from borane radical, Lewis base, and chlorinated substrate libraries
- RDKit and xTB based conformer generation
- Gaussian job generation for ground states, constrained optimizations, transition-state searches, single-point corrections, and IRC validation
- structured database construction in ASE SQLite and Parquet formats
- descriptor extraction, benchmark analysis, BEP analysis, and machine-learning modeling
- reviewer-runnable notebooks for the main paper and revision analyses

The manuscript-facing designed space contains:

- `55` borane radicals
- `386` Lewis bases in the manuscript count
- `179` chlorinated substrates

The released structural database currently contains:

- `BorylXAT-DB.db`: ASE SQLite database
- `BorylXAT-DB.parquet`: flattened Parquet export
- `BorylXAT-DB_qh_update.db`: quasi-harmonic thermochemistry update database used by the QHARM revision notebooks

The ASE database categories are:

| Category | Count |
| --- | ---: |
| `B` | 55 |
| `LB` | 387 |
| `Cl` | 179 |
| `complex_r` | 20010 |
| `complex_p` | 20010 |
| `c_radical` | 179 |
| `ts` | 8980 |

The standalone `LB` count is `387` in the database because `LB_00623` is retained for provenance even though it does not survive the manuscript filtering into stable B-LB complexes or TS entries.

## Dependency Files

To address reproducibility for reviewers and readers, the repository now includes three machine-readable dependency definitions:

- `pyproject.toml`: package metadata plus the core Python stack for the reusable modules and notebooks
- `requirements.txt`: pip-style environment for reviewer-runnable notebook analysis
- `environment.yml`: conda environment, recommended when installing `rdkit`

The Python packages used by the public workflows are:

- `ase`
- `catboost`
- `ipykernel`
- `jupyterlab`
- `matplotlib`
- `morfeus-ml`
- `numpy`
- `openpyxl`
- `pandas`
- `pyarrow`
- `rdkit`
- `scikit-learn`
- `scipy`
- `seaborn`
- `tqdm`
- `xgboost` for some revision baseline notebooks

External software expected only for full raw-calculation reruns:

- Gaussian
- xTB / CREST
- a PBS-style HPC environment for generated batch scripts

## Environment Setup

### Option 1: pip / venv

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

On Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

### Option 2: conda

```bash
conda env create -f environment.yml
conda activate borylxat-db
```

### Option 3: uv

```bash
uv sync --extra revision
```

Use the `revision` extra when you want optional `xgboost` support for the revision baseline notebook.

## Data Files and What They Are For

The repository already includes the main released data files in the current workspace. If they are missing after a fresh clone, download them from Zenodo and place them in the repository root.

| File | Purpose | Needed for |
| --- | --- | --- |
| `BorylXAT-DB.db` | released ASE structural database | `3_Build_DataBase.ipynb`, `5_Modeling.ipynb`, figure regeneration, direct database inspection |
| `BorylXAT-DB.parquet` | flattened reaction/structure export | downstream tabular analysis and quick inspection |
| `BorylXAT-DB_qh_update.db` | QHARM-updated database | QHARM revision notebooks |
| `Data/TS/Borane_all.csv` | curated public reaction table | modeling and BEP workflows |
| `Data/TS/result_filter_extended.csv` | extended reaction-space table used by revision analyses | thermodynamic-filter and coverage notebooks |
| `Data/csvs/Benchmark.csv` | benchmark summary table | `4_Benchmark.ipynb` reviewer-runnable section |
| `Data/csvs/Experiment.csv` | experimental validation input table | `revision_JACS_Experiment.ipynb` |
| `Data/descriptor/*.pkl` | precomputed descriptor maps | `5_Modeling.ipynb`, `5_Modeling_xTB.ipynb`, revision modeling notebooks |

If needed, the two main released database files can be downloaded with:

```bash
curl -L -o BorylXAT-DB.db "https://zenodo.org/records/20134535/files/BorylXAT-DB.db?download=1"
curl -L -o BorylXAT-DB.parquet "https://zenodo.org/records/20134535/files/BorylXAT-DB.parquet?download=1"
```

On Windows PowerShell:

```powershell
Invoke-WebRequest -Uri "https://zenodo.org/records/20134535/files/BorylXAT-DB.db?download=1" -OutFile "BorylXAT-DB.db"
Invoke-WebRequest -Uri "https://zenodo.org/records/20134535/files/BorylXAT-DB.parquet?download=1" -OutFile "BorylXAT-DB.parquet"
```

## Notebook Guide

### Main notebooks

| Notebook | Scope | Reviewer-runnable from released files? | Main required inputs |
| --- | --- | --- | --- |
| `1_Calc_Reactant.ipynb` | reactant preparation, enumeration, xTB conformers, Gaussian setup | no, documents provenance workflow | raw calculation workspace, Gaussian, xTB/CREST |
| `2_Calc_TS.ipynb` | TS guesses, constrained optimization, TS search, IRC workflow | no, documents provenance workflow | raw calculation workspace, Gaussian |
| `3_Build_DataBase.ipynb` | database build plus released-database inspection | yes for inspection sections | `BorylXAT-DB.db`, `BorylXAT-DB.parquet` |
| `4_Benchmark.ipynb` | benchmark analysis | yes for released analysis sections | `Data/csvs/Benchmark.csv` |
| `5_Modeling.ipynb` | DFT-descriptor CatBoost modeling | yes | `Data/TS/Borane_all.csv`, `Data/descriptor/BNdes_new2.pkl`, `Data/descriptor/Cldes_new2.pkl` |
| `5_Modeling_xTB.ipynb` | xTB-descriptor modeling | yes | `Data/TS/Borane_all.csv`, `Data/descriptor/BNdes_xtb.pkl`, `Data/descriptor/Cldes_xtb.pkl` |
| `6_Draw_Figures.ipynb` | consolidated manuscript and SI figures | yes | released CSV, descriptor, database, and output files |

### Revision notebooks added for reviewer response

| Notebook | Purpose | Extra notes |
| --- | --- | --- |
| `revision_ml_baseline_ablation_bep_dependence.ipynb` | ML baselines, ablations, BEP-dependence analysis | `xgboost` is optional but recommended |
| `revision_thermodynamic_filter_BEP_ML.ipynb` | thermodynamic-filter auxiliary-set analysis | uses released/revision CSV outputs |
| `revision_JACS_Experiment.ipynb` | external experimental trend validation | uses `Data/csvs/Experiment.csv` and `output/jacs_experiment/` |
| `revision_qharm_barrier_analysis.ipynb` | RRHO vs QHARM barrier comparison | needs `BorylXAT-DB_qh_update.db` |
| `revision_qharm_modeling.ipynb` | QHARM-updated modeling analysis | prefers QHARM descriptor maps |
| `revision_basis_set_bond_length_comparison.ipynb` | 6-31G(d) vs 6-31+G(d) geometry sensitivity | uses released revision output CSV files |
| `7_QHARM_Thermodynamic_Filter_Coverage.ipynb` | QHARM favorable-domain coverage summary | uses released/revision CSV outputs |

The notebooks use the following review tags:

| Tag | Meaning |
| --- | --- |
| `[REVIEWER-RUNNABLE]` | should run from files included in this repository or the released Zenodo package |
| `[RAW-GAUSSIAN/E:/work]` | depends on the historical raw Gaussian/xTB workspace |
| `[OPTIONAL-DESCRIPTOR-GENERATION]` | expensive descriptor generation that can usually be skipped because precomputed descriptor files are shipped |

## Repository Layout

```text
.
|-- DFTStructureGenerator/     reusable Python modules for parsing, descriptors, database building, and plotting
|-- Data/
|   |-- ChemDraw/              ChemDraw source files for reactant libraries
|   |-- csvs/                  intermediate and released summary CSV tables
|   |-- descriptor/            released descriptor pickle files
|   `-- TS/                    released TS-level CSV tables and coordinates
|-- Figure/                    manuscript and SI figure exports
|-- output/                    revision analysis outputs and executed notebooks
|-- 1_Calc_Reactant.ipynb
|-- 2_Calc_TS.ipynb
|-- 3_Build_DataBase.ipynb
|-- 4_Benchmark.ipynb
|-- 5_Modeling.ipynb
|-- 5_Modeling_xTB.ipynb
|-- 6_Draw_Figures.ipynb
|-- revision_*.ipynb
|-- 7_QHARM_Thermodynamic_Filter_Coverage.ipynb
`-- Database_Structure.md
```

## Core Python Modules

`DFTStructureGenerator/` contains the reusable workflow code:

- `borane_xat_workflow.py`: reaction-site detection, component combination generation, DFT job preparation, TS structure generation, and summary export
- `Build_DataBase.py`: builds `BorylXAT-DB.db` and `BorylXAT-DB.parquet`
- `descriptor.py`: descriptor-map generation and ML feature assembly
- `thermochemistry.py`: RRHO/QHARM database access helpers used by revision notebooks
- `logfile_process.py`: Gaussian log parsing
- `xtb_process.py`: xTB/CREST conformer workflow helpers
- `mol_manipulation.py`: geometry transforms and Gaussian input helpers
- `FormatConverter.py`: conversions among RDKit, XYZ, and Gaussian-related formats
- `draw.py`: plotting helpers used across notebooks
- `project_paths.py`: repository-aware path helpers and the raw-workspace environment-variable hook

## Raw-Calculation Path Notes

The historical raw-calculation workspace was organized under `E:/work/...`. Public reviewer-runnable sections avoid hard dependence on that directory, but full provenance reruns still expect a similar layout.

The reusable path helper is:

- environment variable: `BORYLXAT_RAW_CALC_ROOT`
- default fallback: `E:/work/B_Cl_Nu`

If you need to point the provenance notebooks to a local archive, set:

```powershell
$env:BORYLXAT_RAW_CALC_ROOT = "D:\\path\\to\\raw\\workspace"
```

## Typical Usage

For full provenance reruns:

```text
1_Calc_Reactant.ipynb
2_Calc_TS.ipynb
3_Build_DataBase.ipynb
4_Benchmark.ipynb
5_Modeling.ipynb
5_Modeling_xTB.ipynb
6_Draw_Figures.ipynb
```

For reviewer-style reuse from released files only:

1. Use `3_Build_DataBase.ipynb` for database inspection only.
2. Use `4_Benchmark.ipynb` from the checked-in `Benchmark.csv` loading section onward.
3. Use `5_Modeling.ipynb` or `5_Modeling_xTB.ipynb` from descriptor loading onward.
4. Use the revision notebooks for the specific reviewer-response analyses.
5. Use `6_Draw_Figures.ipynb` to regenerate manuscript and SI figures.

## Programmatic Access Example

```python
from ase.db import connect

db = connect("BorylXAT-DB.db")
ts_rows = list(db.select(category="ts"))
print(len(ts_rows))
print(ts_rows[0].key_value_pairs)
```

## Notes

- `main.py` is only a placeholder; the real workflow lives in the notebooks and `DFTStructureGenerator/`.
- The reviewer-runnable path is intentionally separated from the full Gaussian/xTB provenance path.
- `Data/csvs/Benchmark_Result.csv` was retired; the current benchmark table is `Data/csvs/Benchmark.csv`.
- QHARM revision analyses use `BorylXAT-DB_qh_update.db` and the QHARM descriptor files under `Data/descriptor/`.
