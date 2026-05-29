# BorylXAT-DB: a transition-state database for boryl-radical-mediated C–Cl atom transfer

[中文说明](docs/zh-CN/README.md) | [Database Structure](Database_Structure.md) | [数据库结构说明](docs/zh-CN/Database_Structure.md)

This repository contains the workflow used to build, curate, analyze, and model BorylXAT-DB, a quantum-chemistry dataset for Lewis base-coordinated boryl-radical-mediated C–Cl atom transfer. It is the code companion for the manuscript project on large-scale transition-state data generation, mechanism analysis, and machine-learning assisted reactivity prediction.

Public code repository: <https://github.com/jackie-illiilli/BorylXAT-DB>. Large database files are distributed through Zenodo: <https://doi.org/10.5281/zenodo.20134535>.

The project combines:

- reaction-space enumeration from boryl radical, Lewis base, and chlorinated substrate libraries
- RDKit and xTB based conformer generation
- Gaussian job generation for ground states, constrained optimizations, TS searches, SPE corrections, and IRC validation
- structured database construction in ASE SQLite and Parquet formats
- descriptor extraction, statistical analysis, benchmarking, and ML modeling

## Project Scope

According to the manuscript outline, the study targets a reaction space built from:

- 55 boryl radicals
- 386 Lewis bases in the manuscript description
- 179 chloride substrates

This combinatorial space is filtered before TS calculation. The current generated database files are:

- `BorylXAT-DB.db`: ASE SQLite database with `50057` structures
- `BorylXAT-DB.parquet`: flattened Parquet dataset

The current ASE database contains:

| Category | Count |
| --- | ---: |
| `B` | 55 |
| `LB` | 387 |
| `Cl` | 179 |
| `complex_r` | 20010 |
| `complex_p` | 20010 |
| `c_radical` | 179 |
| `ts` | 9237 |

The manuscript-facing reaction space counts 386 Lewis bases because one standalone Lewis base entry (`LB_00623`) does not form thermodynamically stable B-LB complexes with any borane radical. The database intentionally retains this molecule in the `LB` category for provenance, but it does not appear in the filtered B-LB complex set or TS reaction entries.

## Scientific Goal

The codebase is built to answer three linked questions:

1. Which boryl radical/Lewis base/substrate combinations are thermodynamically feasible?
2. What are the geometric and energetic characteristics of the validated halogen-atom-transfer transition states?
3. Can these DFT results be converted into descriptors and predictive ML models for `ΔG‡` and related trends?

## Workflow Overview

The practical workflow in this repository is notebook driven:

1. Enumerate reactants and reaction sites from curated reactant tables.
2. Generate initial 3D structures with RDKit and optimize conformers with xTB/CREST.
3. Build ground-state borane, Lewis base, chloride, and borane-Lewis base complex structures.
4. Generate TS guesses by combining product-like B-LB-Cl geometries with chloride reactants.
5. Run constrained optimization, TS optimization, SPE correction, and IRC validation.
6. Parse Gaussian outputs and collect energies, charges, spin densities, and geometries.
7. Assemble the final ASE/Parquet databases.
8. Build descriptors, benchmark methods, and train ML models.

## Notebook Guide

The main entry points are the notebooks in repository root:

| Notebook | Role |
| --- | --- |
| `1_Calc_Reactant.ipynb` | Reactant preparation, reaction-site enumeration, xTB conformer sampling, and DFT setup for boranes, Lewis bases, chlorides, and B-LB complexes |
| `2_Calc_TS.ipynb` | TS guess generation, constrained optimization, TS search, SPE correction, IRC analysis, and TS summary generation |
| `3_Build_DataBase.ipynb` | Consolidates parsed outputs into `BorylXAT-DB.db` and `BorylXAT-DB.parquet`; includes reviewer-runnable database inspection and Figure 1 style summary plots |
| `4_Benchmark.ipynb` | DFT method benchmark workflow with separated raw Gaussian-dependent sections and reviewer-runnable result analysis |
| `5_Modeling.ipynb` | Descriptor generation, precomputed-descriptor loading, CatBoost modeling, validation plots, and OOD analysis |
| `6_Draw_Figures.ipynb` | Consolidated manuscript/SI figure-generation notebook; collects figure code previously spread across database-building and molecule-drawing workflows |

`main.py` is only a placeholder; the real workflow lives in the notebooks and the `DFTStructureGenerator` package.

## Notebook Curation and Review Readiness

The notebooks have been curated for manuscript review, with computational-provenance cells kept visible and reviewer-runnable analysis separated where possible.

| Notebook | Curation performed | Reviewer usability |
| --- | --- | --- |
| `1_Calc_Reactant.ipynb` | Rewrote outdated markdown notes, tightened comments around reactant enumeration/optimization, removed unused imports, and removed stale working variables. | Documents the production reactant-optimization workflow. A full rerun requires Gaussian, xTB/CREST, and the historical calculation workspace. |
| `2_Calc_TS.ipynb` | Removed unused imports/variables, moved `complete_target()` into `DFTStructureGenerator/borane_xat_workflow.py`, and moved the target-CSV completion block to the final optional supplement/export section. | Documents TS generation, constrained optimization, TS optimization, SPE, IRC, and summary export. A full rerun requires the original Gaussian/xTB folders. |
| `3_Build_DataBase.ipynb` | Added review tags that distinguish raw Gaussian parsing from database inspection. Figure-producing logic that belongs to the manuscript figure set is now centralized in `6_Draw_Figures.ipynb`. | Reviewers can inspect the released `BorylXAT-DB.db` and `BorylXAT-DB.parquet` without rebuilding from raw logs. |
| `4_Benchmark.ipynb` | Cleaned benchmark comments, grouped setup/method definitions/result analysis, and marked raw structure/input collection separately from checked-in result analysis. | Reviewers can load `Data/csvs/Benchmark_Result.csv` and regenerate the benchmark summary figure without the original Gaussian folders. |
| `5_Modeling.ipynb` | Polished imports and comments, separated descriptor generation from descriptor loading, marked descriptor generation as optional, and added the missing Figure 6B `plt.savefig(...)` call. | Reviewers can use pre-extracted descriptors from `Data/descriptor/` to skip the expensive descriptor-building stage, then rerun model validation and plotting. |
| `6_Draw_Figures.ipynb` | Consolidated manuscript/SI figure-generation code from the database-building and molecule-drawing workflows into one notebook. Former Figure 4 labels and output filenames were corrected to Figure 5. | This is the current reviewer-facing entry point for regenerating manuscript and SI figures from checked-in data products. |

The historical module name `B_N_Cl` is better represented as `borane_xat_workflow`: the workflow now covers borane radical, Lewis-base complex, chloride substrate, XAT transition-state generation, summary export, and reaction annotation rather than a static B/N/Cl data container.

The review tags used in notebooks are:

| Tag | Meaning |
| --- | --- |
| `[REVIEWER-RUNNABLE]` | Cells that should run from the checked-in data files without the original HPC/Gaussian working folders |
| `[RAW-GAUSSIAN/E:/work]` | Provenance or reproduction cells that depend on the historical Gaussian/xTB output tree, especially paths under `E:/work` |
| `[OPTIONAL-DESCRIPTOR-GENERATION]` | Expensive descriptor-generation cells; the repository includes precomputed descriptor files for faster review |

## Repository Structure

```text
.
|-- DFTStructureGenerator/     Reusable Python modules for structure generation, parsing, database building, descriptors, and plotting
|-- Data/
|   |-- ChemDraw/              ChemDraw source files for the reactant libraries
|   |-- csvs/                  Intermediate screening tables before TS calculation, plus reactant SMILES included in the database
|   |-- descriptor/            Saved descriptor files used by the modeling workflow
|   `-- TS/                    The 9237 transition-state records and corresponding TS coordinates
|-- Figure/                    Original manuscript and SI figure outputs
|-- 1_Calc_Reactant.ipynb      Reactant preparation workflow
|-- 2_Calc_TS.ipynb            Transition-state generation and validation workflow
|-- 3_Build_DataBase.ipynb     Database construction and inspection workflow
|-- 4_Benchmark.ipynb          DFT benchmark workflow
|-- 5_Modeling.ipynb           Descriptor-based machine-learning workflow
|-- 6_Draw_Figures.ipynb       Manuscript and SI figure-generation workflow
`-- Database_Structure.md      Field-level description of the released database
```

## Core Package Structure

`DFTStructureGenerator/` contains the reusable workflow code:

- `borane_xat_workflow.py`: reaction-site detection, reactant combination generation, DFT job preparation, TS structure generation, IRC job handling, TS summary export, and reaction AAM SMILES annotation
- `Build_DataBase.py`: converts parsed structure dictionaries into ASE and Parquet databases, and computes `ΔG‡` / `ΔG_rxn`
- `descriptor.py`: builds borane-Lewis base and chloride descriptor maps and converts reaction tables into ML-ready feature matrices
- `logfile_process.py`: Gaussian log parser for energies, structures, imaginary frequencies, IRC trajectories, charge/multiplicity, and failure handling
- `xtb_process.py`: xTB/CREST conformer workflow and PBS script generation
- `mol_manipulation.py`: geometry transforms, TS/IRC input generation, and Gaussian error-recovery helpers
- `FormatConverter.py`: converts between RDKit molecules, XYZ, Gaussian input, and charge tables
- `draw.py`: plotting utilities for scatter plots, distributions, metrics, and correlation maps
- `Tool.py`: small geometry and utility helpers

## Key Data Products

Besides the main database files, the repository already includes several useful intermediate products:

- `Data/csvs/reactants_B.csv`: borane entries (`55`)
- `Data/csvs/reactants_N.csv`: Lewis base entries and reaction sites (`415` rows before final filtering/merging)
- `Data/csvs/reactants_Cl.csv`: chloride substrate entries (`179`)
- `Data/csvs/reactants_B_N.csv`: filtered borane-Lewis base combinations (`20010`)
- `Data/csvs/reactants_B_N_full.csv`: larger combination table before the final subset (`22825`)
- `Data/csvs/Benchmark_Result.csv`: benchmark comparison table
- `Data/descriptor/*.pkl`: serialized descriptor maps used by the modeling notebook
- `Figure/`: manuscript and SI figures exported from the analysis notebooks

Database field details are described in [Database_Structure.md](Database_Structure.md). A Chinese version is available at [docs/zh-CN/Database_Structure.md](docs/zh-CN/Database_Structure.md).

## Figure Outputs

The consolidated figure notebook `6_Draw_Figures.ipynb` is the current entry point for manuscript and SI figure export. During cleanup, the former Figure 4 labels and filenames in this workflow were updated to Figure 5. The checked-in and notebook-generated outputs include:

- Figure 1 and Figure 3 summary plots in `Figure/`
- Figure 5 files: `Figure5_reaction_landscape.png`, `Figure5A_TS_type_distribution.png`, `Figure5B_BCl_BDE_by_B_type.png`, `Figure5B_BCl_BDFE_by_LB_type.png`, `Figure5B_CCl_BDFE_by_hybridization.png`, `Figure5C_BEP_residual.png`, and `Figure5D_TS_geometry.png`
- Benchmark output: `FigureS17_Benchmark_MAE_R2_combined.png`
- Modeling outputs exported when `5_Modeling.ipynb` is rerun: `Figure6B_model_validation.png`, `Figure6C_model_feature_importance.png`, and `FigureS21_descriptor_correlation_map.png`
- SI molecule grids and distribution plots under `Figure/FigureInSI/`

The figure-cleanup pass also renamed the generated Figure 4 image outputs to Figure 5 names, so the current `Figure/` directory matches the manuscript numbering used by `6_Draw_Figures.ipynb`.

## Electronic Structure Setup

From the current notebooks, the main Gaussian settings are:

- geometry optimization and frequencies: `B3LYP/6-31G(d) + D3BJ + SMD(toluene)`
- single-point correction: `wB97X-D/6-311+G(d,p) + SMD(toluene)`
- TS workflow: constrained optimization -> TS optimization -> IRC verification

The code also distinguishes open-shell and closed-shell species through spin multiplicity settings, and can export wavefunction-enabled SPE jobs for downstream charge analysis.

## Environment Requirements

The checked-in `pyproject.toml` is minimal and does not yet list the full scientific stack. Based on the current code, a working environment should include at least:

- Python 3.12+
- `numpy`, `pandas`, `scipy`, `matplotlib`, `tqdm`
- `rdkit`
- `ase`
- `catboost`
- `scikit-learn`
- `seaborn`
- `morfeus`
- `ipykernel`
- `openpyxl` and `pyarrow` for Excel/Parquet workflows

External software expected by the workflow:

- Gaussian
- xTB / CREST
- a PBS-style HPC environment for the generated batch scripts

## Typical Usage

If you want to rerun the full pipeline, the intended order is:

```text
1_Calc_Reactant.ipynb
2_Calc_TS.ipynb
3_Build_DataBase.ipynb
4_Benchmark.ipynb
5_Modeling.ipynb
6_Draw_Figures.ipynb
```

For reviewer-style reuse without the original Gaussian folders, the recommended path is:

1. Use `3_Build_DataBase.ipynb` only for checked-in database inspection and summary analysis cells.
2. Use `4_Benchmark.ipynb` from the benchmark-result loading section onward.
3. Use `5_Modeling.ipynb` from the precomputed descriptor-loading section onward.
4. Use `6_Draw_Figures.ipynb` to regenerate manuscript and SI figures from the checked-in data products.

If you only want the final structured dataset, start from the existing files:

- `BorylXAT-DB.db`
- `BorylXAT-DB.parquet`
- `Data/descriptor/*.pkl`

If the two database files are not present after cloning the code repository, download them from the Zenodo dataset record and place them in the repository root:

```bash
curl -L -o BorylXAT-DB.db "https://zenodo.org/records/20134535/files/BorylXAT-DB.db?download=1"
curl -L -o BorylXAT-DB.parquet "https://zenodo.org/records/20134535/files/BorylXAT-DB.parquet?download=1"
```

On Windows PowerShell:

```powershell
Invoke-WebRequest -Uri "https://zenodo.org/records/20134535/files/BorylXAT-DB.db?download=1" -OutFile "BorylXAT-DB.db"
Invoke-WebRequest -Uri "https://zenodo.org/records/20134535/files/BorylXAT-DB.parquet?download=1" -OutFile "BorylXAT-DB.parquet"
```

For programmatic access to the ASE database:

```python
from ase.db import connect

db = connect("BorylXAT-DB.db")
ts_rows = list(db.select(category="ts"))
print(len(ts_rows))
print(ts_rows[0].key_value_pairs)
```

## Code Availability for Review

The repository is organized so that reviewers can inspect the released dataset and regenerate the main analysis/figure outputs without access to the original private calculation tree.

Reviewer-runnable entry points:

- `3_Build_DataBase.ipynb`: database inspection and summary analysis using the checked-in ASE/Parquet files
- `4_Benchmark.ipynb`: benchmark method analysis from `Data/csvs/Benchmark_Result.csv`
- `5_Modeling.ipynb`: descriptor loading, CatBoost model validation, feature importance, OOD analysis, and Figure 6 exports from checked-in descriptor files
- `6_Draw_Figures.ipynb`: manuscript/SI figure regeneration from checked-in CSV/database products

Full provenance/rerun entry points:

- `1_Calc_Reactant.ipynb`: reactant enumeration, xTB conformer-search setup, Gaussian optimization/SPE setup, and reactant-table parsing
- `2_Calc_TS.ipynb`: TS guess generation, constrained optimization, TS search, SPE correction, IRC validation, and TS CSV completion
- raw sections of `3_Build_DataBase.ipynb` and `4_Benchmark.ipynb`: parsing original Gaussian logs and collecting benchmark structures/results

The full provenance route requires Gaussian, xTB/CREST, a PBS-style HPC environment, and path adjustment for the historical `E:/work` calculation folders. For standard review, use sections marked `[REVIEWER-RUNNABLE]` and skip sections marked `[RAW-GAUSSIAN/E:/work]`. In `5_Modeling.ipynb`, skip `[OPTIONAL-DESCRIPTOR-GENERATION]` unless descriptor regeneration is specifically required; the precomputed descriptor maps in `Data/descriptor/` are intended to save review time.

## Notes for Reuse

- Naming conventions such as `B_00001`, `LB_00001`, `Cl_00001_r`, and `B_00001_LB_00001_Cl_00001` are central to the whole workflow.
- Many scripts assume the historical folder layout used for HPC calculations, so path adjustments may be needed before full reruns.
- The repository contains generated artifacts; reproducing every calculation from scratch requires the original Gaussian/xTB runtime environment.
- For manuscript review, prefer cells marked `[REVIEWER-RUNNABLE]` and use the checked-in database, benchmark table, and descriptor pickle files unless full Gaussian-level provenance regeneration is needed.
- The manuscript-facing Lewis base count is 386, while the database retains 387 standalone `LB` entries because `LB_00623` is kept for provenance but is absent from stable B-LB complexes and TS entries.

## Repository Status

This repository already serves well as:

- a data-packaged computational chemistry project
- a reproducible record of the boryl-radical-mediated C–Cl atom-transfer workflow
- a starting point for mechanism analysis and descriptor-based ML on main-group radical reactions

Potential future cleanup items include expanding `pyproject.toml`, adding a frozen environment file, and exposing the notebook workflow through a small CLI.
