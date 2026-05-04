# Borane Radical C-Cl Activation Database

[中文说明](docs/zh-CN/README.md) | [Database Structure](Database_Structure.md) | [数据库结构说明](docs/zh-CN/Database_Structure.md)

This repository contains the workflow used to build, curate, analyze, and model a quantum-chemistry dataset for Lewis base activated borane radical mediated C-Cl bond activation. It is the code companion for the manuscript project on large-scale transition-state data generation, mechanism analysis, and machine-learning assisted reactivity prediction.

The project combines:

- reaction-space enumeration from borane, Lewis base, and chloride substrate libraries
- RDKit and xTB based conformer generation
- Gaussian job generation for ground states, constrained optimizations, TS searches, SPE corrections, and IRC validation
- structured database construction in ASE SQLite and Parquet formats
- descriptor extraction, statistical analysis, benchmarking, and ML modeling

## Project Scope

According to the manuscript outline, the study targets a reaction space built from:

- 55 borane radicals
- 386 Lewis bases in the manuscript description
- 179 chloride substrates

This combinatorial space is filtered before TS calculation. The current generated database files are:

- `boron_ccl.db`: ASE SQLite database with `50057` structures
- `boron_ccl_dataset.parquet`: flattened Parquet dataset

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

1. Which borane/Lewis base/substrate combinations are thermodynamically feasible?
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
| `3_Build_DataBase.ipynb` | Consolidates parsed outputs into `boron_ccl.db` and `boron_ccl_dataset.parquet`, and prepares analysis figures |
| `4_Benchmark.ipynb` | DFT method benchmarking against a smaller reference set |
| `5_draw_molecule.ipynb` | Molecule drawing helpers for figures and SI |
| `6_Modeling.ipynb` | Descriptor generation and CatBoost modeling for reactivity prediction and OOD analysis |

`main.py` is only a placeholder; the real workflow lives in the notebooks and the `DFTStructureGenerator` package.

## Core Package Structure

`DFTStructureGenerator/` contains the reusable workflow code:

- `B_N_Cl.py`: reaction-site detection, reactant combination generation, DFT job preparation, TS structure generation, IRC job handling, TS summary export, and reaction AAM SMILES annotation
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
6_Modeling.ipynb
```

If you only want the final structured dataset, start from the existing files:

- `boron_ccl.db`
- `boron_ccl_dataset.parquet`
- `Data/descriptor/*.pkl`

For programmatic access to the ASE database:

```python
from ase.db import connect

db = connect("boron_ccl.db")
ts_rows = list(db.select(category="ts"))
print(len(ts_rows))
print(ts_rows[0].key_value_pairs)
```

## Notes for Reuse

- Naming conventions such as `B_00001`, `LB_00001`, `Cl_00001_r`, and `B_00001_LB_00001_Cl_00001` are central to the whole workflow.
- Many scripts assume the historical folder layout used for HPC calculations, so path adjustments may be needed before full reruns.
- The repository contains generated artifacts; reproducing every calculation from scratch requires the original Gaussian/xTB runtime environment.
- The manuscript-facing Lewis base count is 386, while the database retains 387 standalone `LB` entries because `LB_00623` is kept for provenance but is absent from stable B-LB complexes and TS entries.

## Repository Status

This repository already serves well as:

- a data-packaged computational chemistry project
- a reproducible record of the borane radical C-Cl activation workflow
- a starting point for mechanism analysis and descriptor-based ML on main-group radical reactions

Potential future cleanup items include expanding `pyproject.toml`, adding a frozen environment file, and exposing the notebook workflow through a small CLI.
