# Borane Radical Catalyzed C-Cl Activation Database Structure

[中文版本](docs/zh-CN/Database_Structure.md) | [Project README](README.md)

Based on the logic in `Build_DataBase.py`, the database-building workflow generates two core files for structured storage of DFT-derived 3D structures and thermodynamic data. The two formats serve different use cases:

1. **`boron_ccl2.db`** (ASE SQLite database): intended for use with ASE (Atomic Simulation Environment), enabling efficient storage of `Atoms` objects, searchable key-value metadata, and reaction mapping.
2. **`boron_ccl_dataset2.parquet`** (Parquet dataset): a flattened tabular dataset containing coordinates, atomic numbers, and descriptors, suitable for Pandas, PyArrow, and downstream machine-learning workflows.

---

## Included Chemical Categories

All structures are classified into seven categories according to the regex-based naming rule in `key`:

| Category | Description | Open/Closed Shell | Special Data |
| --- | --- | --- | --- |
| `B` | Borane radical catalyst | Open-shell | Spin density |
| `LB` | Lewis base / nucleophile | Closed-shell | — |
| `Cl` | Chloride substrate (`Cl_xxxxx_r`) | Closed-shell | — |
| `complex_r` | Reactant complex (B-LB reactant) | Open-shell | Spin density |
| `complex_p` | Product complex (B-LB product) | Closed-shell | — |
| `ts` | Transition state | — | Imaginary frequency / IRC / barrier |
| `c_radical` | Carbon-radical product (`Cl_xxxxx_p`) | Open-shell | Spin density |

---

## 1. ASE SQLite Database (`boron_ccl2.db`)

In addition to storing ASE `Atoms` objects, the ASE database contains:

- **key-value pairs (`key_value_pairs`)** for efficient filtering and querying
- **data blocks (`data`)** for complex list or dictionary values

### Key-Value Pairs (`key_value_pairs` / `kvp`)

**Fields shared by all species**

| Field | Type | Description |
| --- | --- | --- |
| `category` | str | Structure category |
| `B_id` | float | Borane index |
| `LB_id` | float | Lewis base index |
| `Cl_id` | float | Chloride substrate index |
| `gibbs_hartree` | float | Absolute Gibbs free energy in Hartree |
| `charge` | int | Net system charge |
| `temperature_K` | float | Thermodynamic temperature (`298.15`) |
| `solvent` | str | Implicit solvent (`toluene`) |
| `smiles` | str | SMILES or reaction AAM string |
| `source_key` | str | Original structure identifier |

**Additional fields for non-TS species**

| Field | Type | Description |
| --- | --- | --- |
| `dipole_moment_debye` | float | Dipole moment in Debye |
| `homo_energy_kcal` | float | HOMO energy in kcal/mol |
| `lumo_energy_kcal` | float | LUMO energy in kcal/mol |

**TS-specific fields**

| Field | Type | Description |
| --- | --- | --- |
| `is_transition_state` | bool | Always `True` for TS entries |
| `expected_imaginary_freqs` | int | Expected number of imaginary frequencies (`1`) |
| `imaginary_frequency_cm_1` | float | Imaginary frequency in cm^-1 |
| `barrier_kcal` | float | Activation free energy `ΔG‡` in kcal/mol |
| `delta_g_rxn_kcal` | float | Reaction free energy `ΔG_rxn` in kcal/mol |
| `reactant_complex_key` | str | Key of the reactant complex |
| `reactant_cl_key` | str | Key of the chloride substrate |
| `product_complex_key` | str | Key of the product complex |
| `product_c_radical_key` | str | Key of the carbon-radical product |

### Data Block (`data`)

**Fields shared by all species**

| Field | Type | Description |
| --- | --- | --- |
| `formula` | str | Molecular formula |

**Species other than B, LB, and TS (`Cl_r`, `c_radical`, `complex_r`, `complex_p`)**

| Field | Type | Description |
| --- | --- | --- |
| `hirshfeld_charges` | List[float] | Hirshfeld charges for all atoms |

**Open-shell species: B / `complex_r` / `c_radical`**

| Field | Type | Description |
| --- | --- | --- |
| `spin_densities` | List[float] | Mulliken spin densities for all atoms |

**TS entries**

| Field | Type | Description |
| --- | --- | --- |
| `imaginary_freq_displacement` | List[List[float]] | Imaginary-mode displacement vectors (`N x 3`) |
| `irc_forward_positions` | List[List[float]] | Forward IRC endpoint coordinates (`N x 3`) |
| `irc_reverse_positions` | List[List[float]] | Reverse IRC endpoint coordinates (`N x 3`) |

**Reactants (`complex_r` and `Cl_r`)**

| Field | Type | Description |
| --- | --- | --- |
| `associated_ts_keys` | List[str] | Keys of associated transition states |

---

## 2. Parquet Dataset (`boron_ccl_dataset2.parquet`)

The Parquet file stores a unified tabular `DataFrame` with the following groups of columns.

### Basic Information Columns

| Column | Type | Description |
| --- | --- | --- |
| `key` | str | Unique structure identifier |
| `category` | str | Structure category |
| `B_id`, `LB_id`, `Cl_id` | float | Component indices |
| `smiles` | str | SMILES or reaction AAM string |
| `natoms` | int | Number of atoms |
| `formula` | str | Molecular formula |
| `charge` | int | Total system charge |
| `gibbs_hartree` | float | Gibbs free energy in Hartree |

### 3D Structure Columns

| Column | Type | Description |
| --- | --- | --- |
| `numbers` | List[int] | Atomic numbers |
| `positions` | List[List[float]] | 3D coordinates (`N x 3`) |

### Electronic Descriptor Columns

| Column | Type | Applicable Category | Description |
| --- | --- | --- | --- |
| `hirshfeld_charges` | List[float] | Species other than B, LB, TS | Hirshfeld charges |
| `dipole_moment_debye` | float | All non-TS species | Dipole moment in Debye |
| `spin_densities` | List[float] | B / `complex_r` / `c_radical` | Mulliken spin densities |
| `homo_energy_kcal` | float | All non-TS species | HOMO energy in kcal/mol |
| `lumo_energy_kcal` | float | All non-TS species | LUMO energy in kcal/mol |

### TS-Specific Columns

| Column | Type | Description |
| --- | --- | --- |
| `imaginary_frequency_cm_1` | float | Imaginary frequency in cm^-1 |
| `imaginary_freq_displacement` | List[List[float]] | Imaginary-mode displacement vectors (`N x 3`) |
| `irc_forward_positions` | List[List[float]] | Forward IRC endpoint coordinates |
| `irc_reverse_positions` | List[List[float]] | Reverse IRC endpoint coordinates |
| `barrier_kcal` | float | Activation free energy `ΔG‡` in kcal/mol |
| `delta_g_rxn_kcal` | float | Reaction free energy `ΔG_rxn` in kcal/mol |
| `reactant_complex` | str | Key of the reactant complex |
| `reactant_cl` | str | Key of the chloride substrate |
| `product_complex` | str | Key of the product complex |
| `product_c_radical` | str | Key of the carbon-radical product |
