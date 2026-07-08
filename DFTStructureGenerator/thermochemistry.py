from __future__ import annotations

import sqlite3
import re
from pathlib import Path

import pandas as pd

from DFTStructureGenerator.project_paths import REPO_ROOT, TS_DATA_DIR


RRHO_DB_FILENAME = "BorylXAT-DB.db"
QHARM_DB_FILENAME = "BorylXAT-DB_qh_update.db"
REACTION_KEY_COLUMNS = ["B_Index", "N_Index", "Cl_Index"]


def database_path(use_qharm: bool = True, repo_root: str | Path = REPO_ROOT) -> Path:
    """Return the selected thermochemistry database path (QHARM by default)."""
    filename = QHARM_DB_FILENAME if use_qharm else RRHO_DB_FILENAME
    return Path(repo_root) / filename


def _read_qharm_reactions(db_path: str | Path) -> pd.DataFrame:
    """Read the reaction-level QHARM values stored on transition-state rows."""
    db_path = Path(db_path)
    if not db_path.exists():
        raise FileNotFoundError(f"Thermochemistry database not found: {db_path}")

    query = """
    SELECT
        CAST(MAX(CASE WHEN key = 'B_id' THEN value END) AS INTEGER) AS B_Index,
        CAST(MAX(CASE WHEN key = 'LB_id' THEN value END) AS INTEGER) AS N_Index,
        CAST(MAX(CASE WHEN key = 'Cl_id' THEN value END) AS INTEGER) AS Cl_Index,
        MAX(CASE WHEN key = 'barrier_qharm_kcal' THEN value END) AS "deltaGa_qharm(kcal/mol)",
        MAX(CASE WHEN key = 'delta_g_rxn_qharm_kcal' THEN value END) AS "deltaG_qharm(kcal/mol)"
    FROM number_key_values
    WHERE key IN (
        'B_id', 'LB_id', 'Cl_id', 'is_transition_state',
        'barrier_qharm_kcal', 'delta_g_rxn_qharm_kcal'
    )
    GROUP BY id
    HAVING MAX(CASE WHEN key = 'is_transition_state' THEN value END) = 1
    """
    with sqlite3.connect(db_path) as connection:
        qharm = pd.read_sql_query(query, connection)

    value_columns = ["deltaGa_qharm(kcal/mol)", "deltaG_qharm(kcal/mol)"]
    if qharm.empty:
        raise ValueError(f"No transition-state QHARM records found in {db_path}")
    if qharm.duplicated(REACTION_KEY_COLUMNS).any():
        raise ValueError(f"Duplicate reaction keys found in {db_path}")
    if qharm[value_columns].isna().any().any():
        missing = qharm[value_columns].isna().sum().to_dict()
        raise ValueError(f"Incomplete QHARM values in {db_path}: {missing}")
    return qharm


def load_structure_thermochemistry(db_path: str | Path | None = None) -> pd.DataFrame:
    """Load structure keys, categories, and all stored Gibbs-energy variants."""
    selected_db = Path(db_path) if db_path is not None else database_path(True)
    if not selected_db.exists():
        raise FileNotFoundError(f"Thermochemistry database not found: {selected_db}")

    query = """
    SELECT systems.id, text_values.key, text_values.category,
           number_values.gibbs_hartree, number_values.gibbs_qharm_hartree,
           number_values.gibbs_rrho_hartree, number_values.gibbs_qrrho_hartree,
           number_values.gibbs_full_qrrho_hartree
    FROM systems
    LEFT JOIN (
        SELECT id,
               MAX(CASE WHEN key = 'key' THEN value END) AS key,
               MAX(CASE WHEN key = 'category' THEN value END) AS category
        FROM text_key_values
        WHERE key IN ('key', 'category')
        GROUP BY id
    ) AS text_values ON systems.id = text_values.id
    LEFT JOIN (
        SELECT id,
               MAX(CASE WHEN key = 'gibbs_hartree' THEN value END) AS gibbs_hartree,
               MAX(CASE WHEN key = 'gibbs_qharm_hartree' THEN value END) AS gibbs_qharm_hartree,
               MAX(CASE WHEN key = 'gibbs_rrho_aarontools_hartree' THEN value END) AS gibbs_rrho_hartree,
               MAX(CASE WHEN key = 'gibbs_qrrho_hartree' THEN value END) AS gibbs_qrrho_hartree,
               MAX(CASE WHEN key = 'gibbs_full_qrrho_hartree' THEN value END) AS gibbs_full_qrrho_hartree
        FROM number_key_values
        WHERE key IN (
            'gibbs_hartree', 'gibbs_qharm_hartree', 'gibbs_rrho_aarontools_hartree',
            'gibbs_qrrho_hartree', 'gibbs_full_qrrho_hartree'
        )
        GROUP BY id
    ) AS number_values ON systems.id = number_values.id
    """
    with sqlite3.connect(selected_db) as connection:
        thermo = pd.read_sql_query(query, connection)

    required = ["key", "category", "gibbs_qharm_hartree"]
    if thermo[required].isna().any().any():
        missing = thermo[required].isna().sum().to_dict()
        raise ValueError(f"Incomplete structure thermochemistry in {selected_db}: {missing}")
    thermo["status"] = "ok"
    return thermo


def load_reaction_dataset(
    csv_path: str | Path = TS_DATA_DIR / "Borane_all.csv",
    *,
    use_qharm: bool = True,
    db_path: str | Path | None = None,
    strict: bool = True,
) -> pd.DataFrame:
    """Load the curated reaction table with selectable thermochemistry.

    The public CSV remains the structural/metadata source.  With ``use_qharm=True``,
    reaction barriers and reaction free energies are joined from the QHARM database
    and replace the two historical column names used throughout the notebooks.
    The original RRHO values are retained in ``*_rrho(kcal/mol)`` columns.
    """
    csv_path = Path(csv_path)
    reactions = pd.read_csv(csv_path)
    reactions["deltaGa_rrho(kcal/mol)"] = reactions["deltaGa(kcal/mol)"]
    reactions["deltaG_rrho(kcal/mol)"] = reactions["deltaG(kcal/mol)"]

    selected_db = Path(db_path) if db_path is not None else database_path(use_qharm)
    if use_qharm:
        qharm_columns = ["deltaGa_qharm(kcal/mol)", "deltaG_qharm(kcal/mol)"]
        qharm = _read_qharm_reactions(selected_db)
        csv_has_qharm = all(column in reactions.columns for column in qharm_columns)
        if csv_has_qharm:
            check = reactions[REACTION_KEY_COLUMNS + qharm_columns].merge(
                qharm,
                on=REACTION_KEY_COLUMNS,
                how="left",
                validate="one_to_one",
                suffixes=("_csv", "_db"),
                indicator="_qharm_merge",
            )
            missing_mask = check["_qharm_merge"].ne("both")
            differences = {
                column: float(
                    (
                        check[f"{column}_csv"] - check[f"{column}_db"]
                    ).abs().max()
                )
                for column in qharm_columns
            }
            if strict and (missing_mask.any() or any(value > 1e-8 for value in differences.values())):
                raise ValueError(
                    f"CSV QHARM columns do not match {selected_db}: "
                    f"missing={int(missing_mask.sum())}, max_abs_difference={differences}"
                )
            source = "QHARM CSV columns (validated against database)"
        else:
            reactions = reactions.drop(
                columns=[column for column in qharm_columns if column in reactions.columns]
            )
            reactions = reactions.merge(
                qharm,
                on=REACTION_KEY_COLUMNS,
                how="left",
                validate="one_to_one",
                indicator="_qharm_merge",
            )
            missing_mask = reactions["_qharm_merge"].ne("both")
            if strict and missing_mask.any():
                examples = reactions.loc[missing_mask, REACTION_KEY_COLUMNS].head().to_dict("records")
                raise ValueError(
                    f"QHARM values missing for {int(missing_mask.sum())} CSV reactions; "
                    f"examples: {examples}"
                )
            reactions = reactions.drop(columns="_qharm_merge")
            source = "QHARM database columns"
        reactions["deltaGa(kcal/mol)"] = reactions["deltaGa_qharm(kcal/mol)"]
        reactions["deltaG(kcal/mol)"] = reactions["deltaG_qharm(kcal/mol)"]
    else:
        source = "original RRHO"

    reactions["thermochemistry_source"] = source
    reactions.attrs["thermochemistry_source"] = source
    reactions.attrs["database_path"] = str(selected_db)
    return reactions


def apply_selected_reaction_energy(
    features,
    reactions: pd.DataFrame,
    *,
    column: int = 0,
    reaction_energy_column: str = "deltaG(kcal/mol)",
):
    """Replace a descriptor matrix reaction-energy feature with the selected value."""
    import numpy as np

    selected = np.asarray(features, dtype=float).copy()
    if selected.shape[0] != len(reactions):
        raise ValueError(
            f"Feature/reaction row mismatch: {selected.shape[0]} != {len(reactions)}"
        )
    if reaction_energy_column not in reactions.columns:
        raise KeyError(f"Reaction-energy column not found: {reaction_energy_column}")
    selected[:, column] = reactions[reaction_energy_column].to_numpy(dtype=float)
    return selected


def apply_qharm_component_energies(
    bn_descriptors: dict,
    cl_descriptors: dict,
    *,
    use_qharm: bool = True,
    db_path: str | Path | None = None,
    strict: bool = True,
):
    """Update descriptor-map component free energies for full-space prediction.

    Descriptor entry zero stores a component reaction free energy in Hartree.  This
    function replaces it with the QHARM product-minus-reactant value while leaving
    every electronic/geometric descriptor unchanged.
    """
    bn_selected = {key: list(value) for key, value in bn_descriptors.items()}
    cl_selected = {key: list(value) for key, value in cl_descriptors.items()}
    if not use_qharm:
        return bn_selected, cl_selected

    thermo = load_structure_thermochemistry(db_path)
    gibbs = thermo.set_index("key")["gibbs_qharm_hartree"].to_dict()
    missing = []

    for descriptor_key, values in bn_selected.items():
        match = re.fullmatch(r"B_(\d+)_Nu_(\d+)", descriptor_key)
        if match is None:
            missing.append(descriptor_key)
            continue
        b_id, lb_id = (int(part) for part in match.groups())
        reactant = f"B_{b_id:05d}_LB_{lb_id:05d}_r"
        product = f"B_{b_id:05d}_LB_{lb_id:05d}_p"
        if reactant not in gibbs or product not in gibbs:
            missing.append(descriptor_key)
            continue
        values[0] = gibbs[product] - gibbs[reactant]

    for descriptor_key, values in cl_selected.items():
        match = re.fullmatch(r"Cl_(\d+)(?:_Claid_\d+)?", descriptor_key)
        if match is None:
            missing.append(descriptor_key)
            continue
        cl_id = int(match.group(1))
        reactant = f"Cl_{cl_id:05d}_r"
        product = f"Cl_{cl_id:05d}_p"
        if reactant not in gibbs or product not in gibbs:
            missing.append(descriptor_key)
            continue
        values[0] = gibbs[product] - gibbs[reactant]

    if strict and missing:
        raise KeyError(
            f"QHARM component energies missing for {len(missing)} descriptor keys; "
            f"examples: {missing[:5]}"
        )
    return bn_selected, cl_selected
