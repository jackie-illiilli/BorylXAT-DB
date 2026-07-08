# =============================================================================
# BorylXAT-DB: One-Click Structured Database Building Script (2025 Enhanced Version v2)
# Added: Storage of ΔG‡ / ΔG_rxn + bidirectional reaction mapping in ASE.db
#       + imaginary frequency / IRC / Hirshfeld charges / dipole moment / spin density / frontier orbital energy
# =============================================================================

import os
import re
import pandas as pd
import numpy as np
from ase import Atoms
from ase.db import connect

OUTPUT_ASE_DB = "BorylXAT-DB.db"
OUTPUT_PARQUET = "BorylXAT-DB.parquet"

AMINE_LB_IDS = frozenset(range(0, 69))
PHOSPHINE_LB_IDS = frozenset(range(69, 144))
NHC_LB_IDS = frozenset(range(144, 234))
ARYL_N_LB_IDS = frozenset(range(234, 388))

LB_TYPE_ORDER = ["Amine/Aryl N", "Phosphine", "NHC"]
B_TYPE_ORDER = ["BR3", "R2BH", "RBH2", "BH3"]
CL_SUBSTRATE_TYPE_ORDER = ["CCl4", "CCl3", "CCl2", "CCl"]
CL_CARBON_HYBRIDIZATION_TYPE_ORDER = ["C(sp3)-Cl", "C(sp2)-Cl", "C(sp)-Cl"]
CL_SUBSTRATE_SMARTS = [
    ("CCl4", "Cl*(Cl)(Cl)Cl"),
    ("CCl3", "Cl*(Cl)Cl"),
    ("CCl2", "Cl*Cl"),
]

# Regular patterns (adapted to the latest Cl_xxxxx_r)
PATTERNS = {
    'B':     re.compile(r'^B_(\d{5})$'),
    'LB':     re.compile(r'^LB_(\d{5})$'),
    'Cl':  re.compile(r'^Cl_(\d{5})_r$'),        # reactant substrate
    'complex_r':  re.compile(r'^B_(\d{5})_LB_(\d{5})_r$'),
    'complex_p':  re.compile(r'^B_(\d{5})_LB_(\d{5})_p$'),
    'ts':         re.compile(r'^B_(\d{5})_LB_(\d{5})_Cl_(\d{5})$'),
    'c_radical':  re.compile(r'^Cl_(\d{5})_p$'),
}

def classify_key(key: str):
    for pat, regex in PATTERNS.items():
        m = regex.match(key)
        if m:
            groups = [int(g) for g in m.groups()]
            if pat == 'B':     return 'B',    (groups[0], np.nan, np.nan)
            if pat == 'LB':     return 'LB',    (np.nan, groups[0], np.nan)
            if pat == 'Cl':  return 'Cl', (np.nan, np.nan, groups[0])
            if pat == 'complex_r':  return 'complex_r', (groups[0], groups[1], np.nan)
            if pat == 'complex_p':  return 'complex_p', (groups[0], groups[1], np.nan)
            if pat == 'ts':         return 'ts',        (groups[0], groups[1], groups[2])
            if pat == 'c_radical':  return 'c_radical', (np.nan, np.nan, groups[0])
    return 'unknown', None


def _get_chem():
    from rdkit import Chem

    return Chem


def _as_db(db_or_path):
    if hasattr(db_or_path, "select"):
        return db_or_path
    return connect(db_or_path)


def get_lb_type(lb_id):
    """Classify Lewis base ids using the dataset index ranges."""
    lb_id = int(lb_id)
    if lb_id in AMINE_LB_IDS or lb_id in ARYL_N_LB_IDS:
        return "Amine/Aryl N"
    if lb_id in PHOSPHINE_LB_IDS:
        return "Phosphine"
    if lb_id in NHC_LB_IDS:
        return "NHC"
    return "Other"


def get_boron_type(smiles):
    """Classify boron radical type by the number of hydrogens on boron."""
    Chem = _get_chem()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid boron radical SMILES: {smiles}")

    b_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "B"]
    if not b_atoms:
        raise ValueError(f"No boron atom found in SMILES: {smiles}")

    n_h = b_atoms[0].GetTotalNumHs()
    return {0: "BR3", 1: "R2BH", 2: "RBH2", 3: "BH3"}.get(n_h, f"BH{n_h}")


def get_cl_substrate_type(smiles):
    """Classify chloride substrates by the number of chlorine substituents."""
    Chem = _get_chem()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid chloride substrate SMILES: {smiles}")

    mol = Chem.AddHs(mol)
    for substrate_type, smarts in CL_SUBSTRATE_SMARTS:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
            return substrate_type
    return "CCl"


def get_cl_carbon_hybridization_type(smiles, cl_atom_id=None):
    """Classify a chloride substrate by the hybridization of the C atom bonded to Cl."""
    Chem = _get_chem()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid chloride substrate SMILES: {smiles}")

    if cl_atom_id is None:
        cl_atoms = [atom for atom in mol.GetAtoms() if atom.GetSymbol() == "Cl"]
        if not cl_atoms:
            raise ValueError(f"No chlorine atom found in SMILES: {smiles}")
        cl_atom = cl_atoms[0]
    else:
        cl_atom = mol.GetAtomWithIdx(int(cl_atom_id))
        if cl_atom.GetSymbol() != "Cl":
            raise ValueError(f"Atom {cl_atom_id} is not chlorine in SMILES: {smiles}")

    carbon_neighbors = [
        atom for atom in cl_atom.GetNeighbors()
        if atom.GetSymbol() == "C"
    ]
    if not carbon_neighbors:
        return "Other"

    hybridization = carbon_neighbors[0].GetHybridization()
    hybridization_name = str(hybridization).split(".")[-1].lower()
    if hybridization_name in {"sp3", "sp2", "sp"}:
        return f"C({hybridization_name})-Cl"
    return f"C({hybridization_name})-Cl"


def build_boron_type_map(db_or_path=OUTPUT_ASE_DB):
    """Return {B_id: B_type} for all boron radical rows in the ASE database."""
    db = _as_db(db_or_path)
    return {
        int(row.B_id): get_boron_type(row.smiles)
        for row in db.select(category="B")
    }


def build_lb_type_map(db_or_path=OUTPUT_ASE_DB):
    """Return {LB_id: LB_type} for all Lewis base rows in the ASE database."""
    db = _as_db(db_or_path)
    return {
        int(row.LB_id): get_lb_type(row.LB_id)
        for row in db.select(category="LB")
    }


def build_cl_substrate_type_map(db_or_path=OUTPUT_ASE_DB):
    """Return {Cl_id: Cl_type} for all chloride substrate rows in the ASE database."""
    db = _as_db(db_or_path)
    return {
        int(row.Cl_id): get_cl_substrate_type(row.smiles)
        for row in db.select(category="Cl")
    }


def add_reactant_type_columns(ts_df, db_or_path=OUTPUT_ASE_DB):
    """Add LB_type, B_type, and Cl_type columns to a TS dataframe."""
    typed_df = ts_df.copy()
    typed_df["LB_type"] = typed_df["LB_id"].map(build_lb_type_map(db_or_path))
    typed_df["B_type"] = typed_df["B_id"].map(build_boron_type_map(db_or_path))
    typed_df["Cl_type"] = typed_df["Cl_id"].map(build_cl_substrate_type_map(db_or_path))
    return typed_df


def build_ts_type_distribution_dataframe(db_or_path=OUTPUT_ASE_DB):
    """Build the TS dataframe used by the 8980 TS type-distribution plot."""
    db = _as_db(db_or_path)
    ts_df = pd.DataFrame([
        {
            "key": row.key,
            "B_id": int(row.B_id),
            "LB_id": int(row.LB_id),
            "Cl_id": int(row.Cl_id),
            "barrier_kcal": getattr(row, "barrier_kcal", np.nan),
            "delta_g_rxn_kcal": getattr(row, "delta_g_rxn_kcal", np.nan),
        }
        for row in db.select(category="ts")
    ])
    if ts_df.empty:
        return pd.DataFrame(
            columns=[
                "key",
                "B_id",
                "LB_id",
                "Cl_id",
                "barrier_kcal",
                "delta_g_rxn_kcal",
                "LB_type",
                "B_type",
                "Cl_type",
            ]
        )
    return add_reactant_type_columns(ts_df, db)


def add_bep_type_columns(df):
    """Add reactant-type columns used by BEP plots to a reaction dataframe."""
    typed_df = df.copy()
    typed_df["B_type"] = typed_df["B_smiles"].map(get_boron_type)
    typed_df["LB_type"] = typed_df["N_Index"].map(get_lb_type)
    typed_df["Cl_type"] = typed_df["Cl_smiles"].map(get_cl_substrate_type)
    typed_df["Cl_hybrid_type"] = typed_df.apply(
        lambda row: get_cl_carbon_hybridization_type(
            row["Cl_smiles"],
            row["Cl_Atomid"],
        ),
        axis=1,
    )
    return typed_df

def build_databases(data_dict: dict, rewrite: bool = True):
    """
    Each value in data_dict is a dict containing the following fields:
      Required: smiles, atoms, gibbs_hartree
      Optional (non-TS): hirshfeld_charges, dipole_moment, homo_energy_kcal, lumo_energy_kcal
      Optional (B/complex_r/c_radical): spin_densities
      Optional (TS): imaginary_frequency_cm_1, imaginary_freq_displacement,
                 irc_forward_positions, irc_reverse_positions
    """
    if rewrite and os.path.exists(OUTPUT_ASE_DB):
        os.remove(OUTPUT_ASE_DB)

    db = connect(OUTPUT_ASE_DB)
    records = []
    ts_mapping = {}
    reactant_to_ts = {}

    HARTREE_TO_KCAL = 627.509

    # ===================== First pass: Write all structures =====================
    for key, entry in data_dict.items():
        smiles = entry["smiles"]
        atoms = entry["atoms"]
        gibbs_hartree = entry["gibbs_hartree"]

        if atoms is None and smiles is None:
            continue

        category, ids = classify_key(key)
        if category == 'unknown':
            print(f"Warning: Unrecognized key: {key}")
            continue

        bid, lid, clid = ids
        natoms = len(atoms) if atoms else np.nan
        formula = atoms.get_chemical_formula() if atoms else None
        charge = int(atoms.get_initial_charges().sum()) if atoms else 0

        # ==================== Put searchable fields into kvp ====================
        kvp = {
            "category": category,
            "B_id": float(bid),
            "LB_id": float(lid),
            "Cl_id": float(clid),
            "gibbs_hartree": float(gibbs_hartree),
            "charge": charge,
            "temperature_K": 298.15,
            "solvent": "toluene",
            "smiles": smiles if isinstance(smiles, str) else None,
            "source_key": key,
        }

        # Scalar fields common to non-TS
        if category != 'ts':
            if entry.get("dipole_moment") is not None:
                kvp["dipole_moment_debye"] = float(entry["dipole_moment"])
            if entry.get("homo_energy_kcal") is not None:
                kvp["homo_energy_kcal"] = float(entry["homo_energy_kcal"])
            if entry.get("lumo_energy_kcal") is not None:
                kvp["lumo_energy_kcal"] = float(entry["lumo_energy_kcal"])

        # TS exclusive markers (searchable)
        if category == 'ts':
            kvp.update({
                "is_transition_state": True,
                "expected_imaginary_freqs": 1,
            })
            if entry.get("imaginary_frequency_cm_1") is not None:
                kvp["imaginary_frequency_cm_1"] = float(entry["imaginary_frequency_cm_1"])

        # ==================== Put complex/non-scalar fields into data ====================
        data_dict_for_ase = {
            "formula": formula,
        }

        # Non-TS: Hirshfeld charges (List[float])
        if category != 'ts' and entry.get("hirshfeld_charges") is not None:
            data_dict_for_ase["hirshfeld_charges"] = entry["hirshfeld_charges"]

        # B / complex_r / c_radical: Spin densities (List[float])
        if category in ('B', 'complex_r', 'c_radical') and entry.get("spin_densities") is not None:
            data_dict_for_ase["spin_densities"] = entry["spin_densities"]

        # TS: Imaginary frequency displacement vector + IRC endpoints
        if category == 'ts':
            if entry.get("imaginary_freq_displacement") is not None:
                data_dict_for_ase["imaginary_freq_displacement"] = entry["imaginary_freq_displacement"]
            if entry.get("irc_forward_positions") is not None:
                data_dict_for_ase["irc_forward_positions"] = entry["irc_forward_positions"]
            if entry.get("irc_reverse_positions") is not None:
                data_dict_for_ase["irc_reverse_positions"] = entry["irc_reverse_positions"]

        # Write: kvp for search, data for storing complex information
        db.write(atoms, key=key, key_value_pairs=kvp, data=data_dict_for_ase)

        # Parquet records
        records.append({
            "key": key,
            "category": category,
            "B_id": bid,
            "LB_id": lid,
            "Cl_id": clid,
            "smiles": smiles,
            "natoms": natoms,
            "formula": formula,
            "charge": charge,
            "gibbs_hartree": float(gibbs_hartree),
            "numbers": atoms.numbers.tolist() if atoms else None,
            "positions": atoms.positions.tolist() if atoms else None,
            # Added descriptor columns
            "hirshfeld_charges": entry.get("hirshfeld_charges"),
            "dipole_moment_debye": entry.get("dipole_moment"),
            "spin_densities": entry.get("spin_densities"),
            "homo_energy_kcal": entry.get("homo_energy_kcal"),
            "lumo_energy_kcal": entry.get("lumo_energy_kcal"),
            "imaginary_frequency_cm_1": entry.get("imaginary_frequency_cm_1"),
            "imaginary_freq_displacement": entry.get("imaginary_freq_displacement"),
            "irc_forward_positions": entry.get("irc_forward_positions"),
            "irc_reverse_positions": entry.get("irc_reverse_positions"),
        })

        # Mapping collection
        if category == 'ts':
            b5 = f"{int(bid):05}"
            l5 = f"{int(lid):05}"
            c5 = f"{int(clid):05}"
            ts_mapping[key] = {
                "reactant_complex": f"B_{b5}_LB_{l5}_r",
                "reactant_cl": f"Cl_{c5}_r",
                "product_complex": f"B_{b5}_LB_{l5}_p",
                "product_c_radical": f"Cl_{c5}_p",
            }
            reactant_to_ts.setdefault(f"B_{b5}_LB_{l5}_r", []).append(key)
            reactant_to_ts.setdefault(f"Cl_{c5}_r", []).append(key)

    # ===================== Second pass: Calculate barriers and UPDATE kvp =====================
    df = pd.DataFrame(records)
    updated_rows = []
    for ts_key, mapping in ts_mapping.items():
        r_complex_key = mapping["reactant_complex"]
        r_sub_key = mapping["reactant_cl"]
        p_complex_key = mapping["product_complex"]
        p_rad_key = mapping["product_c_radical"]

        row_ts = db.get(key=ts_key)  # Getting by key is more stable
        g_ts = row_ts.gibbs_hartree

        barrier = np.nan
        delta_g = np.nan

        if r_complex_key in df['key'].values and r_sub_key in df['key'].values:
            g_r_complex = df[df.key == r_complex_key].gibbs_hartree.iloc[0]
            g_r_sub = df[df.key == r_sub_key].gibbs_hartree.iloc[0]
            barrier = (g_ts - g_r_complex - g_r_sub) * HARTREE_TO_KCAL

        if all(k in df['key'].values for k in [p_complex_key, p_rad_key, r_complex_key, r_sub_key]):
            g_p_complex = df[df.key == p_complex_key].gibbs_hartree.iloc[0]
            g_p_rad = df[df.key == p_rad_key].gibbs_hartree.iloc[0]
            delta_g = (g_p_complex + g_p_rad - g_r_complex - g_r_sub) * HARTREE_TO_KCAL


        db.update(row_ts.id,
                  barrier_kcal=barrier,
                  delta_g_rxn_kcal=delta_g,
                  reactant_complex_key=r_complex_key,
                  reactant_cl_key=r_sub_key,
                  product_complex_key=p_complex_key,
                  product_c_radical_key=p_rad_key)
        # Update Parquet records
        updated_rows.append({
            "key": ts_key,
            "barrier_kcal": barrier,
            "delta_g_rxn_kcal": delta_g,
            **mapping
        })

    # Update associated_ts_keys of reactant
    for reactant_key, ts_list in reactant_to_ts.items():
        row = db.get(key=reactant_key)
        db.update(row.id, data={"associated_ts_keys": ts_list})  # Put list in data

    # ===================== Merge to Parquet =====================
    update_df = pd.DataFrame(updated_rows)
    final_df = df.merge(update_df, on="key", how="left")

    final_df.to_parquet(OUTPUT_PARQUET, index=False, compression="zstd")

    print(f"Database building completed!")
    print(f"  ASE Database       : {OUTPUT_ASE_DB}    ({db.count()} records)")
    print(f"  Parquet Dataset    : {OUTPUT_PARQUET} ({len(final_df)} rows)")
    print(f"  Number of TS with successfully calculated barriers: {final_df.barrier_kcal.notna().sum()}")
