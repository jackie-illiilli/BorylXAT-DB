# =============================================================================
# Boron Radical Catalyzed C-Cl Activation Dataset: One-Click Structured Database Building Script (2025 Enhanced Version v2)
# Added: Storage of ΔG‡ / ΔG_rxn + bidirectional reaction mapping in ASE.db
#       + imaginary frequency / IRC / Hirshfeld charges / dipole moment / spin density / frontier orbital energy
# =============================================================================

import os
import re
import pandas as pd
import numpy as np
from ase import Atoms
from ase.db import connect

OUTPUT_ASE_DB = "boron_ccl.db"
OUTPUT_PARQUET = "boron_ccl_dataset.parquet"

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
