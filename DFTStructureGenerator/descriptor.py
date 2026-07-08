from __future__ import annotations

import os
import pickle

import numpy as np
import pandas as pd
from ase.db import connect
from morfeus import BuriedVolume
from tqdm import tqdm

from . import Tool, mol_manipulation
from .thermochemistry import apply_qharm_component_energies


DUPLICATE_N_IDS = [
    9, 43, 285, 310, 314, 345, 346, 347, 348, 349, 350, 351, 352, 353,
    354, 355, 356, 357, 358, 359, 360, 361, 362, 372, 375, 376,
]

DUPLICATE_CL_IDS = [
    624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636,
    637, 638, 639, 640, 642, 644, 645, 652, 653, 654, 655, 656, 657,
    658, 659, 660, 661, 662, 663, 664, 665, 666, 667, 668, 669, 670,
    671, 672, 673, 674, 675, 676, 677, 678, 679, 680, 681, 682, 683,
    684, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696,
    697, 698, 699, 700, 701, 702, 703, 704, 705, 706, 707, 708, 709,
    710, 711, 713, 714, 716, 717, 718, 719, 720, 721, 722,
]

# Backward-compatible aliases for notebook variables.
duplicate_N_id = DUPLICATE_N_IDS
duplicate_Cl_id = DUPLICATE_CL_IDS


def _sum_buried_volume(symbols, positions, center_atom_id, axis_atom_id, radius=6):
    bv = BuriedVolume(
        symbols,
        positions,
        center_atom_id + 1,
        include_hs=1,
        radius=radius,
        z_axis_atoms=[axis_atom_id + 1],
        excluded_atoms=[axis_atom_id + 1],
    )
    bv.octant_analysis()
    return float(np.sum(list(bv.octants["percent_buried_volume"].values())))


def _build_bn_descriptor(row, db):
    b_index = int(row["B_Index"])
    b_atomid = int(row["B_Atomid"])
    b_smiles = row["B_smiles"]
    n_index = int(row["N_Index"])
    n_atomid = int(row["N_Atomid"])
    react_eng = float(row["deltaG_react"])

    b_n_name = f"B_{b_index:05}_Nu_{n_index:05}"
    react_key = f"B_{b_index:05}_LB_{n_index:05}_r"
    prod_key = f"B_{b_index:05}_LB_{n_index:05}_p"

    b_mol = mol_manipulation.smiles2mol(b_smiles)
    b_atomnum = b_mol.GetNumAtoms()

    react_row = db.get(key=react_key)
    prod_row = db.get(key=prod_key)

    react_atoms = react_row.toatoms()
    prod_atoms = prod_row.toatoms()
    react_positions = react_atoms.get_positions()
    prod_positions = prod_atoms.get_positions()
    prod_symbols = prod_atoms.get_chemical_symbols()

    react_spin_densities = react_row.data.get("spin_densities")
    react_hirshfeld = react_row.data.get("hirshfeld_charges")
    prod_hirshfeld = prod_row.data.get("hirshfeld_charges")

    descriptor = [react_eng]
    descriptor += [react_spin_densities[b_atomid]]
    descriptor += [react_hirshfeld[b_atomid]]
    descriptor += [
        Tool.get_atoms_distance(
            react_positions[b_atomid],
            react_positions[n_atomid + b_atomnum - 1],
        )
    ]
    descriptor += [react_row.homo_energy_kcal]

    cl_index = prod_symbols.index("Cl")
    descriptor += [prod_hirshfeld[b_atomid], prod_hirshfeld[cl_index]]
    descriptor += [
        Tool.get_atoms_distance(
            prod_positions[b_atomid],
            prod_positions[cl_index],
        )
    ]
    descriptor += [prod_row.lumo_energy_kcal]
    descriptor += [
        _sum_buried_volume(
            prod_symbols,
            prod_positions,
            b_atomid,
            cl_index,
            radius=6,
        )
    ]
    return b_n_name, descriptor


def _build_cl_descriptor(row, db, duplicate_cl_ids):
    cl_index = int(row["Index"])
    cl_atomid = int(row["Atomid"])
    cl_smiles = row["Smiles"]
    react_eng = float(row["deltaG_react"])

    cl_mol = mol_manipulation.smiles2mol(cl_smiles)
    c_atom_idx = cl_mol.GetAtomWithIdx(cl_atomid).GetNeighbors()[0].GetIdx()

    react_key = f"Cl_{cl_index:05}_r"
    prod_key = f"Cl_{cl_index:05}_p"
    react_row = db.get(key=react_key)
    prod_row = db.get(key=prod_key)

    react_atoms = react_row.toatoms()
    react_positions = react_atoms.get_positions()
    react_symbols = react_atoms.get_chemical_symbols()

    react_hirshfeld = react_row.data.get("hirshfeld_charges")
    prod_spin_densities = prod_row.data.get("spin_densities")
    prod_hirshfeld = prod_row.data.get("hirshfeld_charges")

    descriptor = [react_eng]
    descriptor += [react_hirshfeld[cl_atomid], react_hirshfeld[c_atom_idx]]
    descriptor += [
        Tool.get_atoms_distance(
            react_positions[cl_atomid],
            react_positions[c_atom_idx],
        )
    ]
    descriptor += [react_row.lumo_energy_kcal]
    descriptor += [
        _sum_buried_volume(
            react_symbols,
            react_positions,
            c_atom_idx,
            cl_atomid,
            radius=6,
        )
    ]

    if c_atom_idx > cl_atomid:
        c_atom_idx -= 1
    descriptor += [prod_spin_densities[c_atom_idx]]
    descriptor += [prod_hirshfeld[c_atom_idx]]
    descriptor += [prod_row.homo_energy_kcal]

    cl_name = f"Cl_{cl_index:05}"
    if cl_index in duplicate_cl_ids:
        cl_name = f"Cl_{cl_index:05}_Claid_{cl_atomid:05}"
    return cl_name, descriptor


def build_descriptor_maps(
    db_path="BorylXAT-DB.db",
    bn_csv_path="data/csvs/reactants_B_N.csv",
    cl_csv_path="data/csvs/reactants_Cl.csv",
    duplicate_cl_ids=None,
    show_progress=True,
    use_qharm=False,
    load_bn_cl_map = [],
):
    """Build DFT descriptor maps with selectable component thermochemistry.

    When ``use_qharm`` is true, descriptor entry zero is computed from the
    QHARM Gibbs energies in ``db_path`` rather than the RRHO values stored in
    the reactant CSV tables.  All electronic and geometric entries are
    unchanged.
    """
    if duplicate_cl_ids is None:
        duplicate_cl_ids = DUPLICATE_CL_IDS

    db = connect(db_path)
    bn_csv = pd.read_csv(bn_csv_path)
    cl_csv = pd.read_csv(cl_csv_path)
    if len(load_bn_cl_map):
        bn_map, cl_map = load_descriptor_maps(
        bn_path=load_bn_cl_map[0],
        cl_path=load_bn_cl_map[1],
        )
    else:
        bn_map = {}
        bn_iter = bn_csv.iterrows()
        cl_iter = cl_csv.iterrows()
        if show_progress:
            bn_iter = tqdm(bn_iter, total=len(bn_csv))
            cl_iter = tqdm(cl_iter, total=len(cl_csv))

        for _, row in bn_iter:
            bn_name, descriptor = _build_bn_descriptor(row, db)
            if bn_name not in bn_map:
                bn_map[bn_name] = descriptor

        cl_map = {}
        for _, row in cl_iter:
            try:
                cl_name, descriptor = _build_cl_descriptor(row, db, duplicate_cl_ids)
                cl_map[cl_name] = descriptor
            except Exception as e:
                print(f"Error building descriptor for Cl_{row['Index']:05}: {e}")

    if use_qharm:
        bn_map, cl_map = apply_qharm_component_energies(
            bn_map,
            cl_map,
            use_qharm=True,
            db_path=db_path,
        )
    return bn_map, cl_map


def dataframe_to_descriptors(
    data_csv,
    bn_map,
    cl_map,
    duplicate_cl_ids=None,
    show_progress=True,
    reaction_energy_column=None,
):
    """Combine component maps into reaction descriptors.

    ``reaction_energy_column`` can explicitly supply the selected reaction
    free energy in kcal/mol.  This is the preferred QHARM path; leaving it as
    ``None`` preserves the historical component-map calculation.
    """
    if duplicate_cl_ids is None:
        duplicate_cl_ids = DUPLICATE_CL_IDS
    if reaction_energy_column is not None and reaction_energy_column not in data_csv.columns:
        raise KeyError(f"Reaction-energy column not found: {reaction_energy_column}")

    rows = data_csv.iterrows()
    if show_progress:
        rows = tqdm(rows, total=len(data_csv))

    all_xs = []
    for _, row in rows:
        b_index = int(row["B_Index"])
        n_index = int(row["N_Index"])
        cl_index = int(row["Cl_Index"])
        try:
            cl_atomid = int(row["Cl_Atomid"])
        except:
            cl_atomid = None

        b_n_name = f"B_{b_index:05}_Nu_{n_index:05}"
        cl_name = f"Cl_{cl_index:05}"
        if cl_index in duplicate_cl_ids:
            cl_name = f"Cl_{cl_index:05}_Claid_{cl_atomid:05}"

        des_a = bn_map[b_n_name]
        des_b = cl_map[cl_name]
        if reaction_energy_column is None:
            delta_g = [(des_a[0] + des_b[0]) * 627.5]
        else:
            delta_g = [float(row[reaction_energy_column])]
        all_xs.append(delta_g + des_a[1:] + des_b[1:])

    return np.asarray(all_xs, dtype=float)


def save_descriptor_maps(
    bn_map,
    cl_map,
    bn_path="Data/descriptor/BNdes_new2.pkl",
    cl_path="Data/descriptor/Cldes_new2.pkl",
):
    for path in [bn_path, cl_path]:
        folder = os.path.dirname(path)
        if folder:
            os.makedirs(folder, exist_ok=True)

    with open(bn_path, "wb") as f:
        pickle.dump(bn_map, f)
    with open(cl_path, "wb") as f:
        pickle.dump(cl_map, f)


def load_descriptor_maps(
    bn_path="Data/descriptor/BNdes_new2.pkl",
    cl_path="Data/descriptor/Cldes_new2.pkl",
):
    with open(bn_path, "rb") as f:
        bn_map = pickle.load(f)
    with open(cl_path, "rb") as f:
        cl_map = pickle.load(f)
    return bn_map, cl_map
