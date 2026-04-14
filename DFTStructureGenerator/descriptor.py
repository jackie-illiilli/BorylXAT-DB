from __future__ import annotations

import os
import pickle

import numpy as np
import pandas as pd
from ase.db import connect
from morfeus import BuriedVolume
from tqdm import tqdm

from . import Tool, mol_manipulation


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

MODELING_DESCRIPTOR_NAMES = [
    "deltaG", "B_N Combination Energy", "B_N react deltaG", "B_N B_charge",
    "B_N N_charge", "B_N B_spin", "B_N N_spin", "B_N B_Hirshfield_charge",
    "B_N N_Hirshfield_charge", "B_N B_N_dist", "B_N Alpha-HOMO-1",
    "B_N Alpha-HOMO", "B_N Alpha-LUMO", "B_N Alpha-LUMO+1",
    "B_N Beta-HOMO-1", "B_N Beta-HOMO", "B_N Beta-LUMO",
    "B_N Beta-LUMO+1", "B_N dipole", "B_N MaxOtherCharge",
    "B_N MinOtherCharge", "B_N_Cl B_charge", "B_N_Cl N_charge",
    "B_N_Cl Cl_charge", "B_N_Cl B_Hirshfield_charge",
    "B_N_Cl N_Hirshfield_charge", "B_N_Cl Cl_Hirshfield_charge",
    "B_N_Cl B_N_dist", "B_N_Cl N_Cl_dist", "B_N_Cl B_N_Cl_angel",
    "B_N_Cl HOMO-1", "B_N_Cl HOMO", "B_N_Cl LUMO", "B_N_Cl LUMO+1",
    "B_N_Cl dipole", "B_N_Cl MaxOtherCharge", "B_N_Cl MinOtherCharge",
    "B_N_Cl MaxChargeDist", "B_N_Cl MinOtherChargeDist",
    "B_N_Cl MaxOtherHCharge", "B_N_Cl MinOtherHCharge",
    "B_N_Cl MaxHChargeDist", "B_N_Cl MinHOtherChargeDist",
    "B_N_Cl MaxAromaticRingCharge", "B_N_Cl MinAromaticRingCharge",
    "B_N_Cl MaxAromaticChargeDist", "B_N_Cl MinAromaticChargeDist",
    "B_N_Cl Bv2radius", "B_N_Cl Bv4radius", "B_N_Cl Bv6radius",
    "C_Cl react deltaG", "C_Cl Cl_charge", "C_Cl C_charge",
    "C_Cl Cl_Hirshfield_charge", "C_Cl C_Hirshfield_charge",
    "C_Cl MaxOtherCharge", "C_Cl MinOtherCharge",
    "C_Cl MaxChargeDist", "C_Cl MinOtherChargeDist",
    "C_Cl MaxOtherHCharge", "C_Cl MinOtherHCharge",
    "C_Cl MaxHChargeDist", "C_Cl MinHOtherChargeDist",
    "C_Cl MaxAromaticRingCharge", "C_Cl MinAromaticRingCharge",
    "C_Cl MaxAromaticChargeDist", "C_Cl MinAromaticChargeDist",
    "C_Cl C_Cl_dist", "C_Cl HOMO-1", "C_Cl HOMO", "C_Cl LUMO",
    "C_Cl LUMO+1", "C_Cl dipole", "C_Cl Bv2radius", "C_Cl Bv4radius",
    "C_Cl Bv6radius", "C C_charge", "C C_spin", "C C_Hirshfield_charge",
    "C Alpha-HOMO-1", "C Alpha-HOMO", "C Alpha-LUMO",
    "C Alpha-LUMO+1", "C Beta-HOMO-1", "C Beta-HOMO", "C Beta-LUMO",
    "C Beta-LUMO+1", "C_dipole",
]

# Backward-compatible aliases for notebook variables.
duplicate_N_id = DUPLICATE_N_IDS
duplicate_Cl_id = DUPLICATE_CL_IDS
all_des_name = MODELING_DESCRIPTOR_NAMES


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
    db_path="boron_ccl2.db",
    bn_csv_path="data/csvs/reactants_B_N.csv",
    cl_csv_path="data/csvs/reactants_Cl.csv",
    duplicate_cl_ids=None,
    show_progress=True,
):
    if duplicate_cl_ids is None:
        duplicate_cl_ids = DUPLICATE_CL_IDS

    db = connect(db_path)
    bn_csv = pd.read_csv(bn_csv_path)
    cl_csv = pd.read_csv(cl_csv_path)

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
        cl_name, descriptor = _build_cl_descriptor(row, db, duplicate_cl_ids)
        cl_map[cl_name] = descriptor

    return bn_map, cl_map


def dataframe_to_descriptors(data_csv, bn_map, cl_map, duplicate_cl_ids=None, show_progress=True):
    if duplicate_cl_ids is None:
        duplicate_cl_ids = DUPLICATE_CL_IDS

    rows = data_csv.iterrows()
    if show_progress:
        rows = tqdm(rows, total=len(data_csv))

    all_xs = []
    for _, row in rows:
        b_index = int(row["B_Index"])
        n_index = int(row["N_Index"])
        cl_index = int(row["Cl_Index"])
        cl_atomid = int(row["Cl_Atomid"])

        b_n_name = f"B_{b_index:05}_Nu_{n_index:05}"
        cl_name = f"Cl_{cl_index:05}"
        if cl_index in duplicate_cl_ids:
            cl_name = f"Cl_{cl_index:05}_Claid_{cl_atomid:05}"

        des_a = bn_map[b_n_name]
        des_b = cl_map[cl_name]
        delta_g = [(des_a[0] + des_b[0]) * 627.5]
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
