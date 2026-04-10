from DFTStructureGenerator import B_N_Cl
import numpy as np
from tqdm import tqdm
import pandas as pd
from morfeus import BuriedVolume
from ase.db import connect

DB_PATH = "boron_ccl2.db"
db = connect(DB_PATH)

duplicate_N_id = [9, 43, 285, 310, 314, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 372, 375, 376]
duplicate_Cl_id = [624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636,
       637, 638, 639, 640, 642, 644, 645, 652, 653, 654, 655, 656, 657, 658, 659, 660, 661, 662, 663, 664,
       665, 666, 667, 668, 669, 670, 671, 672, 673, 674, 675, 676, 677,
       678, 679, 680, 681, 682, 683, 684, 685, 686, 687, 688, 689, 690,
       691, 692, 693, 694, 695, 696, 697, 698, 699, 700, 701, 702, 703,
       704, 705, 706, 707, 708, 709, 710, 711, 713, 714, 716, 717, 718, 719, 720, 721, 722]

all_des_name = ["deltaG", "B_N Combination Energy", "B_N react deltaG", "B_N B_charge", "B_N N_charge", "B_N B_spin", "B_N N_spin", 
 "B_N B_Hirshfield_charge", "B_N N_Hirshfield_charge", "B_N B_N_dist", "B_N Alpha-HOMO-1", "B_N Alpha-HOMO", 
 "B_N Alpha-LUMO", "B_N Alpha-LUMO+1", "B_N Beta-HOMO-1", "B_N Beta-HOMO", "B_N Beta-LUMO", "B_N Beta-LUMO+1",
 "B_N dipole", "B_N MaxOtherCharge", "B_N MinOtherCharge", "B_N_Cl B_charge", "B_N_Cl N_charge", "B_N_Cl Cl_charge", 
 "B_N_Cl B_Hirshfield_charge", "B_N_Cl N_Hirshfield_charge", "B_N_Cl Cl_Hirshfield_charge", "B_N_Cl B_N_dist", 
 "B_N_Cl N_Cl_dist", "B_N_Cl B_N_Cl_angel", "B_N_Cl HOMO-1", "B_N_Cl HOMO", "B_N_Cl LUMO", "B_N_Cl LUMO+1", "B_N_Cl dipole",
 "B_N_Cl MaxOtherCharge", "B_N_Cl MinOtherCharge", "B_N_Cl MaxChargeDist", "B_N_Cl MinOtherChargeDist", 
 "B_N_Cl MaxOtherHCharge", "B_N_Cl MinOtherHCharge", "B_N_Cl MaxHChargeDist", "B_N_Cl MinHOtherChargeDist", 
 "B_N_Cl MaxAromaticRingCharge", "B_N_Cl MinAromaticRingCharge", "B_N_Cl MaxAromaticChargeDist", "B_N_Cl MinAromaticChargeDist",
 "B_N_Cl Bv2radius", "B_N_Cl Bv4radius", "B_N_Cl Bv6radius", 
 "C_Cl react deltaG", "C_Cl Cl_charge", "C_Cl C_charge", "C_Cl Cl_Hirshfield_charge", "C_Cl C_Hirshfield_charge", 
 "C_Cl MaxOtherCharge", "C_Cl MinOtherCharge", "C_Cl MaxChargeDist", "C_Cl MinOtherChargeDist", 
 "C_Cl MaxOtherHCharge", "C_Cl MinOtherHCharge", "C_Cl MaxHChargeDist", "C_Cl MinHOtherChargeDist", 
 "C_Cl MaxAromaticRingCharge", "C_Cl MinAromaticRingCharge", "C_Cl MaxAromaticChargeDist", "C_Cl MinAromaticChargeDist",
 "C_Cl C_Cl_dist", 
 "C_Cl HOMO-1", "C_Cl HOMO", "C_Cl LUMO", "C_Cl LUMO+1", "C_Cl dipole", "C_Cl Bv2radius", "C_Cl Bv4radius", "C_Cl Bv6radius", "C C_charge", "C C_spin", "C C_Hirshfield_charge", 
 "C Alpha-HOMO-1", "C Alpha-HOMO", "C Alpha-LUMO", "C Alpha-LUMO+1", "C Beta-HOMO-1", "C Beta-HOMO", "C Beta-LUMO", "C Beta-LUMO+1", "C_dipole"]


selected_des = ['deltaG' 'B_N B_spin' 'B_N B_Hirshfield_charge' 'B_N B_N_dist'
 'B_N Alpha-HOMO' 'B_N_Cl B_Hirshfield_charge'
 'B_N_Cl Cl_Hirshfield_charge' 'B_N_Cl N_Cl_dist' 'B_N_Cl LUMO'
 'B_N_Cl Bv6radius' 'C_Cl Cl_Hirshfield_charge' 'C_Cl C_Hirshfield_charge'
 'C_Cl C_Cl_dist' 'C_Cl LUMO' 'C_Cl Bv6radius' 'C C_spin'
 'C C_Hirshfield_charge' 'C Alpha-HOMO']

B_N_csv = pd.read_csv("data/csvs/reactants_B_N.csv")
B_N_des_map = {}
for line_id, line in tqdm(B_N_csv.iterrows()):

    B_Index = int(line['B_Index'])
    B_Atomid = int(line['B_Atomid'])
    B_smiles = line['B_smiles']
    N_Index = int(line['N_Index'])
    N_Atomid = int(line['N_Atomid'])
    N_smiles = line['N_smiles']
    react_eng = line['deltaG_react']
    B_N_name = f"B_{B_Index:05}_Nu_{N_Index:05}"
    if B_N_name in B_N_des_map:
        continue

    # 从 smiles 获取 B 单体原子数，用于计算 N 在复合物中的原子偏移
    B_mol = B_N_Cl.mol_manipulation.smiles2mol(B_smiles)
    B_atomnum = B_mol.GetNumAtoms()

    # 从 db 读取 complex_r 和 complex_p 记录
    react_key = f"B_{B_Index:05}_LB_{N_Index:05}_r"
    prod_key = f"B_{B_Index:05}_LB_{N_Index:05}_p"
    react_row = db.get(key=react_key)
    prod_row = db.get(key=prod_key)

    react_atoms = react_row.toatoms()
    prod_atoms = prod_row.toatoms()
    react_positions = react_atoms.get_positions()
    prod_positions = prod_atoms.get_positions()
    react_symbols = react_atoms.get_chemical_symbols()
    prod_symbols = prod_atoms.get_chemical_symbols()

    react_spin_densities = react_row.data.get("spin_densities")
    react_hirshfeld = react_row.data.get("hirshfeld_charges")
    prod_hirshfeld = prod_row.data.get("hirshfeld_charges")

    # BN react 描述符
    descriptor = [react_eng]
    descriptor += [react_spin_densities[B_Atomid]]
    descriptor += [react_hirshfeld[B_Atomid]]
    descriptor += [B_N_Cl.Tool.get_atoms_distance(react_positions[B_Atomid], react_positions[N_Atomid + B_atomnum - 1])]
    descriptor += [react_row.homo_energy_kcal]

    # BNCl product 描述符
    Cl_Index = prod_symbols.index("Cl")
    descriptor += [prod_hirshfeld[B_Atomid], prod_hirshfeld[Cl_Index]]
    descriptor += [B_N_Cl.Tool.get_atoms_distance(prod_positions[B_Atomid], prod_positions[Cl_Index])]
    descriptor += [prod_row.lumo_energy_kcal]

    bv = BuriedVolume(prod_symbols, prod_positions, B_Atomid + 1, include_hs=1, radius=6, z_axis_atoms=[Cl_Index + 1], excluded_atoms=[Cl_Index + 1])
    bv.octant_analysis()
    descriptor += [np.sum(list(bv.octants['percent_buried_volume'].values()))]
    
    B_N_des_map[B_N_name] = descriptor

Cl_csv = pd.read_csv("data/csvs/reactants_Cl.csv")
Cl_des_map = {}
for line_id, line in tqdm(Cl_csv.iterrows()):
    Cl_Index = line['Index']
    Cl_Atomid = line['Atomid']
    Cl_smiles = line['Smiles']
    react_eng = line['deltaG_react']

    # 从 smiles 获取 C 原子索引（Cl 的邻居）
    Cl_mol = B_N_Cl.mol_manipulation.smiles2mol(Cl_smiles)
    C_atom_idx = Cl_mol.GetAtomWithIdx(Cl_Atomid).GetNeighbors()[0].GetIdx()

    # 从 db 读取 Cl_r 和 c_radical (Cl_p) 记录
    react_key = f"Cl_{Cl_Index:05}_r"
    prod_key = f"Cl_{Cl_Index:05}_p"
    react_row = db.get(key=react_key)
    prod_row = db.get(key=prod_key)

    react_atoms = react_row.toatoms()
    prod_atoms = prod_row.toatoms()
    react_positions = react_atoms.get_positions()
    prod_positions = prod_atoms.get_positions()
    react_symbols = react_atoms.get_chemical_symbols()

    react_hirshfeld = react_row.data.get("hirshfeld_charges")
    prod_spin_densities = prod_row.data.get("spin_densities")
    prod_hirshfeld = prod_row.data.get("hirshfeld_charges")

    # C-Cl react 描述符
    descriptor = [react_eng]
    descriptor += [react_hirshfeld[Cl_Atomid], react_hirshfeld[C_atom_idx]]

    descriptor += [B_N_Cl.Tool.get_atoms_distance(react_positions[Cl_Atomid], react_positions[C_atom_idx])]
    descriptor += [react_row.lumo_energy_kcal]
    bv = BuriedVolume(react_symbols, react_positions, C_atom_idx + 1, include_hs=1, radius=6, z_axis_atoms=[Cl_Atomid + 1], excluded_atoms=[Cl_Atomid + 1])
    bv.octant_analysis()
    descriptor += [np.sum(list(bv.octants['percent_buried_volume'].values()))]

    # Cl product (C radical) 描述符
    if C_atom_idx > Cl_Atomid: C_atom_idx -= 1
    descriptor += [prod_spin_densities[C_atom_idx]]
    descriptor += [prod_hirshfeld[C_atom_idx]]
    descriptor += [prod_row.homo_energy_kcal]
    Cl_name = f"Cl_{Cl_Index:05}"
    if Cl_Index in duplicate_Cl_id:
        Cl_name = f"Cl_{Cl_Index:05}_Claid_{Cl_Atomid:05}"
    Cl_des_map[Cl_name] = descriptor