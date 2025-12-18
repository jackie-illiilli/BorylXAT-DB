import glob, os, shutil, itertools, copy
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D
import numpy as np
from tqdm import tqdm
import pandas as pd
from matplotlib import pyplot as plt
from . import logfile_process, FormatConverter, xtb_process, mol_manipulation, Tool



def get_react_atomid(smiles, type_='B'):
    """根据常见的几类B、配体、氯确定反应位点
    对于氯和配体，会获取不同化学环境的多个反应位点

    Args:
        smiles (str): SMILES
        type_ (str, optional): in ["B", "Cl", "NHC", "NHCS", "Ar", "N", "P"]. Defaults to 'B'.

    Returns:
        _type_: _description_
    """    
    assert type_ in ["B", "Cl", "NHC", "NHCS", "Ar", "N", "P"]
    if type_ in ["B", "P"]:
        b_mol = Chem.MolFromSmiles(smiles)
        b_atom_idx = [each for each in b_mol.GetAtoms() if each.GetSymbol() == type_][0].GetIdx()
        return [b_atom_idx]
    elif type_ == "NHC":
        nu_mol = Chem.MolFromSmiles(smiles)
        nu_mol = Chem.AddHs(nu_mol)
        nu_atom_idx = [each for each in nu_mol.GetAtoms() if each.GetSymbol() == 'C' and each.GetExplicitValence() == 2][0].GetIdx()
        return [nu_atom_idx]
    elif type_ == "NHCS":
        nu_mol = Chem.MolFromSmiles(smiles)
        nu_atom_idx = [each for each in nu_mol.GetAtoms() if each.GetSymbol() == 'C' and each.GetFormalCharge() == -1][0].GetIdx()
        return [nu_atom_idx]
    else:
        nu_mol = Chem.MolFromSmiles(smiles)
        if type_ == "Ar":
            nu_atom_idx = [each.GetIdx() for each in nu_mol.GetAtoms() if each.GetSymbol() == 'N' and len(each.GetNeighbors()) > 1 and each.GetIsAromatic()]
            if len(nu_atom_idx) == 0:
                nu_atom_idx = [each.GetIdx() for each in nu_mol.GetAtoms() if each.GetSymbol() == 'N' and len(each.GetNeighbors()) > 1]
        elif type_ == "Cl":
            nu_atom_idx = [each.GetIdx() for each in nu_mol.GetAtoms() if each.GetSymbol() == 'Cl']
        else:
            nu_atom_idx = [each.GetIdx() for each in nu_mol.GetAtoms() if each.GetSymbol() == 'N' and len(each.GetNeighbors()) > 1]
        # 挑选出不同化学环境的配位原子，删除具有对称性的配位原子之一
        remain_smiles = []
        remain_atom_ids = []
        for each_idx in nu_atom_idx:
            mol = Chem.MolFromSmiles(smiles)
            mol.GetAtomWithIdx(each_idx).SetAtomicNum(33)
            new_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
            if new_smiles not in remain_smiles:
                remain_smiles.append(new_smiles)
                remain_atom_ids.append(each_idx)
        return remain_atom_ids

def generate_combinations(reactant_file, Bresult_file, Nresult_file, Clresult_file):
    """根据给定的excel，获取反应物、反应位点、编号便于统计

    Args:
        reactant_file (str): 原始Excel的路径
        Bresult_file (str): 要存储B自由基信息的文件路径
        Nresult_file (_type_): 要存储配体信息的文件路径
        Clresult_file (_type_): 要存储氯化物信息的文件路径

    Returns:
        list: 含有多种反应位点的配体的编号
    """    
    duplicate_N = []
    duplicate_Cl = []
    react_csv = pd.read_excel(reactant_file)
    N_smiles = react_csv["Anime"].dropna().to_numpy()
    N_index = react_csv['Anime_Index'].dropna().to_numpy()
    P_smiles = react_csv["Phosphine"].dropna().to_numpy()
    P_index = react_csv['Phosphine_Index'].dropna().to_numpy()
    NHC_smiles = react_csv["NHC"].dropna().to_numpy()
    NHC_index = react_csv['NHC_Index'].dropna().to_numpy()
    NHCS_smiles = react_csv["NHCS"].dropna().to_numpy()
    NHCS_index = react_csv['NHCS_Index'].dropna().to_numpy()
    Ar_smiles = react_csv["Ar"].dropna().to_numpy()
    Ar_index = react_csv['Ar_Index'].dropna().to_numpy()
    B_smiles = react_csv["BR2"].dropna().to_numpy()
    B_index = react_csv['BR2_Index'].dropna().to_numpy()
    Cl_smiles = react_csv["Cl"].dropna().to_numpy()
    Cl_index = react_csv['Cl_Index'].dropna().to_numpy()

    if os.path.exists(Clresult_file):
        result_csv = pd.read_csv(Clresult_file).to_dict()
    else:
        result_csv = {"Smiles":{}, "Index":{}, "Atomid":{}}
    result_csv_id = len(result_csv["Index"])
    for smiles, id in zip(Cl_smiles, Cl_index):
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        if smiles in result_csv["Smiles"]:
            continue
        cl_atom_ids = get_react_atomid(smiles, type_='Cl')
        if len(cl_atom_ids) > 1:
            print(smiles, id, cl_atom_ids)
            duplicate_Cl.append(int(id))
        for cl_atom_id in cl_atom_ids:
            result_csv["Smiles"][result_csv_id] = smiles
            result_csv['Index'][result_csv_id] = int(id)
            result_csv['Atomid'][result_csv_id] = cl_atom_id
            result_csv_id += 1
    result_csv = pd.DataFrame(result_csv)
    result_csv.to_csv(Clresult_file, index=False)
    
    if os.path.exists(Bresult_file):
        result_csv = pd.read_csv(Bresult_file).to_dict()
    else:
        result_csv = {"Smiles":{}, "Index":{}, "Atomid":{}}
    result_csv_id = len(result_csv["Index"])
    for smiles, id in zip(B_smiles, B_index):
        smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        result_csv["Smiles"][result_csv_id] = smiles
        result_csv['Index'][result_csv_id] = int(id)
        result_csv['Atomid'][result_csv_id] = get_react_atomid(smiles, type_='B')[0]
        result_csv_id += 1
    result_csv = pd.DataFrame(result_csv)
    result_csv.to_csv(Bresult_file, index=False)

    if os.path.exists(Nresult_file):
        result_csv = pd.read_csv(Nresult_file).to_dict()
    else:
        result_csv = {"Smiles":{}, "Index":{}, "Atomid":{}}
    result_csv_id = len(result_csv["Index"])
    Nus = [[N_smiles, N_index], [P_smiles, P_index], [NHC_smiles, NHC_index], [NHCS_smiles, NHCS_index], [Ar_smiles, Ar_index]]
    for nu_id, nu in enumerate(Nus):
        for smiles, id in zip(nu[0], nu[1]):
            smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
            n_atom_ids = get_react_atomid(smiles, type_=["N", "P", "NHC", "NHCS", "Ar"][nu_id])
            if len(n_atom_ids) !=1:
                print(smiles, id, n_atom_ids)
                duplicate_N.append(id)
            for n_atom_id in n_atom_ids:
                result_csv["Smiles"][result_csv_id] = smiles
                result_csv['Index'][result_csv_id] = int(id)
                result_csv['Atomid'][result_csv_id] = n_atom_id
                result_csv_id += 1
    result_csv = pd.DataFrame(result_csv)
    result_csv.to_csv(Nresult_file, index=False)
    return duplicate_N, duplicate_Cl

def B_N_Single_Xtb(root_file, Bresult_file, Nresult_file, mol_xtb_name = 'Mol_xtb', mol_name = "Mols"):
    """产生B自由基、配体单体的Xtb优化文件

    Args:
        Bresult_file (str): 存储B自由基信息的文件路径
        Nresult_file (str): 存储配体信息的文件路径
    """    
    Bresult_file=pd.read_csv(Bresult_file)
    Nresult_file=pd.read_csv(Nresult_file)
    mol_xtb_file = os.path.join(root_file, mol_xtb_name)
    old_mol_file = os.path.join(root_file, 'Mols')
    mol_file = os.path.join(root_file, mol_name)
    result_files = [Bresult_file, Nresult_file]
    if not os.path.isdir(mol_file):
        os.mkdir(mol_file)
    if not os.path.isdir(mol_xtb_file):
        os.mkdir(mol_xtb_file)
    all_Bmols = []
    all_Bnames = []
    all_Nmols = []
    all_Nnames = []
    for name, result_file in zip(["B", "N"], result_files):
        for ix, row in result_file.iterrows():
            mol_idx = row['Index']
            mol_name = f"{name}_{mol_idx:05}_r"
            if os.path.exists(os.path.join(old_mol_file, f"{mol_name}.mol")):
                continue
            mol = mol_manipulation.smiles2mol(row['Smiles'], conf_num=1)
            mol_atom_idx = row["Atomid"]
            AllChem.UFFOptimizeMolecule(mol)
            if name == "B":
                b_rwmol = Chem.RWMol(mol)
                h_atom_idx = [each for each in mol.GetAtomWithIdx(mol_atom_idx).GetNeighbors() if each.GetSymbol() == 'H'][0].GetIdx()
                b_rwmol.GetAtomWithIdx(mol_atom_idx).SetNumRadicalElectrons(1)
                b_rwmol.RemoveAtom(h_atom_idx)
                mol = b_rwmol.GetMol()
            Chem.MolToMolFile(mol, os.path.join(mol_file, f"{mol_name}.mol"))
            if name == "B":
                all_Bmols.append(mol)
                all_Bnames.append(mol_name)
            else:
                all_Nmols.append(mol)                                                                             (mol)
                all_Nnames.append(mol_name)
    mol_xtb_file_ = os.path.join(mol_xtb_file, 'B_single')
    xtb_process.xtb_main(all_Bnames, all_Bmols, dir_path=mol_xtb_file_, core=60)
    xtb_process.shift_to_sugan(mol_xtb_file_, 1, 0, 1)
    mol_xtb_file_ = os.path.join(mol_xtb_file, 'N_single')
    xtb_process.xtb_main(all_Nnames, all_Nmols, dir_path=mol_xtb_file_, core=60)
    xtb_process.shift_to_sugan(mol_xtb_file_, 1, 0, 0)

def B_N_Cl_reactant_product_Xtb(root_file, Bresult_file, Nresult_file, Clresult_file, duplicate_N_id, duplicate_Cl_id, mol_xtb_name = 'Mol_xtb', mol_name = "Mols"):
    """产生氯反应物和产物自由基的Xtb优化文件。
    搭建B-N复合物，和Cl-B-N产物的Xtb优化文件

    Args:
        Bresult_file (str): 要存储B自由基信息的文件路径
        Nresult_file (_type_): 要存储配体信息的文件路径
        Clresult_file (_type_): 要存储氯化物信息的文件路径
        duplicate_N_id (list): 含有多种反应位点的配体编号
    """    
    mol_xtb_file = os.path.join(root_file, mol_xtb_name)
    old_mol_file = os.path.join(root_file, 'Mols')
    mol_file = os.path.join(root_file, mol_name)
    if not os.path.isdir(mol_file):
        os.mkdir(mol_file)
    Bresult_file=pd.read_csv(Bresult_file)
    Nresult_file=pd.read_csv(Nresult_file)
    Clresult_file=pd.read_csv(Clresult_file)
    # Generate Cl Reactants / Products:
    Cl_r_mols = []
    Cl_r_names = []
    B_N_r_mols = []
    B_N_r_names = []
    Cl_p_mols = []
    Cl_p_names = []
    Cl_p_d_mols = []
    Cl_p_d_names = []
    B_N_p_mols = []
    B_N_p_names = []
    B_N_p_d_mols = []
    B_N_p_d_names = []
    B_N_r_d_mols = []
    B_N_r_d_names = []
    for ix, row in Clresult_file.iterrows():
        cl_mol_idx = row['Index']
        mol_name = f"Cl_{cl_mol_idx:05}_r"
        cl_mol = mol_manipulation.smiles2mol(row['Smiles'], conf_num=1)
        if os.path.exists(os.path.join(old_mol_file, f"{mol_name}.mol")):
            continue
        cl_atom_idx = row["Atomid"]
        AllChem.UFFOptimizeMolecule(cl_mol)
        Chem.MolToMolFile(cl_mol, os.path.join(mol_file, f"{mol_name}.mol"))
        Cl_r_names.append(mol_name)
        Cl_r_mols.append(cl_mol)
        # Cl Products
        C_atom_idx = cl_mol.GetAtomWithIdx(cl_atom_idx).GetNeighbors()[0].GetIdx()
        cl_rwmol = Chem.RWMol(cl_mol)
        cl_rwmol.GetAtomWithIdx(C_atom_idx).SetAtomMapNum(10)
        cl_rwmol.RemoveAtom(cl_atom_idx)
        new_C_atom = [each for each in cl_rwmol.GetAtoms() if each.GetAtomMapNum() == 10][0]
        new_C_atom.SetNumRadicalElectrons(1)
        new_C_atom.SetAtomMapNum(0)
        cl_product_mol = cl_rwmol.GetMol()
        mol_name = f"Cl_{cl_mol_idx:05}_p"
        if cl_mol_idx in duplicate_Cl_id:
            mol_name = f"Cl_{cl_mol_idx:05}_Claid_{cl_atom_idx:05}_p"
            Cl_p_d_names.append(mol_name)
            Cl_p_d_mols.append(cl_product_mol)
        else:
            Cl_p_names.append(mol_name)
            Cl_p_mols.append(cl_product_mol)
        Chem.MolToMolFile(cl_product_mol, os.path.join(mol_file, f"{mol_name}.mol"))

    # Generate B-N Reactants

    for Bix, Brow in Bresult_file.iterrows():
        for Nix, Nrow in Nresult_file.iterrows():
            b_mol_idx = Brow['Index']
            nu_mol_idx = Nrow['Index']
            b_atom_idx = Brow["Atomid"]
            nu_atom_idx = Nrow['Atomid']
            if nu_mol_idx not in duplicate_N_id:
                mol_name_p = f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_p"
            else:
                mol_name_p = f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_Naid_{nu_atom_idx:05}_p"
            if os.path.exists(os.path.join(old_mol_file, f"{mol_name_p}.mol")):
                continue            
            b_mol = mol_manipulation.smiles2mol(Brow['Smiles'], conf_num=1)
            nu_mol = mol_manipulation.smiles2mol(Nrow['Smiles'], conf_num=1)
            Chem.Kekulize(b_mol, 1)
            Chem.Kekulize(nu_mol, 1)
            try:
                # if nu_mol_idx not in duplicate_N_id:
                #     continue
                # if nu_mol_idx not in [162,165,190,224]:
                #     continue
                # Connect B and Nu
                b_atom_num = b_mol.GetNumAtoms()
                
                h_atom_idx = [each for each in b_mol.GetAtomWithIdx(b_atom_idx).GetNeighbors() if each.GetSymbol() == 'H'][0].GetIdx()
                
                nu_mol = mol_manipulation.move_mol(nu_mol, array=np.array([5, 5, 5]))
                b_nu_mol = Chem.CombineMols(b_mol, nu_mol)
                b_nu_rwmol = Chem.RWMol(b_nu_mol)
                b_nu_rwmol.AddBond(b_atom_idx, nu_atom_idx + b_atom_num, Chem.rdchem.BondType.SINGLE)
                b_nu_rwmol.GetAtomWithIdx(b_atom_idx).SetFormalCharge(-1)
                b_nu_rwmol.GetAtomWithIdx(nu_atom_idx + b_atom_num).SetFormalCharge(1)
                b_nu_product_rwmol = copy.deepcopy(b_nu_rwmol)
                b_nu_rwmol.GetAtomWithIdx(b_atom_idx).SetNumRadicalElectrons(1)
                b_nu_rwmol.GetAtomWithIdx(nu_atom_idx + b_atom_num).SetNumRadicalElectrons(0)
                b_nu_rwmol.RemoveAtom(h_atom_idx)
                b_nu_mol = b_nu_rwmol.GetMol()
                AllChem.SanitizeMol(b_nu_mol)
                AllChem.Kekulize(b_nu_mol)
                AllChem.EmbedMolecule(b_nu_mol)
                AllChem.MMFFOptimizeMolecule(b_nu_mol)
                nb_atom_idx = [each for each in b_nu_mol.GetAtoms() if each.GetSymbol() == 'B'][0].GetIdx()

                # Product
                b_nu_product_rwmol.GetAtomWithIdx(h_atom_idx).SetAtomicNum(17)
                b_nu_product_rwmol.GetAtomWithIdx(b_atom_idx).SetNumRadicalElectrons(0)
                b_nu_product_rwmol.GetAtomWithIdx(nu_atom_idx + b_atom_num).SetNumRadicalElectrons(0)
                b_nu_rwmol.GetAtomWithIdx(b_atom_idx).SetFormalCharge(-1)
                b_nu_rwmol.GetAtomWithIdx(nu_atom_idx + b_atom_num).SetFormalCharge(1)
                b_nu_product_mol = b_nu_product_rwmol.GetMol()
                AllChem.SanitizeMol(b_nu_product_mol)
                AllChem.Kekulize(b_nu_product_mol)
                AllChem.EmbedMolecule(b_nu_product_mol)
                AllChem.MMFFOptimizeMolecule(b_nu_product_mol)

                if nu_mol_idx not in duplicate_N_id:
                    mol_name_p = f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_p"
                    Chem.MolToMolFile(b_nu_product_mol, os.path.join(mol_file, f"{mol_name_p}.mol"), kekulize=False)
                    mol_name_r = f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_r"
                    Chem.MolToMolFile(b_nu_mol, os.path.join(mol_file, f"{mol_name_r}.mol"), kekulize=False)
                    B_N_r_mols.append(b_nu_mol)
                    B_N_r_names.append(mol_name_r)
                    B_N_p_mols.append(b_nu_product_mol)
                    B_N_p_names.append(mol_name_p)
                else:
                    mol_name_p = f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_Naid_{nu_atom_idx:05}_p"
                    Chem.MolToMolFile(b_nu_product_mol, os.path.join(mol_file, f"{mol_name_p}.mol"), kekulize=False)
                    mol_name_r = f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_Naid_{nu_atom_idx:05}_r"
                    Chem.MolToMolFile(b_nu_mol, os.path.join(mol_file, f"{mol_name_r}.mol"), kekulize=False)
                    B_N_r_d_mols.append(b_nu_mol)
                    B_N_r_d_names.append(mol_name_r)
                    B_N_p_d_mols.append(b_nu_product_mol)
                    B_N_p_d_names.append(mol_name_p)
            except:
                with open("test.txt", 'at') as f:
                    f.write(str(b_mol_idx)+ " "+ str(nu_mol_idx) + " " + str(nu_atom_idx) + "Fail to Generate\n")
    all_uhf = [0, 1, 1, 1, 0, 0, 1]
    all_strs = ['Cl_r', 'B_N_r', 'Cl_p', "Cl_p_d", 'B_N_p', 'B_N_p_d', 'B_N_r_d']
    all_names = [Cl_r_names, B_N_r_names, Cl_p_names, Cl_p_d_names, B_N_p_names, B_N_p_d_names, B_N_r_d_names]
    all_mols = [Cl_r_mols, B_N_r_mols, Cl_p_mols, Cl_p_d_mols, B_N_p_mols, B_N_p_d_mols, B_N_r_d_mols]
    for str_, uhf, names, mols in zip(all_strs, all_uhf, all_names, all_mols):
        mol_xtb_file_ = os.path.join(mol_xtb_file, str_)
        if not os.path.isdir(mol_xtb_file_):
            os.mkdir(mol_xtb_file_)
        xtb_process.xtb_main(names, mols, dir_path=mol_xtb_file_, core=60, uhf=1)
        # xtb_process.shift_to_sugan(mol_xtb_file_, 1, 0, uhf)  

def smiles_DFT_calc(root_dir='first_xtb', 
                    mol_dir='mol', 
                    dft_dir='mol_dft', 
                    method="opt freq b3lyp/6-31g* em=gd3bj",
                    conf_limit=3,
                    rmsd_limit=1.5,
                    SpinMultiplicity = None
                    ):
    """通过Xtb结果，优化得到Gaussian优化输入文件

    Args:
        root_dir (str, optional): Xtb的根目录. Defaults to 'first_xtb'.
        mol_dir (str, optional): Mol分子的目录. Defaults to 'mol'.
        dft_dir (str, optional): 要存储Gaussian输入文件的目录. Defaults to 'mol_dft'.
        method (str, optional): Gaussian方法. Defaults to "opt freq b3lyp/6-31g* em=gd3bj".
        conf_limit (int, optional): Xtb读取结构的构象数量限制. Defaults to 3.
        rmsd_limit (float, optional): Xtb读取结构的RMSD限制. Defaults to 1.5.
        SpinMultiplicity (int, optional): 设定的自旋多重度. Defaults to None.
    """                    
    all_files = glob.glob(root_dir + "/*/*/*")
    for xtb_file in all_files:
        if ("crest.out" in xtb_file) or ("best" in xtb_file) or ("crest_conf" in xtb_file):
            pass
        else:
            if os.path.isdir(xtb_file):
                shutil.rmtree(xtb_file)
            else:
                os.remove((xtb_file))
    xtb_dirs = glob.glob(root_dir + "/*/*")
    for i, xtb_dir in enumerate(xtb_dirs):
        mol_name = os.path.split(xtb_dir)[-1][:-2]
        mol_file = mol_dir + f"/{mol_name}.mol" 
        mol = Chem.MolFromMolFile(mol_file, removeHs=False, sanitize=False)
        title = "Singlemol"
        xtb_process.after_xtb(mol,xtb_dir=xtb_dir, save_dir=dft_dir, xtb_title=title, method=method, conf_limit=conf_limit, rmsd_limit=rmsd_limit, SpinMultiplicity=SpinMultiplicity)


def SPE_DFT_calc(target_dir, opt_name="Reactants", eng_name='Reactants_eng', save_chk=None, method="b3lyp/6-311+g(d,p) em=gd3bj"):
    opt_file_dir = os.path.join(target_dir, opt_name)
    eng_dir = os.path.join(target_dir, eng_name)
    log_files = glob.glob(opt_file_dir + "/" + "*.log")
    for log_file in log_files:
        new_log_name = eng_dir + "/" + os.path.split(log_file)[-1].split('.')[0] + ".gjf" 
        opt_log = logfile_process.Logfile(log_file)
        assert len(opt_log.running_positions) != 0
        title, charge, symbol_list, position,= opt_log.title, opt_log.charge, opt_log.symbol_list, opt_log.running_positions[-1]
        title = " ".join(str(each) for each in title)
        if save_chk:
            savechk = os.path.split(new_log_name.strip(".gjf"))[-1]
        else:
            savechk = None
        FormatConverter.block_to_gjf(symbol_list, position, new_log_name, charge, opt_log.multiplicity, title,
                    method=method, savechk=savechk)



def SPE_DFT_calc_wfn(target_dir, opt_name="Reactants", eng_name='Reactants_eng', mol_name="Mols", save_chk=None, method="b3lyp/6-311+g(d,p) em=gd3bj"):
    opt_file_dir = os.path.join(target_dir, opt_name)
    eng_dir = os.path.join(target_dir, eng_name)
    mol_files = glob.glob(os.path.join(target_dir, mol_name, "*.mol"))
    for mol_file in tqdm(mol_files):
        log_files = glob.glob(opt_file_dir + "/" + os.path.split(mol_file)[-1].split(".")[0] + "*.log")
        if len(log_files) == 0:
            continue 
        for log_file in log_files:
            new_log_name = eng_dir + "/" + os.path.split(log_file)[-1].split('.')[0] + ".gjf" 
            opt_log = logfile_process.Logfile(log_file, mol_file_dir=mol_file)
            assert len(opt_log.running_positions) != 0
            wfn_name = os.path.split(log_file)[-1].split('.')[0] + ".wfn"
            opt_log = logfile_process.Logfile(log_file, mol_file_dir=mol_file)
            assert len(opt_log.running_positions) != 0
            title, charge, symbol_list, position,= opt_log.title, opt_log.charge, opt_log.symbol_list, opt_log.running_positions[-1]
            title = " ".join(str(each) for each in title)
            FormatConverter.block_to_gjf(symbol_list, position, new_log_name, charge, opt_log.multiplicity, title,
                        method=method, final_line=wfn_name)


def collection_dft_single(Bresult_path, Nresult_path, Clresult_path, mol_dir, dft_dir, spe_dir, duplicate_Cl_id):
    # Bresult_path ='Data/All_Data/reactants_B.csv'
    # Nresult_path = 'Data/All_Data/reactants_N.csv'
    # Clresult_path = 'Data/All_Data/reactants_Cl.csv'
    result_paths = [Bresult_path, Nresult_path, Clresult_path, Clresult_path]
    types = ['r', 'r', 'r', 'p']
    path_strs = ["B_single", 'N_single', 'Cl_r', 'Cl_p']
    mol_strs = ["B", "N", "Cl", "Cl"]
    for result_path, type_, path_str, mol_str in zip(result_paths, types, path_strs, mol_strs):
        error_reason, E_energy, G_energy, conf_idxs = [],[],[], []
        result_file = pd.read_csv(result_path)
        for line_id, line in tqdm(result_file.iterrows()):
            dft_dir_ = os.path.join(dft_dir, path_str)
            spe_dir_ = os.path.join(spe_dir, path_str)
            smiles = line['Smiles']
            Index = line['Index']
            atom_id = line['Atomid']
            if mol_str == 'Cl' and  Index in duplicate_Cl_id and type_ == 'p':
                dft_dir_ = os.path.join(dft_dir, "Cl_p_d")
                spe_dir_ = os.path.join(spe_dir, "Cl_p_d")
                mol_file = glob.glob(os.path.join(mol_dir, f'{mol_str}_{Index:05}_Claid_{atom_id:05}_{type_}.mol'))
            else:
                mol_file = glob.glob(os.path.join(mol_dir, f'{mol_str}_{Index:05}_{type_}.mol'))
            if len(mol_file) == 0:
                print(smiles, Index, "Is Error")
                # return None
            mol_file = mol_file[0]
            if mol_str == 'Cl' and  Index in duplicate_Cl_id and type_ == 'p':
                opt_files = glob.glob(os.path.join(dft_dir_, f'{mol_str}_{Index:05}_Claid_{atom_id:05}_{type_}*.log'))
            else:
                opt_files = glob.glob(os.path.join(dft_dir_, f'{mol_str}_{Index:05}_{type_}*.log'))
            temp_idx, temp_E, temp_G = [], [], []
            for opt_file in opt_files:
                
                opt = mol_manipulation.logfile_process.Logfile(opt_file)
                spe_files = glob.glob(os.path.join(spe_dir_, os.path.split(opt_file)[-1]))
                if len(spe_files) == 0:
                    print(smiles, Index, "OPT Error")
                    continue
                conf_id = int(opt_file.split('.')[0].split("_")[-1])
                spe_file = spe_files[0]
                spe = mol_manipulation.logfile_process.Logfile(spe_file)
                electric_energy = spe.all_engs[0]
                G_cor = opt.all_engs[-1]
                temp_idx.append(conf_id)
                temp_E.append(electric_energy)
                temp_G.append(G_cor + electric_energy)
            if len(temp_G) == 0:
                error_reason.append("DFT Error")
                E_energy.append(np.nan)
                G_energy.append(np.nan)
                conf_idxs.append(np.nan)
            else:
                min_index = np.argmin(temp_G)
                E_energy.append(temp_E[min_index])
                G_energy.append(temp_G[min_index])
                conf_idxs.append(temp_idx[min_index])
                error_reason.append(np.nan)
        result_file[f"E_energy_{type_}"] = E_energy
        result_file[f"G_energy_{type_}"] = G_energy
        result_file[f"conf_idxs_{type_}"] = conf_idxs
        result_file[f"error_reason_{type_}"] = error_reason
        result_file.to_csv(result_path, index=False)

def collection_dft_ts(ts_csv_path, dft_dir, spe_dir, duplicate_Cl_id = [624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 642, 644, 645], bond_dist=True):
    # Bresult_path ='Data/All_Data/reactants_B.csv'
    # Nresult_path = 'Data/All_Data/reactants_N.csv'
    # Clresult_path = 'Data/All_Data/reactants_Cl.csv'
    error_reason, G_energy, conf_idxs, bond_1, bond_2 = [], [], [], [], []
    result_file = pd.read_csv(ts_csv_path)
    G_energy_react = result_file["G_energy"]
    for line_id, line in tqdm(result_file.iterrows()):
        B_idx = line['B_Index']
        N_idx = line['N_Index']
        Cl_idx = line['Cl_Index']
        cl_atom_idx = line['Cl_Atomid']
        react_name = f"B_{B_idx:05}_Nu_{N_idx:05}_Cl_{Cl_idx:05}"
        if Cl_idx in duplicate_Cl_id:
            react_name = f"B_{B_idx:05}_Nu_{N_idx:05}_Cl_{Cl_idx:05}_Claid_{cl_atom_idx:05}"
        opt_files = glob.glob(os.path.join(dft_dir, f'{react_name}*.log'))
        temp_idx, temp_G = [], []
        for opt_file in opt_files:
            opt = mol_manipulation.logfile_process.Logfile(opt_file)
            title, position = opt.title, opt.running_positions[-1]
            spe_files = glob.glob(os.path.join(spe_dir, os.path.split(opt_file)[-1]))
            if len(spe_files) == 0:
                print(react_name, "SPE Error")
                continue
            file_split = os.path.split(opt_file)[-1].split('.')[0].split("_")
            if len(file_split) == 6:
                conf_id = -1
            else:
                conf_id = int(file_split[-1])
            spe_file = spe_files[0]
            spe = mol_manipulation.logfile_process.Logfile(spe_file)
            electric_energy = spe.all_engs[0]
            G_cor = opt.all_engs[-1]
            temp_idx.append(conf_id)
            temp_G.append(G_cor + electric_energy)
        if len(temp_G) == 0:
            error_reason.append("DFT Error")
            G_energy.append(np.nan)
            conf_idxs.append(np.nan)
            bond_1.append(np.nan)
            bond_2.append(np.nan)
        else:
            min_index = np.argmin(temp_G)
            G_energy.append(temp_G[min_index])
            conf_idxs.append(temp_idx[min_index])
            error_reason.append(np.nan)
            if bond_dist:
                bond_1.append(Tool.get_atoms_distance(position[title[0] - 1], position[title[1] - 1]))
                bond_2.append(Tool.get_atoms_distance(position[title[1] - 1], position[title[2] - 1]))
            else:
                bond_1.append(np.nan)
                bond_2.append(np.nan)

    result_file[f"G_energy_ts"] = G_energy
    result_file[f"conf_idxs_ts"] = conf_idxs
    result_file[f"error_reason_ts"] = error_reason
    result_file[f"deltaGa"] = [627.5 * (a - b) if not np.isnan(a) else np.nan for a, b in zip(G_energy, G_energy_react) ]
    result_file[f"bond_1"] = bond_1
    result_file[f"bond_2"] = bond_2
    result_file.to_csv(ts_csv_path, index=False)

# 复合物整理
def collection_dft_couple(
    mol_dir, dft_dir, spe_dir,duplicate_N_id,
    Bresult_path = 'Data/All_Data/reactants_B.csv',
    Nresult_path = 'Data/All_Data/reactants_N.csv',
    B_N_result_path = 'Data/All_Data/reactants_B_N.csv',
    type_ = 'r',
    ):

    # type_ = "r"
    # Bresult_path = 'Data/All_Data/reactants_B.csv'
    # Nresult_path = 'Data/All_Data/reactants_N.csv'
    B_N_result_path = 'Data/All_Data/reactants_B_N.csv'
    B_N_result_idx = 0
    Bresult_file = pd.read_csv(Bresult_path)
    Nresult_file = pd.read_csv(Nresult_path)
    if type_ == 'r':
        B_N_result_file = {"B_smiles":{}, "B_Index":{}, "B_Atomid":{}, 
                        "N_smiles":{}, "N_Index":{}, "N_Atomid":{},
                        "conf_idxs_r":{}, "conf_idxs_p":{}, "Error_reason":{},
                        "E_energy_r":{}, "G_energy_r":{},"E_energy_p":{}, "G_energy_p":{},  
                        "deltaE_comb(kcal)":{}, "deltaG_comb(kcal)":{}, "deltaE_react":{}, "deltaG_react":{}}
    else:
        B_N_result_file = pd.read_csv(B_N_result_path)
    for Bix, Brow in tqdm(Bresult_file.iterrows()):
        for Nix, Nrow in Nresult_file.iterrows():
            b_smiles = Brow['Smiles']
            n_smiles = Nrow['Smiles']
            b_mol_idx = Brow['Index']
            nu_mol_idx = Nrow['Index']
            b_atom_idx = Brow["Atomid"]
            nu_atom_idx = Nrow['Atomid']
            if type_ == 'r':
                B_N_result_file['B_smiles'][B_N_result_idx] = b_smiles
                B_N_result_file['B_Index'][B_N_result_idx] = b_mol_idx
                B_N_result_file['B_Atomid'][B_N_result_idx] = b_atom_idx
                B_N_result_file['N_smiles'][B_N_result_idx] = n_smiles
                B_N_result_file['N_Index'][B_N_result_idx] = nu_mol_idx
                B_N_result_file['N_Atomid'][B_N_result_idx] = nu_atom_idx
                b_energy_r = Brow['E_energy_r']
                n_energy_r = Nrow['E_energy_r']
                b_g_r = Brow['G_energy_r']
                n_g_r = Nrow['G_energy_r']
            if nu_mol_idx not in duplicate_N_id:
                dft_dir_ = os.path.join(dft_dir, f"B_N_{type_}")
                spe_dir_ = os.path.join(spe_dir, f"B_N_{type_}")
                identify_str = f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_{type_}"
            else:
                dft_dir_ = os.path.join(dft_dir, f"B_N_{type_}_d")
                spe_dir_ = os.path.join(spe_dir, f"B_N_{type_}_d")
                identify_str = f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_Naid_{nu_atom_idx:05}_{type_}"
            mol_file = glob.glob(os.path.join(mol_dir, f'{identify_str}.mol'))
            if len(mol_file) == 0:
                print(b_mol_idx, nu_mol_idx, "OPT Is Error")
                # return None
            mol_file = mol_file[0]
            opt_files = glob.glob(os.path.join(dft_dir_, f'{identify_str}*.log'))
            temp_idx, temp_E, temp_G = [], [], []
            for opt_file in opt_files:
                opt = mol_manipulation.logfile_process.Logfile(opt_file)
                spe_files = glob.glob(os.path.join(spe_dir_, os.path.split(opt_file)[-1]))
                if len(spe_files) == 0:
                    print(b_mol_idx, nu_mol_idx, "SPE Error")
                    continue
                conf_id = int(opt_file.split('.')[0].split("_")[-1])
                spe_file = spe_files[0]
                spe = mol_manipulation.logfile_process.Logfile(spe_file)
                electric_energy = spe.all_engs[0]
                G_cor = opt.all_engs[-1]
                temp_idx.append(conf_id)
                temp_E.append(electric_energy)
                temp_G.append(G_cor + electric_energy)
            if len(temp_G) == 0 or len(opt_files) == 0:
                B_N_result_file['Error_reason'][B_N_result_idx] = "DFT_error"
            else:
                min_index = np.argmin(temp_G)
                B_N_result_file[f"E_energy_{type_}"][B_N_result_idx] = temp_E[min_index]
                B_N_result_file[f"G_energy_{type_}"][B_N_result_idx] = temp_G[min_index]
                B_N_result_file[f"conf_idxs_{type_}"][B_N_result_idx] = temp_idx[min_index]
                if type_ == 'r':
                    B_N_result_file["deltaE_comb(kcal)"][B_N_result_idx] = (temp_E[min_index] - b_energy_r - n_energy_r) * 627.5
                    B_N_result_file["deltaG_comb(kcal)"][B_N_result_idx] = (temp_G[min_index] - b_g_r - n_g_r) * 627.5
                else:
                    B_N_result_file["deltaE_react"][B_N_result_idx] = (temp_E[min_index] - B_N_result_file["E_energy_r"][B_N_result_idx])
                    B_N_result_file["deltaG_react"][B_N_result_idx] = (temp_G[min_index] - B_N_result_file["G_energy_r"][B_N_result_idx])
            B_N_result_idx += 1
    pd.DataFrame(B_N_result_file).to_csv(B_N_result_path, index=False)

def reaction_calc_ts(target_dir, method, om_name="OM", ts_name="TS"):
    om_file_dir = target_dir + "/" + om_name
    ts_dir = target_dir + "/" + ts_name
    if not os.path.isdir(ts_dir):
        os.mkdir(ts_dir)
    om_log_files = glob.glob(om_file_dir + "/*.log")
    for om_log_file in om_log_files:
        try:
            om_log = logfile_process.Logfile(om_log_file)
            # assert om_log.bond_attach
            mol_manipulation.om_to_ts(om_log, 1, new_dir=ts_dir, method=method)     
        except:
            continue     

def reaction_calc_irc(target_dir, ts_name='ts', irc_name='irc', freq_limit = -100, require_all = False):
    # 仅针对频率较低或者振动方向错误的
    ts_file_dir = target_dir + "/" + ts_name
    irc_dir = target_dir + "/" + irc_name
    if not os.path.isdir(irc_dir):
        os.mkdir(irc_dir)
    ts_log_files = glob.glob(ts_file_dir + "/*.log")
    for ts_log_file in ts_log_files:
        print(ts_log_file, end='\r')
        ts_log = logfile_process.Logfile(ts_log_file)
        assert ts_log.bond_attach
        title = [each - 1 for each in ts_log.title]
        position = ts_log.running_positions[-1]
        bond = mol_manipulation.Tool.get_atoms_distance(position[title[0]], position[title[1]])
        if bond > 3.0 or bond < 2.0 or require_all:
            print("%s May have wrong B_Cl bond!!! " % ts_log.file_dir)
        else:
            # continue
            if ts_log.is_right_ts:
                if float(ts_log.first_unreal_freq) <= freq_limit:
                    continue
                else:
                    print("%s May have wrong vibration freq: %.4f!!! " % (ts_log.file_dir, float(ts_log.first_unreal_freq)))
            else:
                print("%s May have wrong vibration direction!!! " % ts_log.file_dir)
        mol_manipulation.ts_to_irc(ts_log, new_dir=irc_dir)

def calc_distribution2(y, eachsize=0.01, title=None, xlab=None, ylab="Count", y_max=None, y_min=None, figure_size = (4,3)):
    if y_max == None:    y_max = np.max(y)
    if y_min == None:    y_min = np.min(y)
    X = np.arange(y_min, y_max + eachsize, eachsize)
    des = [0 for each in X]
    z = (y - y_min)/eachsize
    for each in z:
        try:
            des[int(each)] += 1
        except:
            continue
    des = np.array(des)
    # des = des / len(y)
    
    fig = plt.figure(figsize=figure_size)
    ax = fig.add_subplot(111)
    ax.patch.set_alpha(0.0)
    plt.bar(X, des, width=eachsize/2, color="green")
    plt.xlim(y_min - eachsize, y_max + eachsize)
    plt.ylim(0, np.max(des) * 1.2)
    plt.xlabel(xlab, fontsize=30)
    plt.ylabel(ylab, fontsize=30)
    plt.xticks(fontsize=15)
    plt.yticks(fontsize=15)
    if title != None:
        plt.title = title
    plt.tight_layout()
    plt.savefig('test.svg', format='svg')
    plt.show()
    return des 


# Model 
def normalize_axis(arr, axis=0, mean=[], std=[]):
    """
    对数组中的某一维进行标准化（z-score normalization）
    
    参数：
    arr: ndarray，输入的数组
    axis: int，标准化的维度
    
    返回值：
    normalized_arr: ndarray，标准化后的数组
    """
    if len(mean) == 0 or len(mean) == 0:
        mean = np.mean(arr, axis=axis, keepdims=True)  # 计算均值
        std = np.std(arr, axis=axis, keepdims=True)  # 计算标准差
    normalized_arr = (arr - mean) / std  # 标准化
    normalized_arr = np.nan_to_num(normalized_arr, 0)
    return normalized_arr, mean, std

def half_ood_folds(react_data, n_folds=3, seed=0, ignore_B=0, ignore_N = 0, ignore_Cl = 0):
    np.random.seed(seed)
    all_B_ids = np.unique(react_data["B_Index"].to_numpy())
    all_N_ids = np.unique(react_data["N_Index"].to_numpy())
    all_Cl_ids = np.unique(react_data["Cl_Index"].to_numpy())
    Blist, Nlist, Cllist = np.arange(len(all_B_ids)), np.arange(len(all_N_ids)), np.arange(len(all_Cl_ids))
    np.random.shuffle(Blist)
    np.random.shuffle(Nlist)
    np.random.shuffle(Cllist)
    B_random_list ={id: int(each) % n_folds for id, each in zip(all_B_ids, Blist)}
    N_random_list ={id: int(each) % n_folds for id, each in zip(all_N_ids, Nlist)}
    Cl_random_list ={id: int(each) % n_folds for id, each in zip(all_Cl_ids, Cllist)}
    folds = [[[], []] for _ in range(n_folds)]
    for id, B_Index in enumerate(react_data["B_Index"]):
        N_Index = react_data['N_Index'][id]
        Cl_Index = react_data['Cl_Index'][id]
        B_idx, N_idx, Cl_idx = B_random_list[B_Index], N_random_list[N_Index], Cl_random_list[Cl_Index]
        for idx in range(n_folds):
            if (ignore_B or idx == B_idx) and (ignore_N or idx == N_idx) and (ignore_Cl or idx == Cl_idx):
                folds[idx][1].append(id)
                continue
            if (ignore_B or idx != B_idx) and (ignore_N or idx != N_idx) and (ignore_Cl or idx != Cl_idx):
                folds[idx][0].append(id)
    return folds

def draw_heatmap(x_labels, y_labels, values, title="None", figure_size=(40, 6), min_value = 0.0, max_value = 1):
    import seaborn as sns
    sns.set()
    # desc_labels = ["rdkit_mf", "morgan_mf", "rdkit_des", "modred_des"]
    # model_labels = ["GB", "XGB", "RF", "ET", "AdaB", "Line", "MLP"]
    # model_labels = ["no_product", "with_product", "with_structure"]

    plt.rcParams['font.sans-serif']='Arial'#设置中文显示，必须放在sns.set之后

    uniform_data = values #设置二维矩阵
    f, ax = plt.subplots(figsize=figure_size)
    annot_kws = {"fontsize": 30}
    #heatmap后第一个参数是显示值,vmin和vmax可设置右侧刻度条的范围,
    #参数annot=True表示在对应模块中注释值
    # 参数linewidths是控制网格间间隔
    #参数cbar是否显示右侧颜色条，默认显示，设置为None时不显示
    #参数cmap可调控热图颜色，具体颜色种类参考：https://blog.csdn.net/ztf312/article/details/102474190
    sns.heatmap(uniform_data, ax=ax,vmin=min_value,vmax=max_value,cmap='Blues',linewidths=2,cbar=True, annot=True,annot_kws=annot_kws, fmt='.3f')

    ax.set_title(title, fontsize=40) #plt.title('热图'),均可设置图片标题
    # ax.set_ylabel('descriptor', fontsize=10)  #设置纵轴标签
    # ax.set_xlabel('model', fontsize=10)  #设置横轴标签
    ax.set_xticklabels(x_labels, fontsize=30)
    ax.set_yticklabels(y_labels, fontsize=30)
    # #设置坐标字体方向，通过rotation参数可以调节旋转角度
    label_y =  ax.get_yticklabels()
    plt.setp(label_y, rotation=0, horizontalalignment='right')
    label_x =  ax.get_xticklabels()
    plt.setp(label_x, rotation=0, horizontalalignment='center')
    plt.savefig('test.svg', format='svg')
    plt.show()
    return plt

VAN_DER_WAALS_RADII = {"H": 1.20, "C":1.77, "N": 1.55, "O": 1.50, "F": 1.46, "S": 1.89, "Cl": 1.82, "Si":2.10, "P":1.80}
def Calc_cylinder_bv(symbol_lists, geom, Base_id, Cl_id, radius):
    # 非均匀格点积分
    next_atom = [each for each in range(len(symbol_lists)) if each not in [Base_id, Cl_id]][0]
    distance = np.linalg.norm(geom[Cl_id] - geom[Base_id])
    new_geom_position = mol_manipulation.trfm_rot(geom[Base_id], geom[Cl_id], geom[next_atom], geom, np.array([-distance/2,0,0]))
    new_geom_position = new_geom_position @ mol_manipulation.rotation(np.array([0,1,0]), 1, 0)

    num = 20
    cube_length = radius
    total_points = num * num * num 
    counts = np.zeros(8, dtype=np.int32)

    # 生成均匀的网格点
    x = np.linspace(0.1 -radius, radius - 0.1, num)
    y = x; z = np.linspace(0.1 -distance, distance - 0.1, num)
    points = np.array(np.meshgrid(x, y, z)).T.reshape(-1, 3)
    r = points[:,0] ** 2 + points[:,1] ** 2 <= radius ** 2
    points = points[r]
    points_inside = np.zeros(len(points), dtype=bool)
    for atom_id, (symbol, sphere_center) in enumerate(zip(symbol_lists, new_geom_position)):
        if atom_id in [Base_id, Cl_id]:
            continue
        radii = VAN_DER_WAALS_RADII[symbol]
        translated_points = points - sphere_center[:3]
        distances = np.linalg.norm(translated_points, axis=1)
        points_inside = points_inside | (distances <= radii)
    points_inside_z = points_inside
    return sum(points_inside_z) / len(points)

def shortest_distance(coords):
    """
    计算给定分子中任意两个原子之间的最短距离。

    参数：
    coords: N×3 的 numpy 数组，表示所有原子的坐标。

    返回值：
    分子中任意两个原子之间的最短距离。
    """
    # 计算原子间距离矩阵
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    dist_mat = np.linalg.norm(diff, axis=-1)

    # 忽略对角线元素并取上三角部分
    upper_dist_mat = np.triu(dist_mat, k=1)

    # 找到最小的非零距离值作为最短距离
    min_dist = np.amin(upper_dist_mat[upper_dist_mat > 0])

    return min_dist



def generate_ts_structure(row, model1, model2, B_N_des_map=None, Cl_des_map=None, 
reactant_dir=None,
duplicate_N_id = [9, 43, 285, 310, 314, 345, 346, 347, 348, 349, 350, 351, 352, 353, 354, 355, 356, 357, 358, 359, 360, 361, 362, 372, 375, 376], 
duplicate_Cl_id = [624, 625, 626, 627, 628, 629, 630, 631, 632, 633, 634, 635, 636, 637, 638, 639, 640, 642, 644, 645],
pre_position = [],
):
    def modify_ts_mol(b_nu_product_mol, cl_reactant_mol, distance_B_Cl, distance_C_Cl, bn_cl_dist, c_cl_dist):
        new_cl_reactant_mol = mol_manipulation.move_mol(cl_reactant_mol, np.array([distance_B_Cl + distance_C_Cl - bn_cl_dist / 2 - c_cl_dist / 2, 0, 0]))
        sum_mol = Chem.CombineMols(b_nu_product_mol, new_cl_reactant_mol)
        sum_position = sum_mol.GetConformer(0).GetPositions()
        sum_position[nb_Cl_atom_id] = np.array([distance_B_Cl - bn_cl_dist / 2, 0, 0])
        return sum_position, sum_mol
    dft_dir = os.path.join(reactant_dir, 'GS_OPT')
    Cl_r_file = os.path.join(dft_dir, 'Cl_r')
    B_N_p_file = os.path.join(dft_dir, 'B_N_p')
    B_N_p_d_file = os.path.join(dft_dir, 'B_N_p_d')
    mol_dir = os.path.join(reactant_dir, 'Mols')
    b_mol_idx = row['B_Index']
    nu_mol_idx = row['N_Index']
    cl_mol_idx = row['Cl_Index']
    b_atom_idx = row['B_Atomid']
    nu_atom_idx = row['N_Atomid']
    B_smiles = row['B_smiles']
    N_smiles = row['N_smiles']
    cl_atom_idx = int(row['Cl_Atomid'])
    bn_conf_idx = int(row['B_N_Cl_conf'])
    cl_conf_idx = int(row['Cl_r_conf'])
    react_name = f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_Cl_{cl_mol_idx:05}"
    if cl_mol_idx in duplicate_Cl_id:
        react_name = f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_Cl_{cl_mol_idx:05}_Claid_{cl_atom_idx:05}"
    B_N_name = f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}"
    Cl_name = f"Cl_{cl_mol_idx:05}"
    if cl_mol_idx in duplicate_Cl_id:
        Cl_name = f"Cl_{cl_mol_idx:05}_Claid_{cl_atom_idx:05}"

    if len(pre_position) != 0:
        distance_B_Cl = model1
        distance_C_Cl = model2
    else: 
        des_a = B_N_des_map[B_N_name]
        des_b = Cl_des_map[Cl_name]
        all_Xs = np.array([des_a + des_b])
        delta_G_r = 627.5*(all_Xs[:, 1] + all_Xs[:, 49])
        all_Xs = np.hstack([all_Xs, delta_G_r.reshape(-1, 1)])
        des = all_Xs[0]
        distance_B_Cl = model1.predict(des)
        distance_C_Cl = model2.predict(des)

    # print(distance_B_Cl, distance_C_Cl)
    # if distance_B_Cl > 2.8: distance_B_Cl = 2.8
    # if distance_B_Cl < 2.2: distance_B_Cl = 2.2
    # if distance_C_Cl > 2.2: distance_C_Cl = 2.2
    # if distance_C_Cl < 2.0: distance_C_Cl = 2.0
    # Generate B_N product 
    B_mol = mol_manipulation.smiles2mol(B_smiles)
    N_mol = mol_manipulation.smiles2mol(N_smiles)
    Chem.Kekulize(B_mol, clearAromaticFlags=True)
    Chem.Kekulize(N_mol, clearAromaticFlags=True)
    b_nu_mol = Chem.CombineMols(B_mol, N_mol)
    H_atom_idx = [each.GetIdx() for each in B_mol.GetAtomWithIdx(int(b_atom_idx)).GetNeighbors() if each.GetSymbol() == 'H'][0]
    b_nu_mol.GetAtomWithIdx(H_atom_idx).SetAtomicNum(17)
    b_nu_mol.GetAtomWithIdx(int(b_atom_idx)).SetFormalCharge(-1)
    chg = b_nu_mol.GetAtomWithIdx(int(nu_atom_idx + B_mol.GetNumAtoms())).GetFormalCharge()
    b_nu_mol.GetAtomWithIdx(int(nu_atom_idx + B_mol.GetNumAtoms())).SetFormalCharge(chg + 1)
    rwmol = Chem.RWMol(b_nu_mol)
    rwmol.AddBond(int(b_atom_idx), int(nu_atom_idx + B_mol.GetNumAtoms()), Chem.BondType.SINGLE)
    b_nu_product_mol = rwmol.GetMol()
    if nu_mol_idx in duplicate_N_id:
        # b_nu_product_mol = Chem.MolFromMolFile(os.path.join(mol_dir, f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_Naid_{nu_atom_idx:05}_p.mol"), removeHs=False, sanitize=False)
        b_nu_product_logs = glob.glob(os.path.join(B_N_p_d_file, f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_Naid_{nu_atom_idx:05}_p_{bn_conf_idx:04}.log"))
    else:
        # b_nu_product_mol = Chem.MolFromMolFile(os.path.join(mol_dir, f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_p.mol"), removeHs=False, sanitize=False)
        b_nu_product_logs = glob.glob(os.path.join(B_N_p_file, f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_p_{bn_conf_idx:04}.log"))
    if len(b_nu_product_logs) == 0:
        print(f"error of BN {b_mol_idx, nu_mol_idx}")
        return [None] * 5
    b_nu_product_log = mol_manipulation.logfile_process.Logfile(b_nu_product_logs[0])
    b_nu_atom_list, b_nu_position = b_nu_product_log.symbol_list, b_nu_product_log.running_positions[-1]

    cl_reactant_mol = Chem.MolFromMolFile(os.path.join(mol_dir, f"Cl_{cl_mol_idx:05}_r.mol"), removeHs=False)
    Chem.Kekulize(cl_reactant_mol, clearAromaticFlags=True)
    cl_reactant_logs = glob.glob(os.path.join(Cl_r_file, f"Cl_{cl_mol_idx:05}_r_{cl_conf_idx:04}.log"))
    if len(cl_reactant_logs) == 0:
        return [None] * 5
    cl_reactant_log = mol_manipulation.logfile_process.Logfile(cl_reactant_logs[0])
    cl_atom_list, cl_position = cl_reactant_log.symbol_list, cl_reactant_log.running_positions[-1]

    nb_Cl_atom_id = [each for each in b_nu_product_mol.GetAtomWithIdx(int(b_atom_idx)).GetNeighbors() if each.GetSymbol() == "Cl"][0].GetIdx()
    nb_another_atom_id = [each for each in b_nu_product_mol.GetAtomWithIdx(int(b_atom_idx)).GetNeighbors() if each.GetIdx() != nb_Cl_atom_id][0].GetIdx()
    
    cl_c_atom_id = cl_reactant_mol.GetAtomWithIdx(int(cl_atom_idx)).GetNeighbors()[0].GetIdx()
    cl_another_atom_id = [each for each in cl_reactant_mol.GetAtomWithIdx(int(cl_c_atom_id)).GetNeighbors() if each.GetIdx() != cl_atom_idx][0].GetIdx()
    
    if cl_atom_idx < cl_c_atom_id:
        restrict = [[b_atom_idx + 1, nb_Cl_atom_id + 1], [nb_Cl_atom_id + 1, b_nu_product_mol.GetNumAtoms() + cl_c_atom_id], [b_atom_idx + 1, b_nu_product_mol.GetNumAtoms() + cl_c_atom_id]]
    else:
        restrict = [[b_atom_idx + 1, nb_Cl_atom_id + 1], [nb_Cl_atom_id + 1, b_nu_product_mol.GetNumAtoms() + cl_c_atom_id + 1], [b_atom_idx + 1, b_nu_product_mol.GetNumAtoms() + cl_c_atom_id + 1]]

    if len(pre_position) != 0:
        assert len(pre_position) == b_nu_product_mol.GetNumAtoms() + cl_reactant_mol.GetNumAtoms() - 1
        b_nu_position = pre_position[:b_nu_product_mol.GetNumAtoms()]
        cl_atom_postion = pre_position[nb_Cl_atom_id]
        start_id = 0
        for i, pos in enumerate(cl_position):
            if i == cl_atom_idx:
                start_id = 1
                cl_position[i] = cl_atom_postion
            else:
                cl_position[i] = pre_position[b_nu_product_mol.GetNumAtoms() + i - start_id]
        b_nu_position = np.array(b_nu_position) + 1
        cl_position = np.array(cl_position) + 1

    new_nb_position = mol_manipulation.trfm_rot(b_nu_position[b_atom_idx], b_nu_position[nb_Cl_atom_id], b_nu_position[nb_another_atom_id], b_nu_position)
    new_cl_position = mol_manipulation.trfm_rot(cl_position[cl_atom_idx], cl_position[cl_c_atom_id], cl_position[cl_another_atom_id], cl_position)
    b_nu_product_mol = mol_manipulation.xtb_process.xtb_to_mol(b_nu_product_mol, [b_nu_atom_list], [new_nb_position], 1)
    cl_reactant_mol = mol_manipulation.xtb_process.xtb_to_mol(cl_reactant_mol, [cl_atom_list], [new_cl_position], 1)
    bn_cl_dist = mol_manipulation.Tool.get_atoms_distance(b_nu_position[b_atom_idx], b_nu_position[nb_Cl_atom_id])
    c_cl_dist = mol_manipulation.Tool.get_atoms_distance(cl_position[cl_c_atom_id], cl_position[cl_atom_idx])
    # for distance_B_Cl in np.arange(2.0, 2.9, 0.05):
    #     distance_C_Cl = 4.6 - distance_B_Cl
    sum_position, sum_mol = modify_ts_mol(b_nu_product_mol, cl_reactant_mol, distance_B_Cl, distance_C_Cl, bn_cl_dist, c_cl_dist)
    sum_position = np.delete(sum_position, b_nu_product_mol.GetNumAtoms() + cl_atom_idx, axis=0)
    # if shortest_distance(sum_position) < 0.8:
    #     print(f"{b_mol_idx:05} {nu_mol_idx:05} {cl_mol_idx:05} may have atoms too close ")
    sum_atom_list = b_nu_atom_list + cl_atom_list[:cl_atom_idx] + cl_atom_list[cl_atom_idx + 1:]
    title = " ".join([str(each) for each in [restrict[0][0], restrict[0][1], restrict[1][1]]])
    # B_N_Cl.FormatConverter.block_to_gjf(sum_atom_list, sum_position, os.path.join(ts_dir_, f"B_{b_mol_idx:05}_Nu_{nu_mol_idx:05}_Cl_{cl_mol_idx:05}.gjf"), charge=0, multiplicity=2, title=title,method=TS_METHOD)
    # break
    sum_rwmol = Chem.RWMol(sum_mol)
    sum_rwmol.RemoveBond(int(b_atom_idx), nb_Cl_atom_id)
    sum_rwmol.RemoveAtom(int(b_nu_product_mol.GetNumAtoms() + cl_atom_idx))
    sum_nwmol = sum_rwmol.GetMol()
    sum_nwmol = mol_manipulation.xtb_process.xtb_to_mol(sum_nwmol, [sum_atom_list], [sum_position], 1)
    group_a = list(range(b_nu_product_mol.GetNumAtoms()))
    group_a.remove(nb_Cl_atom_id)
    group_a.remove(b_atom_idx)
    group_b = list(range(b_nu_product_mol.GetNumAtoms(), sum_nwmol.GetNumAtoms()))
    return sum_nwmol, react_name, restrict, group_a, group_b

def generate_ts_structure_simple(row, bond_b_cl, bond_c_cl, title, atoms, pre_position,):
    B_smiles = row['B_smiles']
    N_smiles = row['N_smiles']
    # Generate B_N product 
    B_mol = mol_manipulation.smiles2mol(B_smiles)
    N_mol = mol_manipulation.smiles2mol(N_smiles)
    B_num_atoms = B_mol.GetNumAtoms() + N_mol.GetNumAtoms() - 1
    assert atoms[title[0]] == "B" and atoms[title[1]] == "Cl" and atoms[title[2]] == "C"
    Bs1_position = pre_position[:title[1]]
    Bs2_position = pre_position[title[1] + 1:B_num_atoms]
    Cs_position = pre_position[B_num_atoms:]
    B_position = pre_position[title[0]]
    Cl_position = pre_position[title[1]]
    C_position = pre_position[title[2]]
    this_B_Cl = mol_manipulation.Tool.get_atoms_distance(B_position, Cl_position)
    this_C_Cl = mol_manipulation.Tool.get_atoms_distance(C_position, Cl_position)
    move_B_Cl_array = (B_position - Cl_position) / this_B_Cl * (bond_b_cl - this_B_Cl)
    move_C_Cl_array = (C_position - Cl_position) / this_C_Cl * (bond_c_cl - this_C_Cl)
    new_Bs1_position = Bs1_position + move_B_Cl_array
    new_Bs2_position = Bs2_position + move_B_Cl_array
    new_Cs_position = Cs_position + move_C_Cl_array
    new_position = pre_position
    new_position[:title[1]] = new_Bs1_position
    new_position[title[1] + 1:B_num_atoms] = new_Bs2_position
    new_position[B_num_atoms:] = new_Cs_position
    return new_position
    
def detail_irc(target_dir, old_irc_name="irc", new_irc_name="IRC_detailed", maxpoints=30):
    old_irc_dir = os.path.join(target_dir, old_irc_name)
    new_irc_dir = os.path.join(target_dir, new_irc_name)
    if not os.path.isdir(new_irc_dir):
        os.mkdir(new_irc_dir)
    logs = glob.glob(os.path.join(old_irc_dir, "*.log"))
    for log_file in logs:
        name = os.path.split(log_file)[-1].split(".")[0]
        log = mol_manipulation.logfile_process.Logfile(log_file)
        title = [each - 1 for each in log.title]
        distance_1 = Tool.get_atoms_distance(log.running_positions[-1][title[0]], log.running_positions[-1][title[1]])
        distance_2 = Tool.get_atoms_distance(log.running_positions[-1][title[1]], log.running_positions[-1][title[2]])
        old_distance_1 = Tool.get_atoms_distance(log.running_positions[0][title[0]], log.running_positions[0][title[1]])
        old_distance_2 = Tool.get_atoms_distance(log.running_positions[0][title[1]], log.running_positions[0][title[2]])
        if (distance_1 < 2.1 and distance_2 > 2.1):
            continue
        if (distance_2 < 2.1 and distance_1 > 2.1):
            continue
        if np.abs(distance_1 - old_distance_1) > 0.3 or np.abs(distance_2 - old_distance_2) > 0.3:
            new_method = log.method.replace('maxpoints=15', f'maxpoints={maxpoints}')
            new_log = os.path.join(new_irc_dir, f"{name}.gjf")
            FormatConverter.block_to_gjf(log.symbol_list, log.first_atom_position, new_log, charge=log.charge, multiplicity=log.multiplicity, title=" ".join([str(each) for each in log.title]), method=new_method)
            

def check_irc(target_dir, ts_name = "TS", irc_name = "IRC", ts_fail_name="ts_irc_fail_new"):
    ts_dir_ = os.path.join(target_dir, ts_name)
    irc_dir_ = os.path.join(target_dir, irc_name)
    irc_fail_dir = os.path.join(target_dir, ts_fail_name)
    if not os.path.isdir(irc_fail_dir):
        os.mkdir(irc_fail_dir)
    logs = glob.glob(os.path.join(ts_dir_, "*.log"))
    for log_file in logs:
        pass_a, pass_b = 0,0
        # try:
        name =os.path.split(log_file)[-1].split(".")[0]
        gjfs = glob.glob(os.path.join(irc_dir_, f"{name}*.gjf"))
        if len(gjfs) == 0:
            continue
        irc_log_files = glob.glob(os.path.join(irc_dir_, f"{name}*.log"))
        # if len(irc_log_files) == 0:
        #     continue
        for irc_log_file in irc_log_files:
            log = mol_manipulation.logfile_process.Logfile(irc_log_file)
            title = [each - 1 for each in log.title]
            position = log.running_positions[-1]
            distance_1 = Tool.get_atoms_distance(log.running_positions[-1][title[0]], log.running_positions[-1][title[1]])
            distance_2 = Tool.get_atoms_distance(log.running_positions[-1][title[1]], log.running_positions[-1][title[2]])
            old_distance_1 = Tool.get_atoms_distance(log.running_positions[0][title[0]], log.running_positions[0][title[1]])
            old_distance_2 = Tool.get_atoms_distance(log.running_positions[0][title[1]], log.running_positions[0][title[2]])
            # if distance_1 < 2.1 and distance_2 > 2.1:
            #     pass_a = 1
            # elif distance_2 < 2.1 and distance_1 > 2.1:
            #     pass_b = 1
            if distance_1 < 2.1:
                pass_a = 1
            elif distance_2 < 2.1:
                pass_b = 1
        # except:
        #     pass
        if not pass_a or not pass_b:
            new_file = os.path.join(irc_fail_dir, os.path.split(log_file)[-1])
            shutil.move(log_file, new_file)