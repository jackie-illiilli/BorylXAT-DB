# Some file conversion scripts are written here
import numpy as np
import pandas as pd
import shutil, os

from rdkit import Chem
from rdkit.Chem import AllChem

from . import Tool


def read_gjf(gjf_file, read_details=False):
    with open(gjf_file, "rt") as f:
        lines = f.readlines()
    lines = [line.strip('\n') for line in lines]
    for id, line in enumerate(lines):
        if line.startswith("#"):
            break
    temp_line = lines[id]
    if lines[id + 1] != '':
        temp_line = lines[id] + lines[id + 1]
        id += 1
    if temp_line.startswith("#p"):
        method = temp_line.split("#p ")[1]
    else:
        method = temp_line.split("# ")[1]

    atoms = []
    positions = []
    title = lines[id + 2]
    charge, multi = [int(each) for each in lines[id + 4].split()]
    for line in lines[id + 5:]:
        gjf_list = line.split()
        if len(gjf_list) != 4:
            break
        atoms.append(gjf_list[0])
        positions.append([float(each) for each in gjf_list[1:]])
    if read_details:
        return atoms, positions, method, charge, multi, title
    return atoms, positions

def mol_to_xyz(mol, atom_list=None, position_list=None, file_dir="test.xyz",title=None):
    """The most basic process of converting mol to xyz file

    Args:
        mol (Chem.Mol): 
        file_dir (str, optional): _description_. Defaults to "test.xyz".
        title (str, optional): _description_. Defaults to None.
        

    Returns:
        filename: _description_
    """   
    if mol is None and (atom_list is None or position_list is None):
        print("mol, atom_list/position_list must contain one")
        return None 
    file_dir = file_dir.split(".")[0]
    if mol is not None:
        atom_num = mol.GetNumAtoms()
        atom_list = [atom.GetSymbol() for atom in mol.GetAtoms()]
        position_list = [conf.GetPositions() for conf in mol.GetConformers()]
    else:
        atom_num = len(atom_list)
        
    file_names = []
    for i, position in enumerate(position_list):
        file_name = "%s_%s.xyz" % (file_dir, str(i))
        with open(file_name, "wt") as f:
            f.write(" %d\n" % atom_num)
            f.write("%s\n" % title)
            for j, atom in enumerate(atom_list):
                f.write("%s %.8f %.8f %.8f\n" % (atom, position[j][0],position[j][1],position[j][2]))
        file_names.append(file_name)
    return file_names

def write_xyz_file(filename, molecules):
    """
    Writes the atom lists and positions of multiple molecules to an XYZ file.

    :param filename: output XYZ file name
    :param molecules: A list containing multiple molecules, each molecule is a dictionary, in the following format:
                      {
                          'name': 'molecule name',
                          'atomlist': ['H', 'O', 'C'], # atom type list
                          'positions': [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]] # atomcoordinates list
                      }
    """
    with open(filename, 'w') as f:
        # Write atom information for each molecule
        for molecule in molecules:
            atomlist = molecule['atomlist']
            positions = molecule['positions']
            num_atoms = len(atomlist)

            # Write the title row of the molecule
            f.write(f"{num_atoms}\n")
            f.write(f"{molecule['name']}\n")

            # Write atom information
            for atom, position in zip(atomlist, positions):
                f.write(f"{atom} {' '.join(map(str, position))}\n")

def mol_to_gjf(mol, file_="test_data/mol2gjf.gjf", charge=None, SpinMultiplicity=None, title="Title", method="opt freq b3lyp/6-311g(d,p)", confid=0, ignore_warning=False, freeze=[], difreeze=[], savechk=None, readchk=None, final_line=None):
    """cover mol to gjf, title line with charge

    Args:
        mol (Chem.Mol): _description_
        file_dir (str, optional): file_dir. Defaults to "test_data/mol2gjf.gjf".
        charge (int, optional): 0 or None. Defaults to None.
        title (str, optional): title. Defaults to "Title".
        method (str, optional): method use to calculation. Defaults to "opt freq b3lyp/6-311g(d,p)".
        confid (int, optional): use special conformers with mol. Defaults to 0.
    """    
    if charge == None:
        charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if SpinMultiplicity == None:
        SpinMultiplicity = Tool.GetSpinMultiplicity(mol)
    file_dir, filename = os.path.split(file_)
    if file_dir != "":
        if not os.path.isdir(file_dir):
            os.mkdir(file_dir)
    with open(file_, "wt") as f:
        if savechk != None:
            f.write("%%chk=%s.chk\n" % savechk)
        if readchk != None:
            f.write("%%oldchk=%s.chk\n" % readchk)
        f.write("%nprocshared=28\n%mem=56GB\n#p")
        f.write(" %s\n\n" % method)
        f.write("$$$$%s####%d????\n\n" % (title, charge))
        f.write("%d %d\n" % (int(charge), SpinMultiplicity))
        if SpinMultiplicity != 1 and not ignore_warning:
            print("%s 's SpinMultiplicity != 1, check it" % file_)
        if len(mol.GetConformers()) == 0:
            mol = AllChem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            AllChem.EmbedMultipleConfs(mol, numConfs=1)
        positions = mol.GetConformer(confid).GetPositions()
        for i, atom in enumerate(mol.GetAtoms()):
            f.write(" %s                 " % atom.GetSymbol())
            f.write("%f %f %f\n" % tuple(positions[i]))
        f.write("\n")
        for each in freeze:
            f.write("B %d %d F\n" % tuple(each))
        for each in difreeze:
            f.write("D %d %d %d %d F\n" % tuple(each))
        if final_line != None:
            f.write(final_line)
            f.write("\n")
        f.write("\n\n")

def block_to_gjf(symbol_list, positions, file="test_data/mol2gjf.gjf", charge=0, multiplicity=1, title="Title", method="opt freq b3lyp/6-311g(d,p)", freeze=[], difreeze=[], savechk=None, readchk=None, final_line=None):
    """cover symbol_list and position into gjf, can write link1, as om&ts

    Args:
        symbol_list (_type_): _description_
        positions (_type_): _description_
        file (str, optional): _description_. Defaults to "test_data/mol2gjf.gjf".
        charge (int, optional): _description_. Defaults to 0.
        title (str, optional): _description_. Defaults to "Title".
        method (str, optional): _description_. Defaults to "opt freq b3lyp/6-311g(d,p)".
        method2 (_type_, optional): _description_. Defaults to None.
        freeze (_type_, optional): Freeze atoms, as TS calculation. Defaults to [].
        savechk (str, optional): Name of Chkfile, omit .chk. Defaults to None.
    """    
    file_dir, filename = os.path.split(file)
    if not os.path.isdir(file_dir):
        os.mkdir(file_dir)
    assert len(symbol_list) == len(positions)
    with open(file, "wt") as f:
        if savechk != None:
            f.write("%%chk=%s.chk\n" % savechk)
        if readchk != None:
            f.write("%%oldchk=%s.chk\n" % readchk)
        f.write("%nprocshared=28\n%mem=56GB\n#p")
        f.write(" %s\n\n" % method)
        f.write("$$$$%s####%d????\n\n" % (title, int(charge)))
        f.write("%d %d\n" % (int(charge), int(multiplicity)))
        # positions = mol.GetConformer(0).GetPositions()
        for i, atom in enumerate(symbol_list):
            f.write(" %s                 " % atom)
            f.write("%f %f %f\n" % tuple(positions[i][:3]))
        f.write("\n")
        for each in freeze:
            f.write("B %d %d F\n" % tuple(each))
        for each in difreeze:
            f.write("D %d %d %d %d F\n" % tuple(each))
        if final_line != None:
            f.write(final_line)
            f.write("\n")
        f.write("\n\n")

def read_chg_file(file):
    return_dict = {"index":{}, "symbol":{}, "charge":{}}
    with open(file) as f:
        lines = f.readlines()
    for idx, line in enumerate(lines):
        line_list = line.split()
        assert len(line_list) == 5
        return_dict["index"][idx] = idx
        return_dict["symbol"][idx] = line_list[0]
        return_dict["charge"][idx] = float(line_list[-1])
    return pd.DataFrame(return_dict)


