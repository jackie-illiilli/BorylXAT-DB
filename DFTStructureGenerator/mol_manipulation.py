import os, glob, math, shutil, pickle, random, copy
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import MolDrawing, DrawingOptions
from rdkit.Geometry import Point3D
from . import FormatConverter, xtb_process
from . import Tool, logfile_process, B_N_Cl
op_dir = "../files/"

DrawingOptions.bondLineWidth = 1.8
DrawingOptions.atomLabelFontSize = 14
Draw.DrawingOptions.includeAtomNumbers = True
atom_num = [None, "H", "He", "Li", "Be", "B", "C", "N", "O", "F",
            "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
            "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr"]



def move(a):
    """Convert a 3D rotation-translation matrix 'a' into a 4D rotation-translation matrix

    Args:
        a (iterable): 

    Returns:
        array: size = 4*4
    """    
    x, y, z = a[:3]
    return np.array([[1, 0, 0, x], [0, 1, 0, y], [0, 0, 1, z], [0, 0, 0, 1]])

def rotation(a, sin, cos):
    """Generate a 4D matrix that, when multiplied, rotates along the 'a' axis by sin and cos angles

    Args:
        a (array): 3D
        sin (float): 
        cos (float): 

    Returns:
        array: size 4*4
    """    
    a = np.array(a)
    a = a / np.sqrt(a @ a.T)
    u, v, w = a[:3]
    return np.array([[u * u + (1 - u * u) * cos, u * v * (1 - cos) - w * sin, u * w * (1 - cos) + v * sin, 0],
                     [u * v * (1 - cos) + w * sin, v * v + (1 - v * v)
                      * cos, v * w * (1 - cos) - u * sin, 0],
                     [u * w * (1 - cos) - v * sin, v * w * (1 - cos) +
                      u * sin, w * w + (1 - w * w) * cos, 0],
                     [0, 0, 0, 1]])

def trfm_rot(a, b, c, position=[], center_point=np.array([0, 0, 0])):
    """For points a, b, and c, perform translation and rotation so that the ab vector is parallel to the x-axis, a, b, and c are on the xoy plane, and a and b are symmetric about center_point

    Args:
        a (_type_): 3D coordinates, ultimately located in the negative x-axis direction
        b (_type_): 3D coordinates, ultimately located in the positive x-axis direction
        c (_type_): 3D coordinates, after adjustment its y-axis coordinate is positive and z-axis coordinate is 0
        position (list, optional): If position is given, outputs the transformed position; otherwise, outputs the rotation-translation matrix. Defaults to [].
        center_point (_type_, optional): _description_. Defaults to np.array([0, 0, 0]).

    Returns:
        array: new array (4*4)
        or tr_array, rot1_array, rot2_array, tr_array2
    """    
    zero_point = np.array([0, 0, 0])
    x, y, z = np.array(a), np.array(b), np.array(c)
    x_axis = np.array([1, 0, 0])
    mid_point = (x + y)/2
    tr_array_3d = zero_point - mid_point
    tr_array = move(tr_array_3d)
    x += tr_array_3d
    y += tr_array_3d
    z += tr_array_3d
    tr_array2 = move(center_point - zero_point)
    xy = y - x
    xy = xy / np.sqrt((xy @ xy))
    law_axis = np.cross(xy, x_axis)
    if (law_axis == np.array([0, 0, 0])).all() == False:
        law_axis /= np.sqrt((law_axis @ law_axis))
    sin = np.sqrt(xy[1] * xy[1] + xy[2] * xy[2])
    cos = xy[0]
    rot1_array = rotation(law_axis, sin, cos)
    oz = z - zero_point
    oz = rot1_array[:3, :3] @ oz
    oz[0] = 0
    oz = oz / np.sqrt((oz @ oz))
    sin = oz[2]
    cos = oz[1]
    rot2_array = rotation(-1 * x_axis, sin, cos)
    if len(position) == 0:
        return tr_array, rot1_array, rot2_array, tr_array2
    else:
        if len(position[0]) == 3:
            up_position = np.insert(position, 3, np.ones(len(position)), 1).T
        else:
            up_position = position.T
        up_position = tr_array @ up_position
        up_position = rot1_array @ up_position
        up_position = rot2_array @ up_position
        up_position = tr_array2 @ up_position
        up_position = up_position.T
    return up_position

def rot_mol(mol, axis=np.array([0, 1, 0]), sin=0, cos=-1):
    """Rotate all conformations of the molecule, default 180° around the y-axis

    Args:
        mol (_type_): _description_
        axis (array, optional): Rotation axis. Defaults to np.array([0, 1, 0]).
        sin (int, optional): _description_. Defaults to 0.
        cos (int, optional): _description_. Defaults to -1.

    Returns:
        mol: _description_
    """    
    # Define rotation corner is 180 degree
    react1 = copy.deepcopy(mol)
    for conformer in react1.GetConformers():
        position = conformer.GetPositions()
        up_position = np.insert(position, 3, np.ones(len(position)), 1).T
        matrix = rotation(axis, sin, cos)
        up_position = (matrix @ up_position).T
        for i, c in enumerate(up_position):
            conformer.SetAtomPosition(i, Point3D(c[0], c[1], c[2]))
    return react1


def move_mol(mol, array=np.array([0, 0, 1.5])):
    """Move all conformations of the molecule, default 1.5 units along the z-axis

    Args:
        mol (_type_): _description_
        array (array, optional): _description_. Defaults to np.array([0, 0, 1.5]).

    Returns:
        _type_: _description_
    """    
    react1 = copy.deepcopy(mol)
    for conformer in react1.GetConformers():
        position = conformer.GetPositions()
        up_position = np.insert(position, 3, np.ones(len(position)), 1).T
        matrix = move(array)
        up_position = (matrix @ up_position).T
        for i, c in enumerate(up_position):
            conformer.SetAtomPosition(i, Point3D(c[0], c[1], c[2]))
    return react1

def om_to_ts(log_file:logfile_process.Logfile, write_gjf=True, new_dir="ts",
    method="opt=(calcfc,ts,noeigen) freq b3lyp/6-31g* em=gd3bj"):
    """Convert Gaussian constrained optimization output file into transition state calculation input file or information

    Args:
        log_file (logfile_process.Logfile): log file
        write_gjf (bool, optional): Whether to write the input file directly. Defaults to True.
        new_dir (str, optional): Address for storing transition state input files. Defaults to "ts".
        method (str, optional): #p line. Defaults to "opt=(calcfc,ts,noeigen) freq b3lyp/6-31g* em=gd3bj".

    Returns:
        _type_: _description_
    """    
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)
    new_name = os.path.split(log_file.file_dir)[-1].split(".")[0] + ".gjf"
    newfile = new_dir + "/" + new_name 
    title = log_file.title
    charge = log_file.charge
    symbol_list = log_file.symbol_list
    multiplicity = log_file.multiplicity
    position = log_file.running_positions[-1]
    title = " ".join(str(each) for each in title)
    multiplicity = log_file.multiplicity
    if write_gjf:
        FormatConverter.block_to_gjf(symbol_list, position, newfile, charge, multiplicity, title,
                    method=method,
                    freeze=[])
    
    return title, symbol_list, position, charge, multiplicity


def ts_to_irc(log_file:logfile_process.Logfile, new_dir):
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)
    new_names = [os.path.split(log_file.file_dir)[-1].split(".")[0] + "forward.gjf",
    os.path.split(log_file.file_dir)[-1].split(".")[0] + "reverse.gjf"]
    methods = ['irc=(calcfc,stepsize=30,maxpoints=30,forward, lqa) b3lyp/6-31g(d) em=gd3bj scrf=(smd,solvent=toluene) nosymm', 
    'irc=(calcfc,stepsize=30,maxpoints=30,reverse, lqa) b3lyp/6-31g(d) em=gd3bj scrf=(smd,solvent=toluene) nosymm']

    title = log_file.title
    title = " ".join(str(each) for each in title)
    charge = log_file.charge
    multiplicity = log_file.multiplicity
    symbol_list = log_file.symbol_list
    position = log_file.running_positions[-1]
    for new_name, method in zip(new_names, methods):
        newfile = new_dir + "/" + new_name 
        FormatConverter.block_to_gjf(symbol_list, position, newfile, charge, multiplicity, title,
                    method=method,
                    freeze=[])
    
    return title, symbol_list, position, charge

# def log_to_eng(file_name, new_file_name=None,method="b3lyp/6-311+g(d,p) em=gd3bj"):
#     if new_file_name == None:
#         file_dir, name = os.path.split(file_name)
#         name = name.split(".")[0]
#         newname = name + "_eng"
#         newfile = file_dir + "/" + newname + ".gjf"
#     else:
#         newfile = new_file_name
#     # file_name = "om_test2_0"
#     title, charge, symbol_list, position= read_log(file_name, allow_unreal_freq=1)
#     title = " ".join(str(each) for each in title)
#     block_to_gjf(symbol_list, position, newfile, charge, title,
#                  method=method)

def smiles2mol(smiles, conf_num=20):
    """SMILES to mol conversion, including AddHs and Embed 3D structure generation

    Args:
        smiles (str): _description_

    Returns:
        mol: _description_
    """    
    mol = Chem.MolFromSmiles(smiles)
    if mol == None:
        print(smiles, "can't be read")
        return None
    Hmol = Chem.AddHs(mol)
    AllChem.EmbedMultipleConfs(Hmol, numConfs=conf_num, maxAttempts=100, )
    try:
        AllChem.MMFFOptimizeMolecule(Hmol)
    except:
        pass
    
    return Hmol


def error_improve(target_dir, mol_dir, file_name, dust_bin='dust_bin', improve_dir='improve', bond_attach_std='mol', maxcycles=0, method=None, bond_addition_function=None, bond_ignore_list=None, Inv_dir = 'Inv3', yqc_dir = 'yqc'):
    """Cooperate with log_process module to identify and modify Gaussian output files and handle imaginary frequency issues.
    Modified input files and error/imaginary frequency output files will be moved to specified folders.

    Args:
        target_dir (str): Root directory
        mol_dir (_type_): Directory corresponding to Mol molecules
        file_name (_type_): Folder under the root directory containing Gaussian output files
        dust_bin (str, optional): Trash bin for unsalvageably erroneous output files. Defaults to 'dust_bin'.
        improve_dir (str, optional): Storage path for modified Gaussian input files. Defaults to 'improve'.
    """    
    opt_file_dir = target_dir + "/" + file_name
    dust_bin_dir = target_dir + "/" + dust_bin
    improve_dir_ = target_dir + "/" + improve_dir
    Inv_dir_ = target_dir + "/" + Inv_dir
    yqc_dir_ = target_dir + "/" + yqc_dir
    mol_files = glob.glob(mol_dir + "/*.mol")
    
    for mol_file in mol_files:
        try:
            print("process %s" % mol_file, end='\r')
            log_files = glob.glob(opt_file_dir + "/" + os.path.split(mol_file)[-1].split(".")[0] + "*.log")
            if len(log_files) == 0:
                continue 
            for log_file in log_files:
                fail = 0
                with open(log_file, "r") as f:
                    lines = f.readlines()
                if len(lines) <= 15:
                    fail = 2
                else:
                    try:
                        opt_log= logfile_process.Logfile(log_file, mol_file_dir=mol_file, bond_attach_std=bond_attach_std, bond_addition_function=bond_addition_function, bond_ignore_list=bond_ignore_list)
                        if opt_log.multiplicity <= 0:# or opt_log.S_2 < 0:
                            fail = 1
                        elif not opt_log.bond_attach:
                            fail = 1
                        elif opt_log.file_type == "OM" and opt_log.unreal_freq == 0:
                            print("%s may not be a right OM for unreal freq num of %d" % (opt_log.file_dir, opt_log.unreal_freq))
                            # fail = 1
                        elif opt_log.file_type == "TS":
                            if opt_log.unreal_freq >= 0 and opt_log.unreal_freq != 1: # or not opt_log.is_right_ts:
                                print("%s is not a right TS for unreal freq num of %d" % (opt_log.file_dir, opt_log.unreal_freq))
                                # if opt_log.unreal_freq == 0:
                                fail = 1
                                # else:
                                #     opt_log.is_normal_end = 0
                                #     opt_log.error_reason = ""
                        if opt_log.file_type == "IRC":
                            if opt_log.irc_result == False:
                                fail = 1
                            else:
                                continue
                    except:
                        fail = 2 
                    
                if fail == 1:
                    new_log_name = dust_bin_dir + "/" + os.path.split(log_file)[-1] 
                    # new_log_name =  new_log_name.split(".")[0] + "%s.log" % opt_log.file_type
                    if not os.path.isdir(dust_bin_dir):
                        os.mkdir(dust_bin_dir) 
                    shutil.move(log_file, new_log_name)
                    continue
                if fail == 2:
                    new_log_name = improve_dir_ + "/" + os.path.split(log_file)[-1] 
                    if not os.path.isdir(improve_dir_):
                        os.mkdir(improve_dir_) 
                    shutil.move(log_file, new_log_name)
                    gjf_file = log_file.split(".")[0] + ".gjf"
                    new_gjf_name = improve_dir_ + "/" + os.path.split(gjf_file)[-1] 
                    shutil.move(gjf_file, new_gjf_name)
                    continue
                new_log_name = target_dir + '/' + improve_dir
                savechk = None
                readchk = None
                # if opt_log.file_type == "OM":
                #     savechk = os.path.split(log_file)[-1].split(".")[0]
                # if opt_log.file_type == "TS":
                #     readchk = os.path.split(log_file)[-1].split(".")[0]
                if not opt_log.normal_end:
                    opt_log.solve_error_logfile(new_log_name, Inv_dir=Inv_dir_, yqc_dir=yqc_dir_, savechk=savechk, readchk=readchk, maxcycles=maxcycles, method=method)
                elif opt_log.unreal_freq and opt_log.file_type not in ["OM", "TS"]:
                    opt_log.unreal_freq_improve(new_log_name, savechk=savechk, readchk=readchk, method=method)
        except:
            continue

