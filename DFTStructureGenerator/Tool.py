# Collection of common utility functions
from inspect import BoundArguments
import math, os, shutil
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolTransforms


def find_first_line(fileline, find_str, find_type="start"):
    assert find_type in ['start', 'end', 'all', 'in']
    if find_type == "start":
        line_function = lambda x, y: x.startswith(y)
    elif find_type == 'end':
        line_function = lambda x, y: x.endswith(y)
    elif find_type == 'all':
        line_function = lambda x, y: x == y
    else:
        line_function = lambda x, y: y in x

    for idx, line in enumerate(fileline):
        if line_function(line, find_str):
            return idx, line
    return None, None

def get_atoms_distance(atom_positionA, atom_positionB):
    """[summary]

    Args:
        atom_positionA (array): [description]
        atom_positionA (array): [description]

    Returns:
        array: distance
    """
    return np.sqrt(sum((atom_positionA - atom_positionB) ** 2))

def get_bond_angle(atom_positionA, atom_positionB, atom_positionC):
    conf = Chem.rdchem.Conformer(3)
    all_positions = [atom_positionA, atom_positionB, atom_positionC]
    for i in range(3):
        conf.SetAtomPosition(i, all_positions[i][:3])
    bond_angle = rdMolTransforms.GetAngleDeg(conf, 0, 1, 2)
    return bond_angle


def GetSpinMultiplicity(Mol, CheckMolProp = True):
    """From RDKitUtil.py
    Get spin multiplicity of a molecule. The spin multiplicity is either
    retrieved from 'SpinMultiplicity' molecule property or calculated from
    from the number of free radical electrons using Hund's rule of maximum
    multiplicity defined as 2S + 1 where S is the total electron spin. The
    total spin is 1/2 the number of free radical electrons in a molecule.

    Arguments:
        Mol (object): RDKit molecule object.
        CheckMolProp (bool): Check 'SpinMultiplicity' molecule property to
            retrieve spin multiplicity.

    Returns:
        int : Spin multiplicity.

    """
    
    Name = 'SpinMultiplicity'
    if (CheckMolProp and Mol.HasProp(Name)):
        return int(float(Mol.GetProp(Name)))

    # Calculate spin multiplicity using Hund's rule of maximum multiplicity...
    NumRadicalElectrons = 0
    for Atom in Mol.GetAtoms():
        NumRadicalElectrons += Atom.GetNumRadicalElectrons()

    TotalElectronicSpin = NumRadicalElectrons/2
    SpinMultiplicity = 2 * TotalElectronicSpin + 1
    
    return int(SpinMultiplicity)
