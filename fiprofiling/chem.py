#!/usr/bin/env python3

from rdkit import Chem
from rdkit import DataStructs

def smiles2rdmol(smiles):
    return Chem.MolFromSmiles(smiles)

def smarts2rdmol(smarts):
    return Chem.MolFromSmarts(smarts)

def rdmol_has_substruct_pattern(rdmol, pattern):
    if not rdmol:
        return None
    return rdmol.HasSubstructMatch(pattern)

def sdf2rdmols(sdfpath):
    supplier = Chem.SDMolSupplier(sdfpath)
    for mol in supplier:
        if mol:
            yield mol
