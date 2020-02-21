#!/usr/bin/env python3
from itertools import chain

import numpy as np

from rdkit.Chem import MACCSkeys
from rdkit import DataStructs
from rdkit.DataStructs import cDataStructs as cds

from fiprofiling import chem, util

def rdmol2maccs(rdmol):
    return MACCSkeys.GenMACCSKeys(rdmol)

def boollist2fp(boollist):
    fp = cds.ExplicitBitVect(len(boollist))
    for i, boolvalue in enumerate(boollist):
        if boolvalue:
            fp.SetBit(i)
    return fp

def rdmol_rdpatterns2fp(rdmol, patternvector):
    if not rdmol:
        return None
    boollist = [chem.rdmol_has_substruct_pattern(rdmol, pattern)
                for pattern in patternvector]
    return boollist2fp(boollist)

def bintext2fp(bintext):
    return cds.CreateFromBitString(bintext)

def fp2bintext(fp):
    return cds.BitVectToText(fp)

def hextext2fp(hextext):
    return cds.CreateFromFPSText(hextext)

def fp2hextext(fp):
    return cds.BitVectToFPSText(fp)

def fp_tanimoto_similarity(fp1, fp2):
    return DataStructs.FingerprintSimilarity(fp1, fp2,
                                             metric=DataStructs.TanimotoSimilarity)

def fp2nparray(fp):
    nparray = np.zeros((1,), dtype=np.bool_)
    cds.ConvertToNumpyArray(fp, nparray)
    return nparray

def fp2onbits(fp):
    return tuple(fp.GetOnBits())

def rdmols_smarts2fps(mols, smartsvector):
    patternvector = tuple(chem.smarts2rdmol(smarts) for smarts in smartsvector)
    for mol in mols:
        yield rdmol_rdpatterns2fp(mol, patternvector)

def smiles_smarts2fps(smiles, smartsvector):
    rdmols = (chem.smiles2rdmol(smile) for smile in smiles) # :)
    for fp in rdmols_smarts2fps(rdmols, smartsvector):
        yield fp

def fps2avg_tanimoto(fps):
    fp_count = len(fps)
    pair_count = 0.5*(fp_count**2 - fp_count)
    tanimoto_sum = 0.0
    for i, fp1 in enumerate(fps):
        for fp2 in fps[i+1:]:
            tanimoto_sum += fp_tanimoto_similarity(fp1, fp2)
    return tanimoto_sum / pair_count

def fps2avg_on_bits(fps):
    onbits = 0.0
    fpcount = 0
    for fp in fps:
        fpcount += 1
        onbits += len(fp2onbits(fp))
    return onbits/fpcount

def fps2avg_bit_value(fps):
    # Assumes same size of all fps
    fpiter = iter(fps)
    firstfp = next(fpiter)
    avg_on_bits = fps2avg_on_bits(chain([firstfp], fpiter))
    return avg_on_bits/len(firstfp)
