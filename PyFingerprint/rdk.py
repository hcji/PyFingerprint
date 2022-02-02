# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 15:39:39 2018

@author: hcji
"""

from rdkit.Chem import MACCSkeys
from rdkit import Chem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Avalon import pyAvalonTools

def rdk_fingerprint(smi, fp_type="rdkit", size=1024):
    _fingerprinters = {"rdkit": Chem.rdmolops.RDKFingerprint
                           , "rdk-maccs": MACCSkeys.GenMACCSKeys
                           , "topological-torsion": Torsions.GetTopologicalTorsionFingerprint
                           , "avalon": pyAvalonTools.GetAvalonFP}
    mol = Chem.MolFromSmiles(smi)
    if fp_type in _fingerprinters:
        fingerprinter = _fingerprinters[fp_type]
        fp = fingerprinter(mol)
    elif fp_type == "atom-pair":
        fp = Pairs.GetAtomPairFingerprintAsBitVect(mol)
    elif fp_type == "morgan":
        fp = GetMorganFingerprintAsBitVect(mol, 2, nBits=size)
    else:
        raise IOError('invalid fingerprint type')
    try:
        size = fp.GetNumBits()
        temp = fp.GetOnBits()
    except:
        size = None
        temp = list(fp.GetNonzeroElements().keys())
    res = [i for i in temp]

    return res, size
