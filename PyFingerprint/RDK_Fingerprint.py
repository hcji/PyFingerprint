# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 15:39:39 2018

@author: hcji
"""

import numpy as np
from rdkit.Chem import MACCSkeys
from rdkit import Chem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Avalon import pyAvalonTools

def rdk_fingerprint(smi, fp_type="rdkit", size=1024, output="bit"):
    _fingerprinters = {"rdkit": Chem.rdmolops.RDKFingerprint
                           , "maccs": MACCSkeys.GenMACCSKeys
                           , "TopologicalTorsion": Torsions.GetTopologicalTorsionFingerprint
                           , "Avalon": pyAvalonTools.GetAvalonFP}
    mol = Chem.MolFromSmiles(smi)
    if fp_type in _fingerprinters:
        fingerprinter = _fingerprinters[fp_type]
        fp = fingerprinter(mol)
    elif fp_type == "AtomPair":
        fp = Pairs.GetAtomPairFingerprintAsBitVect(mol, nBits=size)
    elif fp_type == "Morgan":
        fp = GetMorganFingerprintAsBitVect(mol, 2, nBits=size)
    else:
        raise IOError('invalid fingerprint type')
    if output == "bit":
        temp = fp.GetOnBits()
        res = [i for i in temp]
    else:
        res = np.array(fp)
    return res
    