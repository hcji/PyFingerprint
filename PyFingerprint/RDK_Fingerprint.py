# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 15:39:39 2018

@author: hcji
"""

import numpy as np
from rdkit.Chem import MACCSkeys
from rdkit import Chem
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Avalon import pyAvalonTools

def rdk_fingerprint(smi, fp_type="rdkit"):
    _fingerprinters = {"rdkit": Chem.RDKFingerprint
                           , "maccs": MACCSkeys.GenMACCSKeys
                           , "AtomPair": Pairs.GetAtomPairFingerprint
                           , "TopologicalTorsion": Torsions.GetTopologicalTorsionFingerprint
                           , "Avalon": pyAvalonTools.GetAvalonFP
                            }
    mol = Chem.MolFromSmiles(smi)
    
    if fp_type in _fingerprinters:
        fingerprinter = _fingerprinters[fp_type]
        fp = fingerprinter(mol)
    else:
        raise IOError('invalid fingerprint type')    
    vec = np.array(fp)
    
    return vec
    