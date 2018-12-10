# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 10:11:11 2018

@author: hcji
"""

from PyFingerprint.CDK_Fingerprint import cdk_fingerprint
from PyFingerprint.RDK_Fingerprint import rdk_fingerprint
from PyFingerprint.Babel_Fingerprint import ob_fingerprint


def get_fingerprint(smi, fp_type, nbit=None, depth=None):
    if fp_type in ["daylight", "extended", "graph", "pubchem", "estate", "hybridization", "lingo", "klekota-roth", "shortestpath", "signature", "circular"]:
        if nbit is None:
            nbit = 1024
        if depth is None:
            depth = 6
        res = cdk_fingerprint(smi, fp_type, nbit, depth)
    elif fp_type in ["rdkit", "maccs", "AtomPair", "TopologicalTorsion", "Avalon"]:
        res = rdk_fingerprint(smi, fp_type)
    elif fp_type in ["FP2", "FP3", "FP4"]:
        if nbit is None:
            nbit = 307
        res = ob_fingerprint(smi, fp_type, nbit)
    else:
        raise IOError('invalid fingerprint type')
    return res
    