# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 10:11:11 2018

@author: hcji
"""

import numpy as np
from PyFingerprint.CDK_Fingerprint import cdk_fingerprint
from PyFingerprint.RDK_Fingerprint import rdk_fingerprint
try:
    from PyFingerprint.Babel_Fingerprint import ob_fingerprint
except:
    print('Cannot import pybal, openbabel fingerprint cannot be used')


def get_fingerprint(smi, fp_type, nbit=None, depth=None, output='bit'):
    if fp_type in ["daylight", "extended", "graph", "pubchem", "estate", "hybridization", "lingo", "klekota-roth", "shortestpath", "signature", "circular", "Morgan"]:
        if nbit is None:
            nbit = 1024
        if depth is None:
            depth = 6
        res = cdk_fingerprint(smi, fp_type, nbit, depth, output)
    elif fp_type in ["rdkit", "maccs", "AtomPair", "TopologicalTorsion", "Avalon"]:
        res = rdk_fingerprint(smi, fp_type, nbit, output)
    elif fp_type in ["FP2", "FP3", "FP4"]:
        if nbit is None:
            nbit = 307
        res = ob_fingerprint(smi, fp_type, nbit, output)
    else:
        raise IOError('invalid fingerprint type')
    return res

def get_multi_fingerprint(smi, fp_types, nbit=None, depth=None, output='bit'):
    vecs = [get_fingerprint(smi, f, nbit, depth, output='vector') for f in fp_types]
    vec = np.concatenate(vecs)
    if output == 'bit':
        bits = np.where(vec > 0)[0]
        return bits
    else:
        return vec

    