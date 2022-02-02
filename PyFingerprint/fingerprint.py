# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 08:01:30 2022

@author: jihon
"""


import numpy as np
from PyFingerprint.cdk import cdk_fingerprint
from PyFingerprint.rdk import rdk_fingerprint 
from PyFingerprint.babel import ob_fingerprint 

cdktypes = ['standard', 'extended', 'graph', 'maccs', 'pubchem', 'estate', 'hybridization', 'lingo', 
            'klekota-roth', 'shortestpath', 'signature', 'substructure']
rdktypes = ['rdkit', 'morgan', 'rdk-maccs', 'topological-torsion', 'avalon', 'atom-pair']
babeltypes = ['fp2', 'fp3', 'fp4']


class fingerprint:
    
    def __init__(self, bits: list, n: int):
        self.bits = bits
        self.n = n
        
    def __str__(self):
        return (self.bits, self.n)
    
    def check(self):
        for i in self.bits:
            if type(i) is int:
                if i < self.n:
                    pass
                else:
                    raise TypeError
            else:
                raise TypeError    
        
    def to_numpy(self):
        if self.n is None:
            return None
        elif self.n > 10024:
            return None
        v = np.zeros(self.n, dtype=np.uint8)
        for k in self.bits:
            v[k] = 1
        return v


def get_fingerprint(smi, fp_type, nbit=None, depth=None):
    if nbit is None:
        nbit = 1024
    if depth is None:
        depth = 6
    if fp_type in cdktypes:
        res, n = cdk_fingerprint(smi, fp_type, size = nbit, depth = depth)
    elif fp_type in rdktypes:
        res, n = rdk_fingerprint(smi, fp_type, size = nbit)
    elif fp_type in babeltypes:
        if nbit is None:
            nbit = 307
        res, n = ob_fingerprint(smi, fp_type)
    else:
        raise IOError('invalid fingerprint type')
    return fingerprint(res, n)

