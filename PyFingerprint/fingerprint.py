# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 08:01:30 2022

@author: jihon
"""


import numpy as np
from PyFingerprint.cdk import cdk_fingerprint
from PyFingerprint.rdk import rdk_fingerprint 
from PyFingerprint.babel import ob_fingerprint
from PyFingerprint.mol2vec import mol2vec_fingerprint, mol2vec_fingerprints
from PyFingerprint.heteroencoder import hc_fingerprint, hc_fingerprints

cdktypes = ['standard', 'extended', 'graph', 'maccs', 'pubchem', 'estate', 'hybridization', 'lingo', 
            'klekota-roth', 'shortestpath', 'signature', 'substructure']
rdktypes = ['rdkit', 'morgan', 'rdk-maccs', 'topological-torsion', 'avalon', 'atom-pair']
babeltypes = ['fp2', 'fp3', 'fp4']
vectypes = ['mol2vec', 'heteroencoder']


class fingerprint:
    
    def __init__(self, bits: list, values: list, n: int):
        self.bits = bits
        self.values = values
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
        v = np.zeros(self.n)
        for i, k in enumerate(self.bits):
            v[k] = self.values[i]
        return v


def get_fingerprint(smi: str, fp_type: str, nbit=None, depth=None):
    if nbit is None:
        nbit = 1024
    if depth is None:
        depth = 6
    if fp_type in cdktypes:
        bits, n = cdk_fingerprint(smi, fp_type, size = nbit, depth = depth)
        values = [1] * len(bits)
    elif fp_type in rdktypes:
        bits, n = rdk_fingerprint(smi, fp_type, size = nbit)
        values = [1] * len(bits)
    elif fp_type in babeltypes:
        if nbit is None:
            nbit = 307
        bits, n = ob_fingerprint(smi, fp_type)
        values = [1] * len(bits)
    elif fp_type == 'mol2vec':
        values = list(mol2vec_fingerprint(smi))
        n = len(values)
        bits = list(np.arange(n))
    elif fp_type == 'heteroencoder':
        values = list(hc_fingerprint(smi))
        n = len(values)
        bits = list(np.arange(n))
    else:
        raise IOError('invalid fingerprint type')
    return fingerprint(bits, values, n)


def get_fingerprints(smlist: list, fp_type: str, nbit=None, depth=None):
    if fp_type not in vectypes:
        output = [get_fingerprint(smi, fp_type, nbit, depth) for smi in smlist]
    elif fp_type == 'mol2vec':
        vecs = mol2vec_fingerprints(smlist)
        n = vecs.shape[1]
        bits = list(np.arange(n))
        output = [fingerprint(bits, vecs[i,:], n) for i in range(vecs.shape[0])]
    elif fp_type == 'heteroencoder':
        vecs = hc_fingerprints(smlist)
        n = vecs.shape[1]
        bits = list(np.arange(n))
        output = [fingerprint(bits, vecs[i,:], n) for i in range(vecs.shape[0])]        
    return output
