# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 08:01:30 2022

@author: jihon

Updated on Thu Oct 12 20:04:22 2023

@author: Jnelen
"""


import numpy as np
from PyFingerprint.cdk import cdk_fingerprint
from PyFingerprint.rdk import rdk_fingerprint 
from PyFingerprint.babel import ob_fingerprint
from PyFingerprint.mol2vec import mol2vec_fingerprint, mol2vec_fingerprints
try:
    from PyFingerprint.heteroencoder import hc_fingerprint, hc_fingerprints
    hc_enable = True
except:
    hc_enable = False
    

cdktypes = ['standard', 'extended', 'graph', 'maccs', 'pubchem', 'estate', 'hybridization', 'lingo', 'klekota-roth', 'shortestpath', 'cdk-substructure', 'circular', 'cdk-atompairs']
rdktypes = ['rdkit', 'morgan', 'rdk-maccs', 'topological-torsion', 'avalon', 'atom-pair', 'rdk-descriptor']
babeltypes = ['fp2', 'fp3', 'fp4', 'spectrophore']
vectypes = ['mol2vec', 'heteroencoder']

class fingerprint:
    
    def __init__(self, fp: list):
        self.fp = np.array(fp)
        self.bits = np.nonzero(self.fp)[0]
        self.values = self.fp[self.bits]
        self.n = len(self.fp)
        
    def __str__(self):
        return str(list(self.fp))[1:-1]
    
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
        return self.fp

    def getActiveBits(self):
        return self.bits
        
    def numActiveBits(self):
        return len(self.getActiveBits()) 
        
    def getActiveValues(self):
        return dict(zip(self.getActiveBits(), self.values)) 
           
    def to_str(self):          
        return str(list(self.fp))[1:-1]


def get_fingerprint(smi: str, fp_type: str, nbit=None, depth=None):
    if fp_type in cdktypes:
        if nbit is None:
            nbit = 1024
        bits, n = cdk_fingerprint(smi, fp_type, size=nbit, depth=depth)
        fp = np.zeros(n)
        for i, k in enumerate(bits):
            fp[k] = 1
    elif fp_type in rdktypes:
        fp = list(rdk_fingerprint(smi, fp_type, size=nbit))

    elif fp_type in babeltypes:
        fp = ob_fingerprint(smi, fp_type)
    elif fp_type == 'mol2vec':
        fp = list(mol2vec_fingerprint(smi))
    elif fp_type == 'heteroencoder':
        if not hc_enable:
            raise IOError('heteroencoder is not enabled')
        fp = list(hc_fingerprint(smi))
    else:
        raise IOError('invalid fingerprint type')
    return fingerprint(fp)


def get_fingerprints(smlist: list, fp_type: str, nbit=None, depth=None):
    if fp_type not in vectypes:
        output = [get_fingerprint(smi, fp_type, nbit, depth) for smi in smlist]
    elif fp_type == 'mol2vec':
        output = [get_fingerprint(smi, fp_type, nbit, depth) for smi in smlist]
    elif fp_type == 'heteroencoder':
        if not hc_enable:
            raise IOError('heteroencoder is not enabled')
        vecs = hc_fingerprints(smlist)
        n = vecs.shape[1]
        bits = list(np.arange(n))
        output = [fingerprint(bits, vecs[i,:], n) for i in range(vecs.shape[0])]        
    return output
