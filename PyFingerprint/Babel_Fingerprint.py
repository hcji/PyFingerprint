# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:49:26 2018

@author: hcji
"""

import pybel
import numpy as np

def ob_fingerprint(smi, fp_type='FP2', nbit=307):
    mol = pybel.readstring("smi", smi)
    if fp_type == 'FP2':
        fp = mol.calcfp('FP2')
    elif fp_type == 'FP3':
        fp = mol.calcfp('FP3')
    elif fp_type == 'FP4':
        fp = mol.calcfp('FP4')
    bits = fp.bits
    bits = [x for x in bits if x < nbit]
    vec = np.zeros(nbit)
    vec[bits] = 1
    vec = vec.astype(int)
    return vec

