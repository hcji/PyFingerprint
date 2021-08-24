# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:49:26 2018

@author: hcji
"""

from openbabel import pybel
import numpy as np

def ob_fingerprint(smi, fp_type='FP2', nbit=307, output='bit'):
    mol = pybel.readstring("smi", smi)
    if fp_type == 'FP2':
        fp = mol.calcfp('FP2')
    elif fp_type == 'FP3':
        fp = mol.calcfp('FP3')
    elif fp_type == 'FP4':
        fp = mol.calcfp('FP4')
    bits = fp.bits
    bits = [x for x in bits if x < nbit]
    if output == 'bit':
        return bits
    else:
        vec = np.zeros(nbit, dtype=np.uint8)
        vec[bits] = 1
        return vec

def ob_fingerprints(smis, fp_type='FP2', nbit=307, output='bit'):
    output_list = []
    for smi in smis: 
        mol = pybel.readstring("smi", smi)
        if fp_type == 'FP2':
            fp = mol.calcfp('FP2')
        elif fp_type == 'FP3':
            fp = mol.calcfp('FP3')
        elif fp_type == 'FP4':
            fp = mol.calcfp('FP4')
        bits = fp.bits
        bits = [x for x in bits if x < nbit]
        if output == 'bit':
            output_list.append(bits)

        else:
            vec = np.zeros(nbit, dtype=np.uint8)
            vec[bits] = 1
            output_list.append(vec)
    return output_list

