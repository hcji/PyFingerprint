# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:49:26 2018

@author: hcji
"""

from openbabel import pybel

def ob_fingerprint(smi, fp_type='FP2', nbit=1024):
    

    mol = pybel.readstring("smi", smi)
    if fp_type == 'fp2':
        fp = mol.calcfp('FP2')
    elif fp_type == 'fp3':
        fp = mol.calcfp('FP3')
    elif fp_type == 'fp4':
        fp = mol.calcfp('FP4')
    bits = fp.bits
    bits = [x for x in bits if x < nbit]
    return bits, nbit
