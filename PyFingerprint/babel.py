# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 10:49:26 2018

@author: hcji

Updated on Thu Oct 12 20:04:22 2023

@author: Jnelen
"""

from openbabel import openbabel, pybel

import numpy as np

def convertToFP(bits, nbit):
    fp = np.zeros(nbit)
    for i in bits:
        fp[i] = 1
    return list(fp)
    
def ob_fingerprint(smi, fp_type='FP2', depth=0):

    mol = pybel.readstring("smi", smi)
    
    if fp_type == 'fp2':
        bits = mol.calcfp('FP2').bits
        nbit = 1024
        fp = convertToFP(bits, nbit)
    elif fp_type == 'fp3':
        bits = mol.calcfp('FP3').bits
        nbit = 55
        fp = convertToFP(bits, nbit)
    elif fp_type == 'fp4':
        bits = mol.calcfp('FP4').bits
        nbit = 307
        fp = convertToFP(bits, nbit)
    elif fp_type == 'spectrophore':
        mol.addh()
        mol.make3D()       
        spectrophoreCalculator = openbabel.OBSpectrophore()
        spectrophoreCalculator.SetNormalization(depth)
        fp = list(spectrophoreCalculator.GetSpectrophore(mol.OBMol))
    else:
        raise IOError('invalid fingerprint type')       
    return fp
