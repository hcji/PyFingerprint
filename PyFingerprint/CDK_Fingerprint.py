# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 12:46:18 2018

@author: hcji
"""


import numpy as np
from jpype import isJVMStarted, startJVM, getDefaultJVMPath, JPackage
import PyFingerprint

if not isJVMStarted():
    cdk_path = PyFingerprint.__path__[0] + '\\CDK\\cdk-2.2.jar'
    startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % cdk_path)
    cdk = JPackage('org').openscience.cdk
    

def cdk_parser_smiles(smi):
    sp = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())
    try:
        mol = sp.parseSmiles(smi)
    except:
        raise IOError('invalid smiles input')
    return mol


def cdk_fingerprint(smi, fp_type="daylight", size=1024, depth=6):
    if fp_type == 'maccs':
        nbit = 166
    elif fp_type == 'estate':
        nbit = 79
    elif fp_type == 'pubchem':
        nbit = 881
    elif fp_type == 'klekota-roth':
        nbit = 4860
    else:
        nbit = size
        
    _fingerprinters = {"daylight":cdk.fingerprint.Fingerprinter(size, depth)
                            , "extended":cdk.fingerprint.ExtendedFingerprinter(size, depth)
                            , "graph":cdk.fingerprint.GraphOnlyFingerprinter(size, depth)
                            , "maccs":cdk.fingerprint.MACCSFingerprinter()
                            , "pubchem":cdk.fingerprint.PubchemFingerprinter(cdk.silent.SilentChemObjectBuilder.getInstance())
                            , "estate":cdk.fingerprint.EStateFingerprinter()
                            , "hybridization":cdk.fingerprint.HybridizationFingerprinter(size, depth)
                            , "lingo":cdk.fingerprint.LingoFingerprinter(depth)
                            , "klekota-roth":cdk.fingerprint.KlekotaRothFingerprinter()
                            , "shortestpath":cdk.fingerprint.ShortestPathFingerprinter(size)
                            , "signature": cdk.fingerprint.SignatureFingerprinter(depth)
                            , "circular": cdk.fingerprint.CircularFingerprinter()
                            }
    
    mol = cdk_parser_smiles(smi)
    if fp_type in _fingerprinters:
        fingerprinter = _fingerprinters[fp_type]
    else:
        raise IOError('invalid fingerprint type')
        
    fp = fingerprinter.getBitFingerprint(mol).asBitSet()
    bits = []
    idx = fp.nextSetBit(0)
    while idx >= 0:
        bits.append(idx)
        idx = fp.nextSetBit(idx + 1)
    vec = np.zeros(nbit)
    vec[bits] = 1
    vec = vec.astype(int)
    
    return vec

    