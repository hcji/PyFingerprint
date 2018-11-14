# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 12:46:18 2018

@author: hcji
"""

import numpy as np
from jpype import isJVMStarted, startJVM, getDefaultJVMPath, JPackage

if not isJVMStarted():
    cdk_path = 'CDK/cdk-2.2.jar'
    startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % cdk_path)
    cdk = JPackage('org').openscience.cdk

_fingerprinters = {"daylight":cdk.fingerprint.Fingerprinter
                            , "graph":cdk.fingerprint.GraphOnlyFingerprinter
                            , "maccs":cdk.fingerprint.MACCSFingerprinter
                            , "estate":cdk.fingerprint.EStateFingerprinter
                            , "extended":cdk.fingerprint.ExtendedFingerprinter
                            , "hybridization":cdk.fingerprint.HybridizationFingerprinter
                            , "klekota-roth":cdk.fingerprint.KlekotaRothFingerprinter
                            , "pubchem":cdk.fingerprint.PubchemFingerprinter
                            , "substructure":cdk.fingerprint.SubstructureFingerprinter
                            }

def cdk_parser_smiles(smi):
    sp = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())
    try:
        mol = sp.parseSmiles(smi)
    except:
        raise IOError('invalid smiles input')
    return mol


def cdk_fingerprint(smi, fp_type="daylight", nbit=1024):
    if fp_type == 'maccs':
        nbit = 166
    elif fp_type == 'estate':
        nbit = 79
    elif fp_type == 'pubchem':
        nbit = 881
    elif fp_type == 'kr':
        nbit = 4860
    
    mol = cdk_parser_smiles(smi)
    if fp_type in _fingerprinters:
        fingerprinter = _fingerprinters[fp_type]()
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
    return vec


if __name__ == '__main__':
    smi = 'CCCN'
    cdk_parser_smiles(smi)
    