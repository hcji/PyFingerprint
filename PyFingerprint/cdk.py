# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 12:46:18 2018

@author: hcji

Updated on Thu Oct 12 20:04:22 2023

@author: Jnelen
"""

import os
import numpy as np
from jpype import isJVMStarted, startJVM, getDefaultJVMPath, JPackage
import PyFingerprint

if not isJVMStarted():
    cdk_path = os.path.join(PyFingerprint.__path__[0], 'CDK', 'cdk-2.9.jar')
    startJVM(getDefaultJVMPath(), "-ea", "-Djava.class.path=%s" % cdk_path)
    cdk = JPackage('org').openscience.cdk
    

sp = cdk.smiles.SmilesParser(cdk.DefaultChemObjectBuilder.getInstance())
def cdk_parser_smiles(smi):
    try:
        mol = sp.parseSmiles(smi)
    except:
        raise IOError('invalid smiles input')
    return mol

fp_map  = {}

def get_fingerprinter(name, size, depth): 
    ### This was getting made every time!
    ## Checking if the depth is specified. If not, we use the default CDK values where appropriate 
    if depth == None:
	    _fingerprinters = {"standard":lambda : cdk.fingerprint.Fingerprinter(size, 7)
	                       , "cdk-atompairs": lambda : cdk.fingerprint.AtomPairs2DFingerprinter()
	                       , "extended":lambda : cdk.fingerprint.ExtendedFingerprinter(size, 7)
	                       , "graph":lambda : cdk.fingerprint.GraphOnlyFingerprinter(size, 7)
	                       , "maccs":lambda : cdk.fingerprint.MACCSFingerprinter()
	                       , "pubchem":lambda : cdk.fingerprint.PubchemFingerprinter(cdk.silent.SilentChemObjectBuilder.getInstance())
	                       , "estate":lambda : cdk.fingerprint.EStateFingerprinter()
	                       , "hybridization":lambda : cdk.fingerprint.HybridizationFingerprinter(size, 7)
	                       , "lingo":lambda : cdk.fingerprint.LingoFingerprinter()
	                       , "klekota-roth":lambda : cdk.fingerprint.KlekotaRothFingerprinter()
	                       , "shortestpath":lambda : cdk.fingerprint.ShortestPathFingerprinter(size)
	                       , "signature": lambda : cdk.fingerprint.SignatureFingerprinter()
	                       ## circular fingerprint defaults to ECFP6: https://github.com/cdk/cdk/blob/125505c5ea1f69b692183bb0aae65816e7cb44e7/descriptor/fingerprint/src/main/java/org/openscience/cdk/fingerprint/CircularFingerprinter.java
	                       , "circular": lambda : cdk.fingerprint.CircularFingerprinter(4,size)
	                       , "cdk-substructure": lambda : cdk.fingerprint.SubstructureFingerprinter()
	                            }
    ## Use the user-specified settings for the fingerprint generation
    else:
        _fingerprinters = {"standard":lambda : cdk.fingerprint.Fingerprinter(size, depth)
	                       , "cdk-atompairs": lambda : cdk.fingerprint.AtomPairs2DFingerprinter()
	                       , "extended":lambda : cdk.fingerprint.ExtendedFingerprinter(size, depth)
	                       , "graph":lambda : cdk.fingerprint.GraphOnlyFingerprinter(size, depth)
	                       , "maccs":lambda : cdk.fingerprint.MACCSFingerprinter()
	                       , "pubchem":lambda : cdk.fingerprint.PubchemFingerprinter(cdk.silent.SilentChemObjectBuilder.getInstance())
	                       , "estate":lambda : cdk.fingerprint.EStateFingerprinter()
	                       , "hybridization":lambda : cdk.fingerprint.HybridizationFingerprinter(size, depth)
	                       , "lingo":lambda : cdk.fingerprint.LingoFingerprinter(depth)
	                       , "klekota-roth":lambda : cdk.fingerprint.KlekotaRothFingerprinter()
	                       , "shortestpath":lambda : cdk.fingerprint.ShortestPathFingerprinter(size)
	                       , "signature": lambda : cdk.fingerprint.SignatureFingerprinter(depth)
	                       , "circular": lambda : cdk.fingerprint.CircularFingerprinter(depth, size)
	                       , "cdk-substructure": lambda : cdk.fingerprint.SubstructureFingerprinter()
	                            }
	                            
    if name not in _fingerprinters: 
        raise IOError('invalid fingerprint type')

    return _fingerprinters[name]()

def cdk_fingerprint(smi, fp_type="standard", size=1024, depth=None):

    mol = cdk_parser_smiles(smi)
    ## Sanitize input molecules, as is recommended for most fingerprints (especially shortestpath)
    cdk.tools.manipulator.AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(mol)
    cdk.tools.manipulator.AtomContainerManipulator.convertImplicitToExplicitHydrogens(mol)
    
    if fp_type == 'estate':
        nbit = 79
    elif fp_type == 'maccs':
        nbit = 166
    elif fp_type == 'cdk-substructure':
        nbit = 307
    elif fp_type == 'cdk-atompairs':
        nbit = 780
    elif fp_type == 'pubchem':
        nbit = 881
    elif fp_type == 'klekota-roth':
        nbit = 4860      
    elif fp_type == 'signature':
        nbit = None
        print("Signature_FP")
        fingerprinter = cdk.fingerprint.SignatureFingerprinter()
        mol = cdk_parser_smiles(smi)
        print(fingerprinter.getSize())
        print(fingerprinter.getBitFingerprint(mol).getSetbits())
        print(fingerprinter.getBitFingerprint(mol).size())
        print(fingerprinter.getRawFingerprint(mol))
        
    else:
        nbit = size

    

    # Pull from cache if it exists
    if (fp_type, size, depth) in fp_map: 
        fingerprinter = fp_map[(fp_type, size, depth)]
    else:
        fingerprinter = get_fingerprinter(fp_type, size, depth)
        fp_map[(fp_type, size, depth)] = fingerprinter
    
    fp_obj = fingerprinter.getBitFingerprint(mol)
    bits = list(fp_obj.getSetbits())
    return bits, nbit
