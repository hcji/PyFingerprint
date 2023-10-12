# -*- coding: utf-8 -*-
"""
Created on Wed Nov 14 15:39:39 2018

@author: hcji

Updated on Thu Oct 12 20:04:22 2023

@author: Jnelen
"""

from rdkit.Chem import MACCSkeys
from rdkit import Chem
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Avalon import pyAvalonTools
from rdkit.Chem import Descriptors

def calcDescriptors(mol):
   fp = []
   for nm,fn in Descriptors.descList:
       try:
           fp.append(fn(mol))
       except:
           fp.append(None)
   return fp
    
def rdk_fingerprint(smi, fp_type="rdkit", depth=None, size=None):

    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)

    if not depth == None:
	    if fp_type == "avalon":	    
	        size = 512 if size == None else size
	        fp = pyAvalonTools.GetAvalonFP(mol, nBits=size)
	    elif fp_type == "atom-pair":
	        size = 2048 if size == None else size
	        fp = Pairs.GetHashedAtomPairFingerprint(mol, maxLength=depth, nBits=size)
	    elif fp_type == "morgan":
	        size = 2048 if size == None else size
	        fp = GetMorganFingerprintAsBitVect(mol, radius=depth, nBits=size)
	    elif fp_type == "rdkit":
	        size = 2048 if size == None else size
	        fp = Chem.rdmolops.RDKFingerprint(mol, maxPath=depth, fpSize=size)
	    elif fp_type == "rdk-maccs":
	        fp = MACCSkeys.GenMACCSKeys(mol)  
	    elif fp_type == "topological-torsion":
	        size = 2048 if size == None else size
	        fp = Torsions.GetHashedTopologicalTorsionFingerprint(mol, nBits=size, targetSize=size)
	    elif fp_type == "rdk-descriptor":
	        fp = calcDescriptors(mol)          
	    else:
	        raise IOError('invalid fingerprint type')
	        
    else:
	    if fp_type == "avalon":	    
	        size = 512 if size == None else size
	        fp = pyAvalonTools.GetAvalonFP(mol, nBits=size)
	    elif fp_type == "atom-pair":
	        size = 2048 if size == None else size
	        fp = Pairs.GetHashedAtomPairFingerprint(mol, nBits=size)
	    elif fp_type == "morgan":
	        size = 2048 if size == None else size
	        fp = GetMorganFingerprintAsBitVect(mol, radius=3, nBits=size)
	    elif fp_type == "rdkit":
	        size = 2048 if size == None else size
	        fp = Chem.rdmolops.RDKFingerprint(mol, fpSize=size)
	    elif fp_type == "rdk-maccs":
	        fp = MACCSkeys.GenMACCSKeys(mol)  
	    elif fp_type == "topological-torsion":
	        size = 2048 if size == None else size
	        fp = Torsions.GetHashedTopologicalTorsionFingerprint(mol, nBits=size)
	    elif fp_type == "rdk-descriptor":
	        fp = calcDescriptors(mol)
	    else:
	        raise IOError('invalid fingerprint type')

    return fp
