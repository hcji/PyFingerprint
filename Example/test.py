# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 09:51:17 2022

@author: jihon
"""

from PyFingerprint.fingerprint import get_fingerprint

smi = 'CCCCN'
cdktypes = ['standard', 'extended', 'graph', 'maccs', 'pubchem', 'estate', 'hybridization', 'lingo', 
            'klekota-roth', 'shortestpath', 'signature', 'substructure']
rdktypes = ['rdkit', 'morgan', 'rdk-maccs', 'topological-torsion', 'avalon', 'atom-pair']
babeltypes = ['fp2', 'fp3', 'fp4']

output = {}
for f in cdktypes:
    output[f] = get_fingerprint(smi, f)

for f in rdktypes:
    output[f] = get_fingerprint(smi, f)
    
for f in babeltypes:
    output[f] = get_fingerprint(smi, f)

output_np = output.copy()
for k, fp in output.items():
    output_np[k] = fp.to_numpy()