# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 11:42:12 2022

@author: jihon
"""


import numpy as np
import pandas as pd

from PyFingerprint.fingerprint import get_fingerprint

data = pd.read_csv('Example/DrugBank_DTI.csv')
smlist = np.unique([s for s in data['SMILES'].values if s is not np.nan])

fp_type = 'standard'

