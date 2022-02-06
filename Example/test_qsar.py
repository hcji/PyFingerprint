# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 10:43:57 2022

@author: jihon
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tqdm import tqdm
from rdkit import Chem
from PyFingerprint.fingerprint import get_fingerprint

import tensorflow.keras.backend as K
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Dense, Input
from tensorflow.keras import optimizers
from sklearn.model_selection import train_test_split
from lifelines.utils import concordance_index

class MLP:
    def __init__(self, X, Y):
        self.X = X
        self.Y = Y
        self.X_tr, self.X_ts, self.Y_tr, self.Y_ts = train_test_split(X, Y, test_size=0.1)
        
        inp = Input(shape=(X.shape[1],))
        hid = inp
        n = X.shape[1]
        for j in range(3):
            hid = Dense(n, activation="relu")(hid)
            n = int(n * 0.5)
        prd = Dense(1, activation="linear")(hid)
        opt = optimizers.Adam(lr=0.001)
        model = Model(inp, prd)
        model.compile(optimizer=opt, loss='mse', metrics=['mae'])
        self.model = model
        
    def train(self, epochs=10):
        self.model.fit(self.X_tr, self.Y_tr, epochs=epochs)

    def test(self):
        Y_pred = self.model.predict(self.X_ts)[:,0]
        corr = concordance_index(Y_pred, self.Y_ts)
        return corr
    
    def predict(self, X_pd):
        Y_pred = self.model.predict(X_pd)[:,0]
        return Y_pred
    

'''
supp = Chem.SDMolSupplier('d:/BindingDB_All_2D.sdf')

# take EGFR as example
data = []
for m in tqdm(supp):
    if m is None:
        continue
    if m.GetProp('UniProt (SwissProt) Primary ID of Target Chain') != 'P00533':
        continue
    else:
        smi = Chem.MolToSmiles(m)
        ic50 = m.GetProp('IC50 (nM)')
        ki = m.GetProp('Ki (nM)')
        kd = m.GetProp('Kd (nM)')
        data.append([smi, ic50, ki, kd])
data = pd.DataFrame(data)
data.columns = ['smiles', 'Ki', 'Kd', 'IC50']
data.to_csv('Example/qsar_egfr.csv', index = False)
'''

fp_types = ['standard', 'extended', 'maccs', 'pubchem', 'klekota-roth', 'morgan', 'mol2vec']

data = pd.read_csv('Example/qsar_egfr.csv')
uni_smi = np.unique(data['smiles'])

evaluation = {}
for fp_type in fp_types:
    X, Y = [], []
    for smi in tqdm(uni_smi):
        ys = data.loc[data['smiles'] == smi, 'Ki'].values
        yi = np.nanmean(pd.to_numeric(ys, errors='coerce'))
        yi = np.log10(yi)
        if np.isnan(yi):
            continue
        else:
            xi = get_fingerprint(smi, fp_type).to_numpy()
        X.append(xi)
        Y.append(yi)
    X = np.array(X)
    Y = np.array(Y)
    
    scores = []
    for i in range(20):
        model = MLP(X, Y)
        model.train()
        score = model.test()
        scores.append(score)
    evaluation[fp_type] = scores
    

import json
output = json.dumps(evaluation)
f = open("dict.json","w")
f.write(output)
f.close()