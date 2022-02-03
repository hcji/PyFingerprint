# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 13:58:28 2022

@author: jihon
"""

import numpy as np
import pandas as pd

from faiss import IndexBinaryHNSW, IndexFlatL2
from sklearn.metrics import pairwise_distances
from PyFingerprint.fingerprint import get_fingerprint


def cal_similarity(query: list, smlist: list, metric='jaccard', topn=10, fp_type='standard', nbit=None, depth=None):
    qlist = [get_fingerprint(s, fp_type=fp_type, nbit=nbit, depth=depth).to_numpy() for s in query]
    rlist = [get_fingerprint(s, fp_type=fp_type, nbit=nbit, depth=depth).to_numpy() for s in smlist]
    qlist, rlist = np.array(qlist), np.array(rlist)
    dist_matrix = 1 - pairwise_distances(qlist, rlist, metric=metric)
    
    output = {}
    for i, smi in enumerate(query):
        k = np.argsort(-dist_matrix[i,:])[:topn]
        refs = smlist[k]
        score = dist_matrix[i, k]
        output[smi] = dict(zip(refs, score))
    return output
    


def cal_similarity_large(query: list, smlist: list, metric='hamming', topn=10, fp_type='standard', nbit=None, depth=None):
    qlist = [get_fingerprint(s, fp_type=fp_type, nbit=nbit, depth=depth).to_numpy() for s in query]
    rlist = [get_fingerprint(s, fp_type=fp_type, nbit=nbit, depth=depth).to_numpy() for s in smlist]

    if metric == 'hamming':
        qlist, rlist = np.array(qlist).astype('uint8'), np.array(rlist).astype('uint8')
        dim = rlist.shape[1] * 8
        index = IndexBinaryHNSW(dim)
        index.add(rlist)
        D, I = index.search(qlist, topn)
        dist_matrix = 1 - D / rlist.shape[1]
    elif metric == 'euclidean':
        qlist, rlist = np.array(qlist).astype('float32'), np.array(rlist).astype('float32')
        dim = rlist.shape[1]
        index = IndexFlatL2(dim)
        index.add(rlist)
        D, I = index.search(qlist, topn)
        dist_matrix = 1 - D

    output = {}
    for i, smi in enumerate(query):
        k = np.argsort(-dist_matrix[i,:])[:topn]
        refs = smlist[k]
        score = dist_matrix[i, k]
        output[smi] = dict(zip(refs, score))
    return output    