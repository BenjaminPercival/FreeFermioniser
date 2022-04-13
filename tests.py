#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 15:22:35 2022

@author: viktor
"""
import itertools
import numpy as np

InputBasisFile = "Input/InBasis.txt"
with open(InputBasisFile, "r") as InBasis:
     Basis = np.loadtxt(InBasis)

lcms=[2,2]
NSects=4
NBV=2
AllSectorBC = np.zeros((NSects,Basis.shape[1]))
AllSector = np.zeros((NSects,NBV))
for i,t in enumerate(itertools.product(*[range(j) for j in lcms])):
    AllSectorBC[i,:] = sum([Basis[i,:] * t[i] for i in range(len(t))])
    AllSector[i,:] = t

print(AllSectorBC)