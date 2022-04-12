#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 17:08:20 2022

@author: wmnc67
"""
#from FFerGetModelDetails.py import *
import sys 
from z3 import * #pip install z3-solver
import numpy as np
import sys
import itertools
from SpectrumFunctions import NSSec, SSec
from FFerGetModelDetails import getNumBVs, BasDotProds, IsModInvG, IsModInvB
import psyco  # pip install psychopy
psyco.full()

#################################################################################
#INPUT DATA
InputBasisFile = "Input/InBasis.txt"  # File Containing the Basis Vectors
InputGSOFile = "Input/InGSO.txt"  # File Containing the GSO Coefficients
#OutputFile = "OutputModel.txt"  # Destination for the q-Exp. Data
##################################################################################


#load input data
with open(InputBasisFile, "r") as InBasis:
     Basis = np.loadtxt(InBasis)
with open(InputGSOFile, "r") as InGSO:
     GSO = np.loadtxt(InGSO, dtype=complex)

#extract variables from input
NumBVs=getNumBVs(Basis)  # num basis vecs
BP=BasDotProds(Basis,GSO) # basis dot prods as numpy array
ModInvG=IsModInvG(NumBVs,BP,GSO) #GGSO phase matric MI Boolean 
Nis=LCMs(Basis) #LCMs of basis vecs
ModInvB=IsModInvB(NumBVs,BP,Basis,Nis)
NSec=NumbSecs(NumBVs,Nis)

#FreeFermioniser create output file

#print(readBasis()) # read input basis and output as in form 1,S,e1, etc + check in appropriate form- 44 entries

#print(IsModInv()) #if not then bottom out
#print(supercurrent()) #if not then bottom out
#print(IsSymmetric()) #if not then bottom out

print(NSSec(NumBVs))
print(SSec())#in this function account for possible Stilde
#print("Model's spacetime SUSY is N=", ModSusy())
#print(enhancements())
#print(ObservableMassless())
#print(OnShellTachyons())
#print(Exotics())
#print(Hidden())
