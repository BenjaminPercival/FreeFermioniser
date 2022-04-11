#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 17:08:20 2022

@author: wmnc67
"""
#from FFerGetModelDetails.py import *
import SpectrumFunctions
#FreeFermioniser create output file

#print(readBasis()) # read input basis and output as in form 1,S,e1, etc + check in appropriate form- 44 entries

#print(IsModInv()) #if not then bottom out
#print(supercurrent()) #if not then bottom out
#print(IsSymmetric()) #if not then bottom out

print(SpectrumFunctions.NSSec())
print(SpectrumFunctions.SSec())#in this function account for possible Stilde
#print("Model's spacetime SUSY is N=", ModSusy())
#print(enhancements())
#print(ObservableMassless())
#print(OnShellTachyons())
#print(Exotics())
#print(Hidden())
