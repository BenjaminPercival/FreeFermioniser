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
from FFerGetModelDetails import getNumBVs, BasDotProds, IsModInvG, IsModInvB,LCMs,readBasis,NumbSecs,\
GetGGSOMSecBas,GetAllSecUnRedBC,GetAllSecRedBC,GetAllSecs,MasslessUnRedSecs,MasslessSecs,MSecDeltas #GetGGSOMSecMSec
#import psyco  # pip install psychopy
#psyco.full()
import timeit
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
Nis=LCMs(Basis) #LCMs of basis vecs
#ModInvB=IsModInvB(NumBVs,BP,Basis,Nis)
#NSec=NumbSecs(NumBVs,Nis)
#AllUnRedScsBC=GetAllSecUnRedBC(NumBVs,Basis,NSec,Nis)
#AllRedScsBC=GetAllSecRedBC(AllUnRedScsBC)
#AllScs=GetAllSecs(NumBVs,NSec,Nis)
#MSecBCUnRed=MasslessUnRedSecs(AllUnRedScsBC)
#MSecs=MasslessSecs(AllRedScsBC,AllScs,NSec,Basis) #hstacked secs and secBCs
#MSecDelts=MSecDeltas(MSecs,NumBVs)
#MSecGSOs=GetGGSO(MSecs,MSecBCUnRed,NumBVs,MSecDelts,Basis)
t1=timeit.default_timer()
if __name__=='__main__':
    
    ModInvG=IsModInvG(NumBVs,BP,GSO) #GGSO phase matric MI Boolean 
    #print("BP is: ", BP)
    if ModInvG is True:
        #print("GSO mat is MI")
        InpBasisForm=readBasis(Basis)
        if InpBasisForm is True:
            #print("basis in good form")
            ModInvB=IsModInvB(NumBVs,BP,Basis,Nis)
            if ModInvB is True:
                print("Input is in good order")
                NSec=NumbSecs(NumBVs,Nis)
                #t2=timeit.default_timer()
                #print(t2-t1)
                AllUnRedScsBC=GetAllSecUnRedBC(NumBVs,Basis,NSec,Nis)
                #t3=timeit.default_timer()
                #print(t3-t2)
                AllRedScsBC=GetAllSecRedBC(AllUnRedScsBC)
                #t4=timeit.default_timer()
                #print(t4-t3)
                AllScs=GetAllSecs(NumBVs,NSec,Nis)
                #t5=timeit.default_timer()
                #print(t5-t4)
                MSecBCUnRed=MasslessUnRedSecs(AllUnRedScsBC)
                #t6=timeit.default_timer()
                #print(t6-t5)
                MSecs=MasslessSecs(AllRedScsBC,AllScs,NSec,Basis) #hstacked secs and secBCs
                #t7=timeit.default_timer()
                #print(t7-t6)
                #print(MSecs)
                MSecDelts=MSecDeltas(MSecs,NumBVs)
                MSecGSOsMB=GetGGSOMSecBas(MSecs,MSecBCUnRed,NumBVs,MSecDelts,Basis,GSO)
                #MSecGSOsMM=GetGGSOMSecMSec(MSecs,MSecBCUnRed,NumBVs,MSecDelts,Basis,GSO)
                #print(MSecGSOsMM)
                #FreeFermioniser create output file
                
                #print(readBasis()) # read input basis and output as in form 1,S,e1, etc + check in appropriate form- 44 entries
                
                #print(IsModInv()) #if not then bottom out
                #print(supercurrent()) #if not then bottom out
                #print(IsSymmetric()) #if not then bottom out
                
                #print(NSSec(NumBVs))
                #print(SSec())#in this function account for possible Stilde
                #print("Model's spacetime SUSY is N=", ModSusy())
                #print(enhancements())
                #print(ObservableMassless())
                #print(OnShellTachyons())
                #print(Exotics())
                #print(Hidden())
