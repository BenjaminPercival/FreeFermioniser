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
from SpectrumFunctions import UnprojectedSecs, NSSec, Massless40, printUnProjdSecs, Massless48, Massless44, Massless04, Massless08
from FFerGetModelDetails import getNumBVs, BasDotProds, IsModInvG, IsModInvB,LCMs,readBasis,NumbSecs,supercurrent,\
GetGGSOMSecBas,GetAllSecUnRedBC,GetAllSecRedBC,GetAllSecs,MasslessUnRedSecs,MasslessSecs,MSecDeltas,MSectVacEs,GetGGSOMSecMSec
#import psyco  # pip install psychopy
#psyco.full()
import timeit
import os
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

#delete previous output file (if there)
if os.path.exists("OutputModel.txt"):
    os.remove("OutputModel.txt")
    print("Old output file has been deleted successfully")
else:
    pass

#extract variables from input
NumBVs=getNumBVs(Basis)  # num basis vecs
BP=BasDotProds(Basis,GSO) # basis dot prods as numpy array
Nis=LCMs(Basis) #LCMs of basis vecs

t1=timeit.default_timer()

if __name__=='__main__':
    
    ModInvG=IsModInvG(NumBVs,BP,GSO) 
    #print("BP is: ", BP)
    if ModInvG is True: #GGSO phase matric MI Boolean 
        #print("GSO mat is MI")
        InpBasisForm=readBasis(Basis)
        if InpBasisForm is True: #basis in good form- Symmetric, right dimensions..
            #print("basis in good form")
            SCurrent=supercurrent(NumBVs,Basis)
            if SCurrent is True: #supercurrent constraint satisfied
                ModInvB=IsModInvB(NumBVs,BP,Basis,Nis)
                if ModInvB is True: #basis satisfies modular invariance rules
                    print("Input is all consistent and as desired")
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
                    MSecGSOsMM=GetGGSOMSecMSec(MSecs,MSecBCUnRed,NumBVs,MSecDelts,GSO)
                    #print(MSecGSOsMM)
                    MScVEs=MSectVacEs(MSecs,NumBVs)
                    #print(MScVEs)
                    UnProjdSecs=UnprojectedSecs(NumBVs,MSecs,MSecGSOsMM,MSecDelts,Basis,MScVEs)
                    printUnProjdSecs(MSecs,UnProjdSecs,NumBVs)
                    #NSSec(NumBVs,Basis)
                    Massless40(Basis,NumBVs,MSecGSOsMB,UnProjdSecs)
                    Massless48(Basis,NumBVs,MSecGSOsMB,UnProjdSecs)
                    Massless44(Basis,NumBVs,MSecGSOsMB,UnProjdSecs)
                    Massless04(Basis,NumBVs,MSecGSOsMB,UnProjdSecs)
                    Massless08(Basis,NumBVs,MSecGSOsMB,UnProjdSecs)
