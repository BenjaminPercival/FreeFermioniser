#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 16:01:41 2022

@author: wmnc67
"""

import numpy as np
import sys
#from z3 import *
import itertools
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

#define an object:  model with the characteristics: number basis vecs, secs, dot prods,..
NumBVs = Basis.shape[0]

def DProd(vi,vj):
   DP = 0.5*np.dot(vi[0:20],vj[0:20]) - 0.5*np.dot(vi[20:32],vj[20:32]) \
          - np.dot(vi[32:44],vj[32:44])
   return DP

# Dot product of basis vectors
BP = np.zeros((len(GSO[0]),len(GSO[0])))
for i in range(len(GSO[0])):
    for k in range(len(GSO[0])):
        BP[i][k] = DProd(Basis[i],Basis[k])
  

# Check if basis GSOs are modular invariant
n1 = n2 = 0
for i in range(NumBVs):
    for k in range(NumBVs):
        TGSO1 = -np.around(np.exp(+1j*np.pi*BP[i][k]/4)*GSO[i][0])
        TGSO2 = np.around(np.exp(1j*np.pi*BP[i][k]/2)*np.conj(GSO[k][i]))
        if i == k and TGSO1 != GSO[i][k]:
            n1=n1+1
        elif i != k and TGSO2 != GSO[i][k]:
            n2=n2+1
if n1 != 0 or n2 != 0:
    print("~ Error: Basis GSOs not modular invariant!")
    #sys.exit()
else:
    print("~ Basis GSOs are MI!")

#basis=['One','S']
#basisBCs=[[1 for i in range(44)],[1 if i<4 else 0 for i in range(44)]]
#secs=[[1,0],[1,1],[0,0],[0,1]]
#secsBCs=[[1 for i in range(44)],[0 if i<4 else 1 for i in range(44)],\
#         [0 for i in range(44)],[1 if i<4 else 0 for i in range(44)]]
#secsVacEs=[[1*np.dot(secBC[0:4],secBC[0:4])+0.5*np.dot(secBC[4:16],secBC[4:16]),\
#            0.5*np.dot(secBC[16:28],secBC[16:28]) \
#          + np.dot(secBC[28:44],secBC[28:44])] for secBC in secsBCs]
#print(secsVacEs)

#get number of sectors
CBasis = np.zeros(NumBVs)
for i in range(NumBVs):
    for k in range(len(Basis[0][:])):
        if Basis[i][k] % 1 != 0:
            CBasis[i] = CBasis[i]+1
NumSec = 1;
for i in range(NumBVs):
    if CBasis[i] == 0:
        NumSec = NumSec*2
        CBasis[i] = int(2)
    elif CBasis[i] != 0:
        NumSec = NumSec*4
        CBasis[i] = int(4)


#get inner L and R prods to restrict to only M<=0 sec
def vacE(sec):
    a_vL=0.5*np.dot(sec[0:20],sec[0:20])
    a_vR=0.5*np.dot(sec[20:32],sec[20:32]) \
          + np.dot(sec[32:44],sec[32:44])
    return [a_vL,a_vR]

#Want only M<=0 secs but get all sectors first
AllSectorBC = np.zeros((NumSec,Basis.shape[1]))
AllSector = np.zeros((NumSec,NumBVs))
rngs = CBasis.astype(int)
for i,t in enumerate(itertools.product(*[range(i) for i in rngs])):
    AllSectorBC[i,:] = sum([Basis[i,:] * t[i] for i in range(len(t))])
    AllSector[i,:] = t
#now extract just potentially masless ones
NumMSecs = 0
for i in range(NumSec):
    VacESec = vacE(AllSectorBC[i])
    if VacESec[0]>4 or VacESec[1]>8:
        pass
    else:
       NumMSecs+=1 

MSectorBC = np.zeros((NumMSecs,Basis.shape[1]))
MSector = np.zeros((NumMSecs,NumBVs))
rngs = CBasis.astype(int)
for i,t in enumerate(itertools.product(*[range(i) for i in rngs])):
    VacESec = vacE(AllSectorBC[i])
    if VacESec[0]>4 or VacESec[1]>8:
        pass
    else:
        MSectorBC[i,:] = sum([Basis[i,:] * t[i] for i in range(len(t))])
        MSector[i,:] = t
SectorUnRed = MSectorBC.copy()
SectorUnRed[0][:] = 2
MSectorBC = MSectorBC % 2 # would think this needs repeating until in (-1,1] range
MSector[0][0] = 2


# Deltas of the Sectors
SDelta = np.zeros((NumMSecs,1))
for i in range(NumMSecs):
    if MSectorBC[i][0] == 1:
        SDelta[i] = -1
    elif MSectorBC[i][0]  == 0:
        SDelta[i] = 1
        
# GSO phases for the Sectors
SecGSO = np.ones((NumMSecs,NumBVs),dtype=np.complex64)
#maybe can just do the sectors with the basis vecs GSOs
for i in range(NumMSecs):
    for j in range(NumBVs):
        SGSO1 = (SDelta[i]**(np.sum(Basis[j])-1) * SDelta[j]**(np.sum(MSector[i])-1))
        SGSO2 = np.around(np.exp(1j*np.pi*DProd((MSectorBC[i][:]-SectorUnRed[i][:]),Basis[j][:])/2))
        SGSO3 = 1
        for k in range(NumBVs):
            for l in range(NumBVs):
                TSGSO3 = GSO[k][l]**(MSector[i][k]*MSector[j][l])
                SGSO3 = SGSO3 * TSGSO3
        SecGSO[i][j] = SGSO1 * SGSO2 * SGSO3


    

