#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 16:01:41 2022

@author: wmnc67
"""

import numpy as np
import sys
#from z3 import *
import itertools as it


#define an object:  model with the characteristics: number basis vecs, secs, dot prods,..
def getNumBVs(bas):
    NumBV = bas.shape[0]
    return NumBV


def DProd(vi,vj):
   DP = 0.5*np.dot(vi[0:20],vj[0:20]) - 0.5*np.dot(vi[20:32],vj[20:32]) \
          - np.dot(vi[32:44],vj[32:44])
   return DP

# Dot product of basis vectors
def BasDotProds(bas,GMat):
    BDP=np.array([[DProd(bas[i],bas[k]) for k in np.xrange(GMat.shape[1])] for i in np.xrange(GMat.shape[1])])
    return BDP
  
def printr(thing):
    f = open('OutputModel.txt','a') 
    
    old_stdout = sys.stdout  #  store the default system handler to be able to restore it 

    sys.stdout = f 
    print(thing)
    f.close()
    sys.stdout=old_stdout
    
def printr2(str1,thing):
    f = open('OutputModel.txt','a') 
    
    old_stdout = sys.stdout  #  store the default system handler to be able to restore it 

    sys.stdout = f 
    print(str1)
    print(thing)
    f.close()
    sys.stdout=old_stdout 


# Check if basis GSOs are modular invariant
def IsModInvG(NBV,BDPs,Gmat):
    IsMI=True
    while IsMI is True:
        for i in range(NBV):
            if Gmat[i][i]!=-np.around(np.exp(+1j*np.pi*BDPs[i][i]/4)*Gmat[i][0]):
                printr2("Diagonal GGSO phases not Modular Invariant for i =", i)
                IsMI=False
                break
        for jk in it.combinations(range(NBV), 2):
            el=list(jk)
            if Gmat[el[0]][el[1]] != np.around(np.exp(1j*np.pi*BDPs[el[0]][el[1]]/2)*np.conj(Gmat[el[1]][el[0]])):
                IsMI=False
                printr2("Diagonal GGSO phases not Modular Invariant for el =", el)
                break
                
    #print("~ Basis GSOs are MI!")
    return IsMI



def readBasis(bas):
    RightForm=True
    NCols=np.size(bas, 1)
    if NCols!=44:
        printr("There are not 44 elements per basis vector")
        RightForm=False
    if sum(bas[0,:])!=44:
        printr("1 is not the first basis vector")
        RightForm=False
        
    IsSymmetric=True
    while IsSymmetric is True:
        for i in range(bas.shape[0]):
            for j in range(12):
                if bas[i][4+j]!=bas[i][16+j]: # y^i and w^i's symmetric
                    IsSymmetric=False
                    break
    if IsSymmetric is False:
        RightForm=False
    
    ValidBCs=True
    while ValidBCs is True:
        for i in range(bas.shape[0]):
            for j in range(NCols):
                if j<28:
                    if bas[i][j]!=0 and bas[i][j]!=1: # psi,chi, y,w 's not real BCs
                        printr("Psi^mu, chi12, chi34, chi56 and internal fermions must have real BCs= 0 or 1")
                        ValidBCs=False
                        break
                else:
                    if bas[i][j]!=0 and bas[i][j]!=1 and  bas[i][j]!=0.5: 
                        printr("For the complex gauge fermions we require BCs=0,0.5 or 1")
                        ValidBCs=False
                        break
    if ValidBCs is False:
        RightForm=False
    return RightForm
    

def LCMs(bas):
    Ni=[2]
    for i in range(1,bas.shape[0]):
        Ni.append(2)
        for k in range(np.size(bas, 1)):
                
            if bas[i][k]==0.5: # psi,chi, y,w 's not real BCs
                Ni[i]=4
            else:
                pass
    return Ni
    
def IsModInvB(NBV,BDPs,bas,lcms):
    BasisMI=True
    

    while BasisMI is True:
        for i in range(NBV):
            if lcms[i]*BDPs[i][i]%8!=0:
                BasisMI=False
                printr2("Basis breaks modular invariance for basis vector number =", i)
                break
        for jk in it.combinations(range(NBV), 2):
            el=list(jk)
            Nij=max(lcms[el[0]],lcms[el[1]])
            if Nij*BDPs[el[0]][el[1]]%4!=0:
                BasisMI=False
                printr2("Basis breaks modular invariance for el =", el)
                break
    return BasisMI
    
        
def NumbSecs(NBV,lcms):
    N2s=lcms.count(2)
    N4s=lcms.count(4)
    NSecs=2**N2s*4**N4s
    return NSecs
    
#get inner L and R prods to restrict to only M<=0 sec
def vacE(sec):
    a_vL=0.5*np.dot(sec[0:16],sec[0:16])
    a_vR=0.5*np.dot(sec[16:28],sec[16:28]) \
          + np.dot(sec[28:44],sec[28:44])
    return [a_vL,a_vR]

#Want only M<=0 secs but get all sectors first
def GetAllSecs(NBV,bas,NSects,lcms):
    AllSectorBC = np.zeros((NSects,bas.shape[1]))
    AllSector = np.zeros((NSects,NBV))
    
    
# =============================================================================
# rngs = CBasis.astype(int)
# for i,t in enumerate(itertools.product(*[range(i) for i in rngs])):
#     AllSectorBC[i,:] = sum([Basis[i,:] * t[i] for i in range(len(t))])
#     AllSector[i,:] = t
# #now extract just potentially masless ones
# NumMSecs = 0
# for i in range(NumSec):
#     VacESec = vacE(AllSectorBC[i])
#     if VacESec[0]>4 or VacESec[1]>8:
#         pass
#     else:
#        NumMSecs+=1 
# 
# MSectorBC = np.zeros((NumMSecs,Basis.shape[1]))
# MSector = np.zeros((NumMSecs,NumBVs))
# rngs = CBasis.astype(int)
# for i,t in enumerate(itertools.product(*[range(i) for i in rngs])):
#     VacESec = vacE(AllSectorBC[i])
#     if VacESec[0]>4 or VacESec[1]>8:
#         pass
#     else:
#         MSectorBC[i,:] = sum([Basis[i,:] * t[i] for i in range(len(t))])
#         MSector[i,:] = t
# SectorUnRed = MSectorBC.copy()
# SectorUnRed[0][:] = 2
# MSectorBC = MSectorBC % 2 # would think this needs repeating until in (-1,1] range
# MSector[0][0] = 2
# =============================================================================


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


    

