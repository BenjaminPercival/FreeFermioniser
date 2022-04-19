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
#self parameter is a reference to the current instance of the class- used to access variables in class. can name some other than self but has to be the first parameter of any function in class:
#The __init__() function called automatically every time the class used to create a new object.

#class FFModel:
    #def __init__(self, bas, gsos):
        #self.bas=bas
        #self.gsos=gsos
        
        
def getNumBVs(bas):
    NumBV = bas.shape[0]
    return NumBV


def DProd(vi,vj):
   DP = np.dot(vi[0:4],vj[0:4])+0.5*np.dot(vi[4:16],vj[4:16]) - 0.5*np.dot(vi[16:28],vj[16:28]) \
          - np.dot(vi[28:44],vj[28:44])
   return DP

def DProdLR(vi,vj):
   DPL = np.dot(vi[0:4],vj[0:4])+0.5*np.dot(vi[4:16],vj[4:16]) 
   DPR = 0.5*np.dot(vi[16:28],vj[16:28]) + np.dot(vi[28:44],vj[28:44])
   return [DPL,DPR]

# Dot product of basis vectors
def BasDotProds(bas,GMat):
    BDP=np.array([[DProd(bas[i],bas[k]) for k in range(GMat.shape[1])] for i in range(GMat.shape[1])])
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
    if NBV==Gmat.shape[0]:
        while IsMI is True:
            for i in range(1,NBV):
                if Gmat[i][i]!=-np.around(np.exp(+1j*np.pi*BDPs[i][i]/4)*Gmat[i][0]):
                    printr2("Diagonal GGSO phases not Modular Invariant for i =", i)
                    print("Diagonal GGSO phases not Modular Invariant for i =", i)
                    IsMI=False
                    break
            for jk in it.combinations(range(NBV), 2):
                el=list(jk)
                if Gmat[el[0]][el[1]] != np.around(np.exp(1j*np.pi*BDPs[el[0]][el[1]]/2)*np.conj(Gmat[el[1]][el[0]])):
                    IsMI=False
                    printr2("Diagonal GGSO phases not Modular Invariant for el =", el)
                    print("Diagonal GGSO phases not Modular Invariant for el =", el)
                    #print("this phase is:", Gmat[el[0]][el[1]])
                    #print("but the e^ipi i.j/2 * phase flipped is:", np.around(np.exp(1j*np.pi*BDPs[el[0]][el[1]]/2)*np.conj(Gmat[el[1]][el[0]]))) 
                    break
            break
    else:
        IsMI=False
        print("Number of basis vectors not matched by size of GGSO matrix")
                
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
        break
    if ValidBCs is False:
        RightForm=False
        
    return RightForm
    
def supercurrent(NBV,bas):
    check6=[False for i in range(NBV)]
    ConsistentSCurrent=False
    for i in range(NBV):
        if bas[i][0]==(bas[i][1]+bas[i][4]+bas[i][10])%2:
            if bas[i][0]==(bas[i][1]+bas[i][5]+bas[i][11])%2:
                if bas[i][0]==(bas[i][2]+bas[i][6]+bas[i][12])%2:
                    if bas[i][0]==(bas[i][2]+bas[i][7]+bas[i][13])%2:
                        if bas[i][0]==(bas[i][3]+bas[i][8]+bas[i][14])%2:
                            if bas[i][0]==(bas[i][3]+bas[i][9]+bas[i][15])%2:
                                
                                check6[i]=True
    if check6.count(True)==6:
        ConsistentSCurrent=True
    return ConsistentSCurrent
    
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
        break
    return BasisMI
    
        
def NumbSecs(NBV,lcms):
    N2s=lcms.count(2)
    N4s=lcms.count(4)
    NSecs=2**N2s*4**N4s
    return NSecs
    
#get inner L and R prods to restrict to only M<=0 sec
def vacE(sec):
    a_vL=np.dot(sec[0:4],sec[0:4])+0.5*np.dot(sec[4:16],sec[4:16])
    a_vR=0.5*np.dot(sec[16:28],sec[16:28])+ np.dot(sec[28:44],sec[28:44])
    return [a_vL,a_vR]

#Want only M<=0 secs but get all sectors first
def GetAllSecUnRedBC(NBV,bas,NSects,lcms):
    AllSecUnRedBC = np.zeros((NSects,bas.shape[1]))
    #AllSector = np.zeros((NSects,NBV))
    for i,t in enumerate(it.product(*[range(j) for j in lcms])):
        AllSecUnRedBC[i,:] = sum([bas[i,:] * t[i] for i in range(len(t))])
        #AllSector[i,:] = t
    #SectorUnRed = AllSectorBC.copy()
    #SectorUnRed[0][:] = 2
    #AllSectorBC = AllSectorBC % 2
    #AllSector[0][0] = 2 #NS sec
    return AllSecUnRedBC

def GetAllSecRedBC(SecsUnRed):
    
    SecsUnRed[0][:] = 2
    AllSectorBC = SecsUnRed % 2
    
    return AllSectorBC


def GetAllSecs(NBV,NSects,lcms):
    
    AllSector = np.zeros((NSects,NBV))
    for i,t in enumerate(it.product(*[range(j) for j in lcms])):
        AllSector[i,:] = t
    AllSector[0][0] = 2 #NS sec
    
    return AllSector   

def MasslessUnRedSecs(AllUnRed):
    
    AllUnRed[0][:] = 2
    AllSectorBC = AllUnRed % 2
    NumSec=AllUnRed.shape[0]
    MSecBCsUnRed=[]
    for i in range(NumSec):
        VacESec = vacE(AllSectorBC[i])
        if VacESec[0]>4 or VacESec[1]>8:
            pass
        else:
            MSecBCsUnRed.append(AllSectorBC[i])
    MSecBCsUnRed=np.array(MSecBCsUnRed) 
   
    return MSecBCsUnRed
            

def MasslessSecs(allSecBcRed,AllSecs,NumSec,bas):
    #NMsecs=0
    MSecBCs=[]
    MSecs=[]
    for i in range(NumSec):
        VacESec = vacE(allSecBcRed[i])
        if VacESec[0]>4 or VacESec[1]>8:
            pass
        else:
            MSecBCs.append(allSecBcRed[i])
            MSecs.append(AllSecs[i])
            
            #NMsecs+=1         
    MSecBCs=np.array(MSecBCs)
    MSecs=np.array(MSecs)
    
    return np.hstack((MSecs, MSecBCs))
    
    

# Deltas of the Sectors
def MSecDeltas(MSec,NBV):
    NumMSecs=MSec.shape[0]
    SDelta = np.zeros((NumMSecs,1))
    for i in range(NumMSecs):
        if MSec[i][NBV] == 1:
            SDelta[i] = -1
        elif MSec[i][NBV]  == 0:
            SDelta[i] = 1
    return SDelta
        
# GSO phases for the Sectors
def GetGGSOMSecBas(MSec,MSecUnRed,NBV,delts,bas,ggso):
    NumMSecs=MSec.shape[0]
    SecGSO = np.ones((NumMSecs,NBV),dtype=np.complex64)
    #maybe can just do the sectors with the basis vecs GSOs
    for i in range(NumMSecs):
        for j in range(NBV):
            SGSO1 = (delts[i]**(np.sum(bas[j])-1) * delts[j]**(np.sum(MSec[i][NBV:])-1))
            SGSO2 = np.around(np.exp(1j*np.pi*DProd((MSec[i][NBV:]-MSecUnRed[i][:]),bas[j][:])/2))
            SGSO3 = 1
            for k in range(NBV):
                for l in range(NBV):
                    TSGSO3 = ggso[k][l]**(MSec[i][k]*MSec[j][l])
                    SGSO3 = SGSO3 * TSGSO3
            SecGSO[i][j] = SGSO1 * SGSO2 * SGSO3
    return SecGSO

def MSectVacEs(MSects,NBV):
    NumMSecs=MSects.shape[0]
    MSecVacEs=[vacE(MSects[i][NBV:]) for i in range(NumMSecs)]
    return MSecVacEs


def GetGGSOMSecMSec(MSec,MSecUnRed,NBV,delts,ggso):
    NumMSecs=MSec.shape[0]
    SecGSO = np.ones((NumMSecs,NumMSecs),dtype=np.complex64)
    #maybe can just do the sectors with the basis vecs GSOs
    for i in range(NumMSecs):
        for j in range(NumMSecs):
            SGSO1 = (delts[i]**(np.sum(MSec[j][NBV:])-1) * delts[j]**(np.sum(MSec[i][NBV:])-1))
            SGSO2 = np.around(np.exp(1j*np.pi*DProd(MSec[i][NBV:]-MSecUnRed[i][:],MSec[j][NBV:]-MSecUnRed[j][:])))
            SGSO3 = 1
            for k in range(NBV):
                for l in range(NBV):
                    TSGSO3 = ggso[k][l]**(MSec[i][k]*MSec[j][l])
                    SGSO3 = SGSO3 * TSGSO3
            SecGSO[i][j] = SGSO1 * SGSO2 * SGSO3
    return SecGSO


