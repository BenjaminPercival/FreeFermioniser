#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 17:25:50 2022

@author: wmnc67
"""

from FFerGetModelDetails import DProdLR, printr, printr2 #DProd, 
from z3 import *
import numpy as np
import sys
np.set_printoptions(threshold=sys.maxsize)
#NSSec()
#SSec() #in this function account for possible Stilde
#"Model's spacetime SUSY is N=", ModSusy()
#enhancements()
#ObservableMassless()
#OnShellTachyons()
#Exotics()
#Hidden()

def ConvertRFtoStr(RFInds):
    RFStrs=['Psi^mu','Chi12','Chi34','Chi56','y1yb1','y2yb2','y3yb3','y4yb4','y5yb5','y6yb6',\
            'w1wb1','w2wb2','w3wb3','w4wb4','w5wb5','w6wb6','y/wbar error','y/wbar error',\
            'y/wbar error','y/wbar error','y/wbar error','y/wbar error','y/wbar error',\
            'y/wbar error','y/wbar error','y/wbar error','y/wbar error','y/wbar error',\
            'Psi1','Psi2','Psi3','Psi4','Psi5','Eta1','Eta2','Eta3',\
            'Phi1','Phi2','Phi3','Phi4','Phi5','Phi6','Phi7','Phi8']
    RFS=[RFStrs[i] for i in RFInds]
    return RFS
    

def UnprojectedSecs(NBV,MSects,MMSectGSOs,deltMSects,bas,MScVacEs):
    NumMSecs=MSects.shape[0]
    #MSecVacEs=[vacE(MSects[i][NBV:]) for i in range(NumMSecs)]
    #MSecs48UnprojBool=[]
    MSecs48Unproj=[]
    MSec48UnProjInds=[]
    #MSecs44UnprojBool=[]
    MSecs44Unproj=[]
    MSec44UnProjInds=[]
    #MSecs08UnprojBool=[]
    MSecs08Unproj=[]
    MSec08UnProjInds=[]
    MSecs04Unproj=[]
    MSec04UnProjInds=[]
    MSecs40Unproj=[]
    for i in range(1,NumMSecs):#don't need to do NS
        vacEi=MScVacEs[i]
        SecUnProjd=True
        #print("For sector: ", MSects[i])
        #print("vacE is: ", vacEi)
        if vacEi==[4,8]:
            for j in range(1,NumMSecs):#don't need 1
                if DProdLR(MSects[i][NBV:],MSects[j][NBV:])==[0,0]:
                    LHS=deltMSects[i]*MMSectGSOs[i][j]
                    if LHS==-1:
                        SecUnProjd=False
                        #print("project by vec j=", MSects[j])
            if SecUnProjd is True:
                MSec48UnProjInds.append(i)
                MSecs48Unproj.append(MSects[i])
                
        elif vacEi==[4,4]:
            oscillsProjd44=[]
            for j in range(1,NumMSecs):#don't need 1
                if DProdLR(MSects[i][NBV:],MSects[j][NBV:])==[0,0]:
                    for k in range(28,44):
                        LHS=deltMSects[i]*MMSectGSOs[i][j]*(-1)**(MSects[j][NBV+k])
                        if LHS==-1:
                            oscillsProjd44.append(k)
            oscillsUnProjd44=[osc for osc in range(28,44) if osc not in oscillsProjd44]
            if len(oscillsUnProjd44)!=0:
                MSec44UnProjInds.append(i)
                MSecs44Unproj.append(MSects[i])
                MSecs44Unproj.append(oscillsUnProjd44)
                
                
        elif vacEi==[0,8]:
            oscillsProjd08=[]
            print("(0,8) sec: ", MSects[i][:NBV])
            for j in range(1,NumMSecs):#don't need 1
                if DProdLR(MSects[i][NBV:],MSects[j][NBV:])==[0,0]:
                    for k in range(0,4):
                        LHS=deltMSects[i]*MMSectGSOs[i][j]*(-1)**(MSects[j][NBV+k])
                        if LHS==-1:
                            oscillsProjd08.append(k)
                        if k==0:
                            if LHS==-1:
                                print("Projected by sec: ", MSects[j][:NBV])
                                print("GSO of sec i with sec j: ", MMSectGSOs[i][j])
                                print("print osc bit: ", (-1)**(MSects[j][NBV+k]))
            oscillsUnProjd08=[osc for osc in range(0,4) if osc not in oscillsProjd08]
            if len(oscillsUnProjd08)!=0:
                MSec08UnProjInds.append(i)
                MSecs08Unproj.append(MSects[i])
                MSecs08Unproj.append(oscillsUnProjd08)
                
        elif vacEi==[0,4]:
            oscillsProjd04=[]
            for j in range(1,NumMSecs):#don't need 1
                if DProdLR(MSects[i][NBV:],MSects[j][NBV:])==[0,0]:
                    for k in range(0,4):
                        for l in range(28,44):
                            LHS=deltMSects[i]*MMSectGSOs[i][j]*(-1)**(MSects[j][NBV+k])*(-1)**(MSects[j][NBV+l])
                            if LHS==-1:
                                oscillsProjd04.append([k,l])
                                #oscillsProjd04.append(l)
            oscillsUnProjd04=[[osck,oscl] for osck in range(28,44) for oscl in range(28,44) if [osck,oscl] not in oscillsProjd04]
            #oscillsUnProjd04=[osc for osc in range(28,44) if osc not in oscillsProjd04]
            if len(oscillsUnProjd04)!=0:
                MSec04UnProjInds.append(i)
                MSecs04Unproj.append(MSects[i])
                MSecs04Unproj.append(oscillsUnProjd04)
                
        elif vacEi==[4,0]:
            oscillsProjd40=[]
            for j in range(1,NumMSecs):#don't need 1
                if DProdLR(MSects[i][NBV:],MSects[j][NBV:])==[0,0]:
                    for k in range(28,44):
                        for l in range(28,44):
                            LHS=deltMSects[i]*MMSectGSOs[i][j]*(-1)**(MSects[j][NBV+k])
                            if LHS==-1:
                                oscillsProjd40.append([k,l])
                                #oscillsProjd40.append(l)
                                
            
            oscillsUnProjd40=[[osck,oscl] for osck in range(28,44) for oscl in range(28,44) if [osck,oscl] not in oscillsProjd40]
            if len(oscillsUnProjd40)!=0:
                MSecs40Unproj.append(MSects[i])
                MSecs40Unproj.append(oscillsUnProjd40)
                
        
    return [MSecs48Unproj, MSec48UnProjInds, MSecs44Unproj, MSec44UnProjInds, MSecs08Unproj, MSec08UnProjInds, MSecs04Unproj, MSec04UnProjInds, MSecs40Unproj]
            
            
def printUnProjdSecs(MScts,UnProjdScs,NBV):
    
    MSecs48=UnProjdScs[0]#[MSecs48Unproj, MSecs44Unproj, MSecs08Unproj, MSecs04Unproj, MSecs40Unproj]
    MSecs44=UnProjdScs[2]
    MSecs08=UnProjdScs[4]
    MSecs04=UnProjdScs[6]
    MSecs40=UnProjdScs[8]
    
    
    printr("Massless Sectors:")
    printr2("NS Sector:", MScts[0][:NBV])
    printr("Sectors of Vacuum Energy=(4,0) (S Sector?):")
    for i in range(len(MSecs40)):
        if i%2==0:
            printr2("Sector: ", MSecs40[i][:NBV])
        else:
            #printr2("oscillators: ", MSecs40[i])
            printr("oscillators: ...")
    printr("Sectors of Vacuum Energy=(4,8):")
    for sec48 in MSecs48:
        printr(sec48[:NBV])
    printr("Sectors of Vacuum Energy=(4,4):")
    for i in range(len(MSecs44)):
        if i%2==0:
            printr2("Sector: ", MSecs44[i][:NBV])
        else:
            printr2("oscillators: ", MSecs44[i][:NBV])
    printr("Sectors of Vacuum Energy=(0,8):")
    for i in range(len(MSecs08)):
        if i%2==0:
            printr2("Sector: ", MSecs08[i][:NBV])
        else:
            printr2("oscillators: ", MSecs08[i][:NBV])
    printr("Sectors of Vacuum Energy=(0,4):")
    for i in range(len(MSecs04)):
        if i%2==0:
            printr2("Sector: ", MSecs04[i][:NBV])
        else:
            printr2("oscillators: ", MSecs04[i][:NBV])            
    
    return None


    
def NSSec(NBvs,bas):
    #secBC=[0 for i in range(44)]
    printr("NS sector:")# 1 left and two right-moving NS oscillators or 1 dX^mu")
    #oscillsL=[index for index, char in enumerate(secBC) if index<4 if char==0]
    #oscillsR=[index for index, char in enumerate(secBC) if index>16 if char==0]
    printr("Along with the Gravitational states: psi^mu dXbar^mu |NS>, we have:")
    printr("[L Osc, R Osc1, Freq 1, R Osc2, Freq 2]")
    for OscInd0 in range(0,4):
        for OscInd1 in range(28,44):
            for freq1 in range(2):
                for OscInd2 in range(28,44):
                    for freq2 in range(2):
                        for j in range(NBvs):
                        #for basisBC in Basis:
                            OscBit=(-1)**(bas[j][OscInd0])*(-1)**(freq1*bas[j][OscInd1])*(-1)**(freq2*bas[j][OscInd2])
                            if OscBit==1: #LHS=+1 necessarily.
                                #print("L osc, R oscs: ")
                                printr([OscInd0,OscInd1,freq1,OscInd2,freq2])
    

def Massless40(bas,NBV,ggso,UnProjdScs):  #S Sector
    
    #secBC=[1 if i<4 else 0 for i in range(44)]
    MSecs40=UnProjdScs[8]
    NSecs40=len(MSecs40)
    NSUSY=10
    
    for i in range(NSecs40):
        if i%2==0:
            if np.sum(MSecs40[i][NBV:NBV+4])==4:
            
                printr("S sector:") # two right-moving NS oscillators or 1 dX^mu")
                for oscS in MSecs40[i+1]:
                    #print(oscS)
                    b=[Bool('b%s' % (i)) for i in range(4)]
                    s=Solver()
                    OscInd1=oscS[0]
                    OscInd2=oscS[1]
                    for freq1 in range(2):
                        for freq2 in range(2):
                            for i in range(NBV):
                                LHS=(-1)*ggso[1][i]
                                RRs=[]
                                for j in range(4):
                                    if bas[i][j]==1:
                                        RRs.append(j)
                                if len(RRs)==0 and ggso[1][i]==1: # projects
                                    NSUSY=0
                                    printr("S gravitini projected- Model is non-supersymmetric")
                                
                                
                                    #printr2("Oscills: ") 
                    
                                OscBit=(-1)**(freq1*bas[i][OscInd1])*(-1)**(freq2*bas[i][OscInd2])
                                fRR=[b[k] for k in RRs] #ramond fermions to constrain
                                #print("fRR: ", fRR)
                                if LHS*OscBit==1:
                                    s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==0) 
                                else:
                                    s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==1) 
                                #SUSYcount=0
                            printr([OscInd1,freq1,OscInd2,freq2])
                            while s.check() == sat: 
                                #if s.check() == sat: 
                                m = s.model()
        # =============================================================================
        #                                 for v in vars:
        #                                   print("%s = %5s" % (v, m.eval(v, model_completion = True))),
        #                                 print
        #                                 s.add(Or([p != v for p, v in [(v, m.eval(v, model_completion = True)) for v in vars]]))
        # =============================================================================
                                stateB=[m.eval(b[i], model_completion = True) for i in range(4)]
                                #state=[ -1 if item.sexpr()=='true' else 0 for item in stateB]
                                
                                printr(stateB)
                                
                                SUSYcount+=1
                                #if countr%100==0:
                                #    print(countr)
                                s.add(Not(And([v() == m.eval(v, model_completion = True) for v in m]))) 
                                    
                            s.reset()
                                        
                         
            else:
                printr("Other (4,0) sector(?)- error?")
    
    #s.reset()
    if NSUSY!=0:
        Sket=[Bool('S%s' % (i)) for i in range(4)]
    #gravitino states:
        s=Solver()
        s.reset()
        for i in range(NBV):
            LHS=(-1)*ggso[1][i]
            RRs=[]
            for j in range(4):
                if bas[i][j]==1:
                    RRs.append(i)
            if len(RRs)!=0:
                #print("RRs: ", RRs)
                SRR=[Sket[k] for k in range(len(Sket))] #ramond fermions to constrain
                
                if ggso[1][i]==1:
                    s.add(Sum([If(SRR[i],1,0) for i in range(len(SRR))])%2==1) 
                else:
                    s.add(Sum([If(SRR[i],1,0) for i in range(len(SRR))])%2==0)
                #if s.check() == sat: 
        printr("Gravitino: dXbar^mu |S >")
                #NumOscStates+=1
        countr=0
        while s.check() == sat: 
        #if s.check() == sat: 
            
            m = s.model () 
            
            stateB=[m[Sket[i]] for i in range(len(Sket))]
            #state=[ -1 if item.sexpr()=='true' else 0 for item in stateB]
            printr(stateB)
            countr+=1
            #if countr%100==0:
            #    print(countr)
            s.add(Not(And([v() == m[v] for v in m]))) 
        s.reset()
        NSUSY=countr/2
        printr2("Number of gravitini", NSUSY)
            
    return NSUSY

#def enhancements():
    
def Massless48(bas,NBV,ggsoMBs,UnProjdScs):
    MSecs48=UnProjdScs[0]
    MSecs48Inds=UnProjdScs[1]
    NSecs48=len(MSecs48)
    trU11=0
    trU12=0
    trU13=0
    printr("Massless Sectors Vacuum Energy= (4,8)-> spinorials")
    for i in range(NSecs48):
         
        #RFsX=[]
        RFsX=[index for index, char in enumerate(MSecs48[i][NBV:]) if char == 1]
        RFs=[index for index in RFsX if index in range(16) or index in range(28,44)]#symmetric pairings assumed
        #if BC==1:
            #print(secBC.index(BC))
            #RFs.append(secBC.index(BC))
        #printr(RFs)
        #RFlen=len(RFs)
        b=[Bool('b%s' % (i)) for i in range(44)]
        f = [b[i] for i in RFs] 
        #print("f:", f)
        s=Solver()
        
        #GGSO eqn
    
        #print("For sec: ", MSecs48[i][:NBV])
        #print(LHS)
        for j in range(NBV):#just need to loop through basis vecs
            #print("secProj:", secProj)
            #totalRR=DProd(secBC,secProj)
            LHS=(-1)**(MSecs48[i][NBV])*ggsoMBs[MSecs48Inds[i]][j]
            RRsX=[]
            
            for k in range(44):
                if bas[j][k]==1 and MSecs48[i][NBV+k]==1:
                    RRsX.append(k)
            RRs=[index for index in RRsX if index in range(16) or index in range(28,44)]#symmetric pairings assumed
            #print("RRs: ", RRs)
            fRR=[b[l] for l in RRs] #ramond fermions to constrain
            #print("fRR: ", fRR)
            if np.real(LHS)==1:#assuming 
                s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==0) 
            elif np.real(LHS)==-1:
                s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==1) 
            else:
                print("complex GGSO phases not implemented yet / error in LHS variable")
            #print("for basis vec: ", j)
            #print("LHS is: ", LHS)
        #print(s)
        #print(s.check())
        stateArray=[]
        if s.check() == sat:
            printr(MSecs48[i][:NBV])
            #printr(RFs)
            printr(np.asarray(ConvertRFtoStr(RFs)))
            stateArray.append(RFs)
            if MSecs48[i][NBV]==1: #fermionic
                if np.sum(MSecs48[i][NBV:NBV+4])==2 and np.sum(MSecs48[i][NBV+28:NBV+36])==6:
                    #F^1,2,3 - 16/16bars uniquely(?)- practically yes but, e.g., {psi1234,eta23} allowed currently
                    printr("Sector gives fermion generation(s):")
        
        
        while s.check() == sat: 
            #m = s.model () 
            
            #stateB=[m.eval(f[i], model_completion = True) for i in range(len(RFs))]
            #state=[ -1 if item.sexpr()=='true' else 0 for item in stateB]
            #printr(stateB)
            #countr+=1
            #if countr%100==0:
            #    print(countr)
            #s.add(Not(And([v() == m[v] for v in m]))) 
            #s.add(Not(And([v() == m.eval(v, model_completion = True) for v in m]))) 
            m = s.model()
            #for v in f:
            stateB=[m.eval(f[i], model_completion = True) for i in range(len(RFs))]
            state=[ -1 if item.sexpr()=='true' else 0 for item in stateB]
            #printr(state)
            stateArray.append(state)
            #print
            s.add(Or([p != v for p, v in [(v, m.evaluate(v, model_completion = True)) for v in f]]))
        if len(stateArray)!=0:
            stateArr=np.asarray(stateArray)
            printr(stateArr)
            for j in range(len(RFs)):
                if stateArr[0][j]==33:
                    trU11+=np.sum(stateArr[1:][j])
                elif stateArr[0][j]==34:
                    trU12+=np.sum(stateArr[1:][j])
                elif stateArr[0][j]==35:
                    trU13+=np.sum(stateArr[1:][j])
        s.reset()
    return [trU11,trU12,trU13]

def Massless44(bas,NBV,ggsoMBs,UnProjdScs):
    MSecs44=UnProjdScs[2]
    TwistedHiggs=False
    MSecs44Inds=UnProjdScs[3]
    NSecs44=len(MSecs44)
    #print(NSecs44)
    trU11=0
    trU12=0
    trU13=0
    printr("Massless Sectors Vacuum Energy= (4,4)-> Vectorials with R-moving oscillator")
    for i in range(NSecs44):
        if i%2==0:
            RFsX=[index for index, char in enumerate(MSecs44[i][NBV:]) if char == 1]
            RFs=[index for index in RFsX if index in range(16) or index in range(28,44)]#symmetric pairings assumed
            #if BC==1:
                #print(secBC.index(BC))
                #RFs.append(secBC.index(BC))
            #print("for sec: ", MSecs44[i][:NBV])
            b=[Bool('b%s' % (n)) for n in range(44)]
            f = [b[m] for m in RFs] 
            #print("f:", f)
            s=Solver()
            #oscRs=MSecs44[i+1]#already collected oscills in unproj function...needs to be correct- freq not done?
            for oscR in range(28,44):
                for freq in range(2):
                    #print(oscR)
                    Oscbit=complex((-1)**(freq*MSecs44[i][NBV+oscR]))
                    #print(Oscbit)
                    for j in range(NBV):
                        LHS=(-1)**(MSecs44[i][NBV])*ggsoMBs[MSecs44Inds[int(i/2)]][j] #better to write the Oscbit as a complex number
                        #print(LHS)
                        RRsX=[]
                        
                        for k in range(44):
                            if bas[j][k]==1 and MSecs44[i][NBV+k]==1:
                                RRsX.append(k)
                        RRs=[index for index in RRsX if index in range(16) or index in range(28,44)]#symmetric pairings assumed
            
                        #print("RRs: ", RRs)
                        fRR=[b[l] for l in RRs] #ramond fermions to constrain
                        #print("fRR: ", fRR)
                        if np.real(LHS)==np.real(Oscbit) and np.imag(LHS)==np.imag(Oscbit):
                            s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==0) 
                        else:
                            s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==1) 
                    #be good to check for higgs here
                        #if oscR==37:
                            #print("oscR=37 and basis vec: ", j)
                            #print("LHS is: ", LHS)
                            #print("Oscbit is: ", Oscbit)
                    #print(s.check())
                    stateArray=[]
                    if s.check() == sat: 
                        printr(MSecs44[i][:NBV])
                        #printr(RFs)
                        printr(np.asarray(ConvertRFtoStr(RFs)))
                        stateArray.append(RFs)
                        printr([oscR,freq])
                    
                    while s.check() == sat: 
                        m = s.model()
                        #for v in f:
                        stateB=[m.eval(f[i], model_completion = True) for i in range(len(RFs))]
                        state=[ -1 if item.sexpr()=='true' else 0 for item in stateB]
                        stateArray.append(state)
                        #print
                        s.add(Or([p != v for p, v in [(v, m.evaluate(v, model_completion = True)) for v in f]]))
                    if len(stateArray)!=0:
                        stateArr=np.asarray(stateArray)
                        printr(stateArr)
                        for j in range(len(RFs)):
                            if stateArr[0][j]==33:
                                trU11+=np.sum(stateArr[1:][j])
                            elif stateArr[0][j]==34:
                                trU12+=np.sum(stateArr[1:][j])
                            elif stateArr[0][j]==35:
                                trU13+=np.sum(stateArr[1:][j])
                    s.reset()
    return [trU11,trU12,trU13,TwistedHiggs]



def Massless04(bas,NBV,ggsoMBs,UnProjdScs):
    MSecs04=UnProjdScs[4]
    Enhancement04=False
    MSecs04Inds=UnProjdScs[5]
    NSecs04=len(MSecs04)
    #print(NSecs44)
    trU11=0
    trU12=0
    trU13=0
    printr("Massless Sectors Vacuum Energy= (0,4)-> L and R-moving oscillators, possible enhancments")
    for i in range(NSecs04):
        if i%2==0:
            RFsX=[index for index, char in enumerate(MSecs04[i][NBV:]) if char == 1]
            RFs=[index for index in RFsX if index in range(16) or index in range(28,44)]#symmetric pairings assumed
            #if BC==1:
                #print(secBC.index(BC))
                #RFs.append(secBC.index(BC))
            #print("for sec: ", MSecs04[i][:NBV])
            b=[Bool('b%s' % (n)) for n in range(44)]
            f = [b[m] for m in RFs] 
            #print("f:", f)
            s=Solver()
            #oscRs=MSecs44[i+1]#already collected oscills in unproj function...needs to be correct- freq not done?
            for oscL in range(0,4):
                for oscR in range(28,44):
                    for freq in range(2):
                        #print(oscR)
                        Oscbit=complex((-1)**(MSecs04[i][NBV+oscL])*(-1)**(freq*MSecs04[i][NBV+oscR]))
                        #print(Oscbit)
                        for j in range(NBV):
                            LHS=(-1)**(MSecs04[i][NBV])*ggsoMBs[MSecs04Inds[int(i/2)]][j] #better to write the Oscbit as a complex number
                            #print(LHS)
                            RRsX=[]
                            
                            for k in range(44):
                                if bas[j][k]==1 and MSecs04[i][NBV+k]==1:
                                    RRsX.append(k)
                            RRs=[index for index in RRsX if index in range(16) or index in range(28,44)]#symmetric pairings assumed
            
                            #print("RRs: ", RRs)
                            fRR=[b[l] for l in RRs] #ramond fermions to constrain
                            #print("fRR: ", fRR)
                            if np.real(LHS)==np.real(Oscbit) and np.imag(LHS)==np.imag(Oscbit):
                                s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==0) 
                            else:
                                s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==1) 
                        #be good to check for higgs here
                            #if oscR==37:
                                #print("oscR=37 and basis vec: ", j)
                                #print("LHS is: ", LHS)
                                #print("Oscbit is: ", Oscbit)
                        #print(s.check())
                        stateArray=[]
                        if s.check() == sat: 
                            printr(MSecs04[i][:NBV])
                            #printr(RFs)
                            printr(np.asarray(ConvertRFtoStr(RFs)))
                            stateArray.append(RFs)
                            printr([oscL,oscR,freq])
                        
                        while s.check() == sat: 
                            m = s.model()
                            #for v in f:
                            stateB=[m.eval(f[i], model_completion = True) for i in range(len(RFs))]
                            state=[ -1 if item.sexpr()=='true' else 0 for item in stateB]
                            stateArray.append(state)
                            #print
                            s.add(Or([p != v for p, v in [(v, m.evaluate(v, model_completion = True)) for v in f]]))
                        if len(stateArray)!=0:
                            stateArr=np.asarray(stateArray)
                            printr(stateArr)
                            for j in range(len(RFs)):
                                if stateArr[0][j]==33:
                                    trU11+=np.sum(stateArr[1:][j])
                                elif stateArr[0][j]==34:
                                    trU12+=np.sum(stateArr[1:][j])
                                elif stateArr[0][j]==35:
                                    trU13+=np.sum(stateArr[1:][j])
                        s.reset()
    return [trU11,trU12,trU13,Enhancement04]

def Massless08(bas,NBV,ggsoMBs,UnProjdScs):
    MSecs08=UnProjdScs[6]
    Enhancement08=False
    MSecs08Inds=UnProjdScs[7]
    NSecs08=len(MSecs08)
    trU11=0
    trU12=0
    trU13=0
    
    #print(NSecs44)
    printr("Massless Sectors Vacuum Energy= (0,8)-> L-moving Oscillator, possible enhancements")
    for i in range(NSecs08):
        if i%2==0:
            
            RFsX=[index for index, char in enumerate(MSecs08[i][NBV:]) if char == 1]
            RFs=[index for index in RFsX if index in range(16) or index in range(28,44)]#symmetric pairings assumed
            #if BC==1:
                #print(secBC.index(BC))
                #RFs.append(secBC.index(BC))
            print("for sec: ", MSecs08[i][:NBV])
            b=[Bool('b%s' % (n)) for n in range(44)]
            f = [b[m] for m in RFs] 
            #print("f:", f)
            s=Solver()
            #oscRs=MSecs44[i+1]#already collected oscills in unproj function...needs to be correct- freq not done?
            for oscL in range(0,4):
                
                #print(oscR)
                Oscbit=complex((-1)**(MSecs04[i][NBV+oscL]))
                #print(Oscbit)
                for j in range(NBV):
                    LHS=(-1)**(MSecs04[i][NBV])*ggsoMBs[MSecs08Inds[int(i/2)]][j] #better to write the Oscbit as a complex number
                    #print(LHS)
                    RRsX=[]
                    
                    for k in range(44):
                        if bas[j][k]==1 and MSecs08[i][NBV+k]==1:
                            RRsX.append(k)
                    RRs=[index for index in RRsX if index in range(16) or index in range(28,44)]#symmetric pairings assumed
            
                    #print("RRs: ", RRs)
                    fRR=[b[l] for l in RRs] #ramond fermions to constrain
                    #print("fRR: ", fRR)
                    if np.real(LHS)==np.real(Oscbit) and np.imag(LHS)==np.imag(Oscbit):
                        s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==0) 
                    else:
                        s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==1) 
                #be good to check for higgs here
                    #if oscR==37:
                        #print("oscR=37 and basis vec: ", j)
                        #print("LHS is: ", LHS)
                        #print("Oscbit is: ", Oscbit)
                #print(s.check())
                
                stateArray=[]
                if s.check() == sat: 
                    printr(MSecs08[i][:NBV])
                    #printr(RFs)
                    printr(np.asarray(ConvertRFtoStr(RFs)))
                    stateArray.append(RFs)
                    #stateArray.append(ConvertRFtoStr(RFs))
                    printr([oscL,oscR,freq])
                
                while s.check() == sat: 
                    m = s.model()
                    #for v in f:
                    stateB=[m.eval(f[i], model_completion = True) for i in range(len(RFs))]
                    state=[ -1 if item.sexpr()=='true' else 0 for item in stateB]
                    stateArray.append(state)
                    #print
                    s.add(Or([p != v for p, v in [(v, m.evaluate(v, model_completion = True)) for v in f]]))
                if len(stateArray)!=0:
                    stateArr=np.asarray(stateArray)
                    printr(stateArr)
                    for j in range(len(RFs)):
                        if stateArr[0][j]==33:
                            trU11+=np.sum(stateArr[1:][j])
                        elif stateArr[0][j]==34:
                            trU12+=np.sum(stateArr[1:][j])
                        elif stateArr[0][j]==35:
                            trU13+=np.sum(stateArr[1:][j])
                s.reset()
    return [trU11,trU12,trU13,Enhancement08]

#from operator import add
#better doing with numpy 
def sumU1s(U1s48,U1s44,U1s04,U1s08):
    U1sAll=list(map(lambda w,x,y,z:w+x+y+z, U1s48, U1s44,U1s04,U1s08))
    printr2("Tr U(1)_1,2,3= ", U1sAll)
    return U1sAll




