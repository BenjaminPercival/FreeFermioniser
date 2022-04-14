#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  7 17:25:50 2022

@author: wmnc67
"""

from FFerGetModelDetails import DProd, printr, printr2

#NSSec()
#SSec() #in this function account for possible Stilde
#"Model's spacetime SUSY is N=", ModSusy()
#enhancements()
#ObservableMassless()
#OnShellTachyons()
#Exotics()
#Hidden()

def UnprojectedSecs(NBV,MSects,MSectGSOs,deltMSects,bas,MScVacEs):
    NumMSecs=MSects.shape[0]
    #MSecVacEs=[vacE(MSects[i][NBV:]) for i in range(NumMSecs)]
    #MSecs48UnprojBool=[]
    MSecs48Unproj=[]
    #MSecs44UnprojBool=[]
    MSecs44Unproj=[]
    #MSecs08UnprojBool=[]
    MSecs08Unproj=[]
    MSecs04Unproj=[]
    MSecs40Unproj=[]
    for i in range(1,NumMSecs):#don't need to do NS
        vacEi=MScVacEs[i]
        SecUnProjd=True
        print("For sector: ", MSects[i])
        print("vacE is: ", vacEi)
        if vacEi==[4,8]:
            for j in range(1,NBV):#don't need 1
                if DProd(MSects[i][NBV:],bas[j][:])==0:
                    LHS=deltMSects[i]*MSectGSOs[i][j]
                    if LHS==-1:
                        SecUnProjd=False
                        print("project by basis vec j=", j)
            if SecUnProjd is True:
                MSecs48Unproj.append(MSects[i])
        elif vacEi==[4,4]:
            oscillsProjd44=[]
            for j in range(1,NBV):#don't need 1
                if DProd(MSects[i][NBV:],bas[j][:])==0:
                    for k in range(28,44):
                        LHS=deltMSects[i]*MSectGSOs[i][j]*(-1)**(bas[j][k])
                        if LHS==-1:
                            oscillsProjd44.append(k)
            oscillsUnProjd44=[osc for osc in range(28,44) if osc not in oscillsProjd44]
            if len(oscillsUnProjd44)!=0:
                MSecs44Unproj.append(MSects[i])
                MSecs44Unproj.append(oscillsUnProjd44)
                
                
        elif vacEi==[0,8]:
            oscillsProjd08=[]
            for j in range(1,NBV):#don't need 1
                if DProd(MSects[i][NBV:],bas[j][:])==0:
                    for k in range(0,4):
                        LHS=deltMSects[i]*MSectGSOs[i][j]*(-1)**(bas[j][k])
                        if LHS==-1:
                            oscillsProjd08.append(k)
            oscillsUnProjd08=[osc for osc in range(0,4) if osc not in oscillsProjd08]
            if len(oscillsUnProjd08)!=0:
                MSecs08Unproj.append(MSects[i])
                MSecs08Unproj.append(oscillsUnProjd08)
                
        elif vacEi==[0,4]:
            oscillsProjd04=[]
            for j in range(1,NBV):#don't need 1
                if DProd(MSects[i][NBV:],bas[j][:])==0:
                    for k in range(0,4):
                        for l in range(28,44):
                            LHS=deltMSects[i]*MSectGSOs[i][j]*(-1)**(bas[j][k])*(-1)**(bas[j][l])
                            if LHS==-1:
                                oscillsProjd04.append([k,l])
                                #oscillsProjd04.append(l)
            oscillsUnProjd04=[[osck,oscl] for osck in range(28,44) for oscl in range(28,44) if [osck,oscl] not in oscillsProjd04]
            #oscillsUnProjd04=[osc for osc in range(28,44) if osc not in oscillsProjd04]
            if len(oscillsUnProjd04)!=0:
                MSecs04Unproj.append(MSects[i])
                MSecs04Unproj.append(oscillsUnProjd04)
                
        elif vacEi==[4,0]:
            oscillsProjd40=[]
            for j in range(1,NBV):#don't need 1
                if DProd(MSects[i][NBV:],bas[j][:])==0:
                    for k in range(28,44):
                        for l in range(28,44):
                            LHS=deltMSects[i]*MSectGSOs[i][j]*(-1)**(bas[j][k])
                            if LHS==-1:
                                oscillsProjd40.append([k,l])
                                #oscillsProjd40.append(l)
            
            oscillsUnProjd40=[[osck,oscl] for osck in range(28,44) for oscl in range(28,44) if [osck,oscl] not in oscillsProjd40]
            if len(oscillsUnProjd40)!=0:
                MSecs40Unproj.append(MSects[i])
                MSecs40Unproj.append(oscillsUnProjd40)
                
        
    return [MSecs48Unproj, MSecs44Unproj, MSecs08Unproj, MSecs04Unproj, MSecs40Unproj]
            
            
def printUnProjdSecs(MScts,UnProjdScs,NBV):
    
    MSecs48=UnProjdScs[0]#[MSecs48Unproj, MSecs44Unproj, MSecs08Unproj, MSecs04Unproj, MSecs40Unproj]
    MSecs44=UnProjdScs[1]
    MSecs08=UnProjdScs[2]
    MSecs04=UnProjdScs[3]
    MSecs40=UnProjdScs[4]
    
    
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

    
def NSSec(NBvs):
    #secBC=[0 for i in range(44)]
    printr("NS sector:")# 1 left and two right-moving NS oscillators or 1 dX^mu")
    #oscillsL=[index for index, char in enumerate(secBC) if index<4 if char==0]
    #oscillsR=[index for index, char in enumerate(secBC) if index>16 if char==0]
    
    for OscInd0 in range(0,4):
        for OscInd1 in range(27,44):
            for freq1 in range(2):
                for OscInd2 in range(27,44):
                    for freq2 in range(2):
                        for j in range(NBvs):
                        #for basisBC in Basis:
                            OscBit=(-1)**(Basis[j][OscInd0])*(-1)**(freq1*Basis[j][OscInd1])*(-1)**(freq2*Basis[j][OscInd2])
                            if OscBit==1: #LHS=+1 necessarily.
                                #print("L osc, R oscs: ")
                                printr([OscInd0,OscInd1,freq1,OscInd2,freq2])
    printr("Gravitational states: psi^mu dXbar^mu |NS>")

def SSec(): 
    secBC=[1 if i<4 else 0 for i in range(44)]
    NSUSY=10
    #if S in basis:
    printr("S sector:") # two right-moving NS oscillators or 1 dX^mu")
    RFs=[]
    RFs=[index for index, char in enumerate(secBC) if char == 1]
        #if BC==1:
            #print(secBC.index(BC))
            #RFs.append(secBC.index(BC))
    #print(RFs)
    #RFlen=len(RFs)
    oscills=[index for index, char in enumerate(secBC) if index>27 if char==0]
    b=[Bool('b%s' % (i)) for i in range(44)]
    f = [b[i] for i in RFs] 
    #print("f:", f)
    s=Solver()
    s.reset()
    #GGSO eqn
    for i in range(NumBVs):
        LHS=(-1)*GSO[1][i]
        for OscInd1 in range(27,44):
            for freq1 in range(2):
                for OscInd2 in range(27,44):
                    for freq2 in range(2):

                        OscBit=(-1)**(freq1*Basis[i][OscInd1])*(-1)**(freq2*Basis[i][OscInd2])
                        RRs=[]
                        for j in range(4):
                            if Basis[i][j]==1:
                                RRs.append(i)
                        if len(RRs)==0 and GSO[1][i]==1: # projects gravitino
                            NSUSY=0
                            printr("Model is non-supersymmetric")
                        else:
                        #print("RRs: ", RRs)
                            fRR=[b[k] for k in RRs] #ramond fermions to constrain
                            #print("fRR: ", fRR)
                            if LHS*OscBit==1:
                                s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==0) 
                            else:
                                s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==1) 
                    
                #print(s.check())
                            #if s.check() == sat: 
                                #print("Right-moving NS Oscillators: ")
                                #printr([OscInd1,OscInd2])
                                #NumOscStates+=1
                
                    
                            while s.check() == sat: 
                            #if s.check() == sat: 
                                m = s.model()
# =============================================================================
#                                 for v in vars:
#                                   print("%s = %5s" % (v, m.evaluate(v, model_completion = True))),
#                                 print
#                                 s.add(Or([p != v for p, v in [(v, m.evaluate(v, model_completion = True)) for v in vars]]))
# =============================================================================
                               
                                
                                stateB=[m.evaluate(f[i], model_completion = True) for i in range(len(RFs))]
                                #state=[ -1 if item.sexpr()=='true' else 0 for item in stateB]
                                
                                printr(stateB)
                                #countr+=1
                                #if countr%100==0:
                                #    print(countr)
                                s.add(Not(And([v() == m.evaluate(v, model_completion = True) for v in m]))) 
                            s.reset()
    #s.reset()
    if NSUSY!=0:
    #gravitino states:
        for i in range(NumBVs):
            LHS=(-1)*GSO[1][i]
            RRs=[]
            for j in range(4):
                if Basis[i][j]==1:
                    RRs.append(i)
            if len(RRs)!=0:
                #print("RRs: ", RRs)
                fRR=[b[k] for k in RRs] #ramond fermions to constrain
                s=Solver()
                s.reset()
                if GSO[1][i]==1:
                    s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==1) 
                else:
                    s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==0)
                #if s.check() == sat: 
            printr("Gravitino: dXbar^mu |S >")
                #NumOscStates+=1
            countr=0
            while s.check() == sat: 
            #if s.check() == sat: 
                
                m = s.model () 
                
                stateB=[m[f[i]] for i in range(len(RFs))]
                #state=[ -1 if item.sexpr()=='true' else 0 for item in stateB]
                printr(stateB)
                countr+=1
                #if countr%100==0:
                #    print(countr)
                s.add(Not(And([v() == m[v] for v in m]))) 
            s.reset()
            printr2("Number of gravitini", countr)

#def enhancements():
    
def ObservableMassless():
    for secBC in MSectorBC:
        
        Ind=MSectorBC.index(secBC)
        sec=MSector[Ind]
        secVacE = vacE(secBC)
        
        if secVacE[0]==4 and secVacE[1]==8:
            print("massless spinorial sector")
            
            print("sector:", sec)
            printr(sec)
            print("secBC:", secBC)
            
            #RFs=[] #Ramond fermions
            #for BC in secBC:
            RFs=[]
            RFs=[index for index, char in enumerate(secBC) if char == 1]
                #if BC==1:
                    #print(secBC.index(BC))
                    #RFs.append(secBC.index(BC))
            print(RFs)
            RFlen=len(RFs)
            b=[Bool('b%s' % (i)) for i in range(44)]
            f = [b[i] for i in RFs] 
            #print("f:", f)
            s=Solver()
            s.reset()
            #GGSO eqn
            
            LHS=(-1)**secBC[0]*SecGSO(sec[0],sec[1])
            #print(LHS)
            for secProj in secsBCs:#just need to loop through basis vecs
                #print("secProj:", secProj)
                #totalRR=DProd(secBC,secProj)
                RRs=[]
                for i in range(44):
                    if secProj[i]==1 and secBC[i]==1:
                        RRs.append(i)
                #print("RRs: ", RRs)
                fRR=[b[k] for k in RRs] #ramond fermions to constrain
                #print("fRR: ", fRR)
                if LHS==1:
                    s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==0) 
                else:
                    s.add(Sum([If(fRR[i],1,0) for i in range(len(fRR))])%2==1) 
            #print(s)
            print(s.check())
            if s.check() == sat: 
                m = s.model () 
                
                stateB=[m[f[i]] for i in range(len(RFs))]
                #state=[ -1 if item.sexpr()=='true' else 0 for item in stateB]
                printr(stateB)
                #countr+=1
                #if countr%100==0:
                #    print(countr)
                s.add(Not(And([v() == m[v] for v in m]))) 
