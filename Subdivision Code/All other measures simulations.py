# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:06:56 2020

@author: garry





"""

from simuPOP import *
import math
import simuPOP as sim
from simuPOP.utils import importPopulation, export
from simuPOP.sampling import drawRandomSample
from simuPOP.demography import *
import numpy as np
from numpy import nanmean
import csv


with open('MainTextData.csv', 'a', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(["Gen","D0rP1","D1rP1","D2rP1","D0cP1","D1cP1","D2cP1","D0rP2","D1rP2","D2rP2","D0cP2","D1cP2","D2cP2","D0rPT","D1rPT","D2rPT","D0cPT","D1cPT","D2cPT","D0aP1","D1aP1","D2aP1","D0aP2","D1aP2","D2aP2","D0aPT","D1aPT","D2aPT"])


def BillsRareAlleleIdea(pop):
    import statistics
################################### POP 1 #################################
    stat(pop, alleleFreq=ALL_AVAIL, subPops= [0])
    
    D0rP1 = []
    D1rP1 = []
    D2rP1 = []

    D0cP1 = []
    D1cP1 = []
    D2cP1 = []
    
    
    D0aP1 = []
    D1aP1 = []
    D2aP1 = []
        
    for i in range(LociNumber):
        D0cListP1 = []
        D1cListP1 = []
        D2cListP1 = []
        D0rListP1 = []
        D1rListP1 = []
        D2rListP1 = []
        
        D0aListP1 = []
        D1aListP1 = []
        D2aListP1 = []
        
        
        if pop.dvars().alleleFreq[i][0]*pop.dvars().alleleFreq[i][1]  == 0:
            yeah = "nah"
        else:
            if pop.dvars().alleleFreq[i][0] < 0.1 or pop.dvars().alleleFreq[i][1] < 0.1:
                for A in range(len(alProportions)):
                    D0rListP1.append(1)
                    D1rListP1.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2rListP1.append(pop.dvars().alleleFreq[i][A]**2)
                    
                    D0aListP1.append(1)
                    D1aListP1.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2aListP1.append(pop.dvars().alleleFreq[i][A]**2)


            else:
                for A in range(len(alProportions)):
                    D0cListP1.append(1)
                    D1cListP1.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2cListP1.append(pop.dvars().alleleFreq[i][A]**2)
                    
                    D0aListP1.append(1)
                    D1aListP1.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2aListP1.append(pop.dvars().alleleFreq[i][A]**2)
        
        if D0rListP1 == [] or D1rListP1 == [] or D2rListP1 == []:
            pass
        else:     
            D0rP1.append(1+ sum(D0rListP1)-1)
            D1rP1.append(math.exp(-sum(D1rListP1)))
            D2rP1.append(1/ (1 - (1 - sum(D2rListP1))))
        
        if D0cListP1 == [] or D1cListP1 == [] or D2cListP1 == []:
            pass
        else:                 
            D0cP1.append(1 + sum(D0cListP1)-1)
            D1cP1.append(math.exp(-sum(D1cListP1)))
            D2cP1.append(1/ (1 - (1 - sum(D2cListP1))))
            
        if D0aListP1 == [] or D1aListP1 == [] or D2aListP1 == []:
            pass
        else: 
            D0aP1.append(1 + sum(D0aListP1)-1)
            D1aP1.append(math.exp(-sum(D1aListP1)))
            D2aP1.append(1/ (1 - (1 - sum(D2aListP1))))
################################### POP 2 #################################
    stat(pop, alleleFreq=ALL_AVAIL, subPops= [1])
    
    D0rP2 = []
    D1rP2 = []
    D2rP2 = []

    D0cP2 = []
    D1cP2 = []
    D2cP2 = []

    D0aP2 = []
    D1aP2 = []
    D2aP2 = []  
      
    for i in range(LociNumber):
        D0cListP2 = []
        D1cListP2 = []
        D2cListP2 = []
        D0rListP2 = []
        D1rListP2 = []
        D2rListP2 = []
        
        D0aListP2 = []
        D1aListP2 = []
        D2aListP2 = []
        
        if pop.dvars().alleleFreq[i][0]*pop.dvars().alleleFreq[i][1]  == 0:
            pass
        else:
            if pop.dvars().alleleFreq[i][0] < 0.1 or pop.dvars().alleleFreq[i][1] < 0.1:
                for A in range(len(alProportions)):
                    D0rListP2.append(1)
                    D1rListP2.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2rListP2.append(pop.dvars().alleleFreq[i][A]**2)
                    
                    D0aListP2.append(1)
                    D1aListP2.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2aListP2.append(pop.dvars().alleleFreq[i][A]**2)
            else:
                for A in range(len(alProportions)):
                    D0cListP2.append(1)
                    D1cListP2.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2cListP2.append(pop.dvars().alleleFreq[i][A]**2)
                    
                    D0aListP2.append(1)
                    D1aListP2.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2aListP2.append(pop.dvars().alleleFreq[i][A]**2)

        if D0rListP2 == [] or D1rListP2 == [] or D2rListP2 == []:
            pass
        else:           
            D0rP2.append(1+ sum(D0rListP2)-1)
            D1rP2.append(math.exp(-sum(D1rListP2)))
            D2rP2.append(1/ (1 - (1 - sum(D2rListP2))))
            
        if D0cListP2 == [] or D1cListP2 == [] or D2cListP2 == []:
            pass
        else:           
            D0cP2.append(1 + sum(D0cListP2)-1)
            D1cP2.append(math.exp(-sum(D1cListP2)))
            D2cP2.append(1/ (1 - (1 - sum(D2cListP2))))
            
        if D0aListP2 == [] or D1aListP2 == [] or D2aListP2 == []:
            pass
        else:         
            D0aP2.append(1 + sum(D0aListP2)-1)
            D1aP2.append(math.exp(-sum(D1aListP2)))
            D2aP2.append(1/ (1 - (1 - sum(D2aListP2))))
                    
################################### TOTAL POPS #################################    
    stat(pop, alleleFreq=ALL_AVAIL, subPops= ALL_AVAIL)
    
    D0rPT = []
    D1rPT = []
    D2rPT = []

    D0cPT = []
    D1cPT = []
    D2cPT = []
    
    D0aPT = []
    D1aPT = []
    D2aPT = []
        
    for i in range(LociNumber):
        D0cListPT = []
        D1cListPT = []
        D2cListPT = []
        D0rListPT = []
        D1rListPT = []
        D2rListPT = []
        
        D0aListPT = []
        D1aListPT = []
        D2aListPT = []
        
        if pop.dvars().alleleFreq[i][0]*pop.dvars().alleleFreq[i][1]  == 0:
            pass
        else:
            if pop.dvars().alleleFreq[i][0] < 0.1 or pop.dvars().alleleFreq[i][1] < 0.1:
                for A in range(len(alProportions)):
                    D0rListPT.append(1)
                    D1rListPT.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2rListPT.append(pop.dvars().alleleFreq[i][A]**2)
                    
                                        
                    D0aListPT.append(1)
                    D1aListPT.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2aListPT.append(pop.dvars().alleleFreq[i][A]**2)

            else:
                for A in range(len(alProportions)):
                    D0cListPT.append(1)
                    D1cListPT.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2cListPT.append(pop.dvars().alleleFreq[i][A]**2)
                    
                    D0aListPT.append(1)
                    D1aListPT.append(pop.dvars().alleleFreq[i][A]* math.log(pop.dvars().alleleFreq[i][A]))
                    D2aListPT.append(pop.dvars().alleleFreq[i][A]**2)
                    
        if D0rListPT == [] or D1rListPT == [] or D2rListPT == []:
            pass
        else:        
            D0rPT.append(1+ sum(D0rListPT)-1)
            D1rPT.append(math.exp(-sum(D1rListPT)))
            D2rPT.append(1/ (1 - (1 - sum(D2rListPT))))
            
        if D0cListPT == [] or D1cListPT == [] or D2cListPT == []:
            pass
        else:                   
            D0cPT.append(1 + sum(D0cListPT)-1)
            D1cPT.append(math.exp(-sum(D1cListPT)))
            D2cPT.append(1/ (1 - (1 - sum(D2cListPT))))
       
        if D0aListPT == [] or D1aListPT == [] or D2aListPT == []:
            pass
        else:         
            D0aPT.append(1 + sum(D0aListPT)-1)
            D1aPT.append(math.exp(-sum(D1aListPT)))
            D2aPT.append(1/ (1 - (1 - sum(D2aListPT))))

    


    Output = [
     statistics.mean(D0rP1) if D0rP1 != [] else 'No Value'
    ,statistics.mean(D1rP1) if D1rP1 != [] else 'No Value'
    ,statistics.mean(D2rP1) if D2rP1 != [] else 'No Value'
    ,statistics.mean(D0cP1) if D0cP1 != [] else 'No Value'
    ,statistics.mean(D1cP1) if D1cP1 != [] else 'No Value'
    ,statistics.mean(D2cP1) if D2cP1 != [] else 'No Value'
    ,statistics.mean(D0rP2) if D0rP2 != [] else 'No Value'
    ,statistics.mean(D1rP2) if D1rP2 != [] else 'No Value'
    ,statistics.mean(D2rP2) if D2rP2 != [] else 'No Value'
    ,statistics.mean(D0cP2) if D0cP2 != [] else 'No Value'
    ,statistics.mean(D1cP2) if D1cP2 != [] else 'No Value'
    ,statistics.mean(D2cP2) if D2cP2 != [] else 'No Value'
    ,statistics.mean(D0rPT) if D0rPT != [] else 'No Value'
    ,statistics.mean(D1rPT) if D1rPT != [] else 'No Value'
    ,statistics.mean(D2rPT) if D2rPT != [] else 'No Value'
    ,statistics.mean(D0cPT) if D0cPT != [] else 'No Value'
    ,statistics.mean(D1cPT) if D1cPT != [] else 'No Value'
    ,statistics.mean(D2cPT) if D2cPT != [] else 'No Value'
    ,statistics.mean(D0aP1) if D0aP1 != [] else 'Fixed'
    ,statistics.mean(D1aP1) if D1aP1 != [] else 'Fixed'
    ,statistics.mean(D2aP1) if D2aP1 != [] else 'Fixed'
    ,statistics.mean(D0aP2) if D0aP2 != [] else 'Fixed'
    ,statistics.mean(D1aP2) if D1aP2 != [] else 'Fixed'
    ,statistics.mean(D2aP2) if D2aP2 != [] else 'Fixed'
    ,statistics.mean(D0aPT) if D0aPT != [] else 'Fixed'
    ,statistics.mean(D1aPT) if D1aPT != [] else 'Fixed'
    ,statistics.mean(D2aPT) if D2aPT != [] else 'Fixed']
        
    with open('MainTextData.csv', 'a', newline='') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow([pop.dvars().gen]+ Output)
        
    return True   
            

def calcDValuesHETFSTONLY(pop):
    #####'Calculate all D values,Alpha 0a,1a,2a and Beta 0b,1b,2b'
       
    #Alpha Measures and Total measures to be used later for Beta measures

    PerLociTotalPopHetDList = []
    listofFixedLoci = []
    stat(pop, alleleFreq=ALL_AVAIL, subPops= ALL_AVAIL)

    
    for i in range(LociNumber):
        listofalleleFreq = []
        for A in range(len(alProportions)):
            listofalleleFreq.append(pop.dvars().alleleFreq[i][A])

        if len(listofalleleFreq) == 1:
            listofFixedLoci.append(i)
        else:                
            ##Q = 2
            Hlist = []
            for f in listofalleleFreq:
                Hlist.append(f**2)                 
            Het  =  1 - sum(Hlist)
            PerLociTotalPopHetDList.append(Het)
            
        
        
        
    #same as before for Subpop 1       
    stat(pop, alleleFreq=ALL_AVAIL, subPops= [0])
    PerLociPop1HetDList = []

    for i in range(LociNumber):
        listofalleleFreq = []
        for A in range(len(alProportions)):
            listofalleleFreq.append(pop.dvars().alleleFreq[i][A])
            
        if i in listofFixedLoci:
            fixd = "yes"
        else:
            ##Q = 2
            Hlist = []
            for f in listofalleleFreq:
                Hlist.append(f**2)   
            Het  =  1 - sum(Hlist)
            PerLociPop1HetDList.append(Het)
    
    #same as before for Subpop 2    
    stat(pop, alleleFreq=ALL_AVAIL, subPops= [1])  
    PerLociPop2HetDList = []
    for i in range(LociNumber):
        listofalleleFreq = []
        
        for A in range(len(alProportions)):
            listofalleleFreq.append(pop.dvars().alleleFreq[i][A])
        
        if i in listofFixedLoci:
            fixd = "yes"
        else:        
            #this removes any 0 freq alleles, to then gives a len() of all lesf over alleles (Number of alleles)
            #placed here, it only does it if the total is nott fixed, s0 q0>1
            ##Q = 2
            Hlist = []
            for f in listofalleleFreq:
                Hlist.append(f**2)   
            Het  =  1 - sum(Hlist)
            PerLociPop2HetDList.append(Het)
    

    T = pop.subPopSizes()[1] +pop.subPopSizes()[0]
    P1 = pop.subPopSizes()[0]
    P2 = pop.subPopSizes()[1]
    AvgBetweenPopD2_AllLoci =[]
    for A in range(len(PerLociPop1HetDList)):
        AvgBetweenPopD2_AllLoci.append((PerLociPop1HetDList[A]*P1 + PerLociPop2HetDList[A]*P2)/(P1+P2))

    #AvgBetweenPopD2_AllLoci = np.average([PerLociPop1HetDList ,PerLociPop2HetDList], axis=0)
    D2a = np.average(PerLociTotalPopHetDList)
    D2as = np.average(AvgBetweenPopD2_AllLoci)

    stat(pop, structure= ALL_AVAIL, vars=['F_st', 'G_st'])
    
    Fst = pop.dvars().F_st
    
    if pop.dvars().gen < BurnInTime:
        PopState = 'BEFORE'
    elif pop.dvars().gen > BurnInTime:
        PopState = 'AFTER'
    else:
        PopState = 'DURING'
    
    
    with open("B.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow([pop.dvars().gen,D2a,D2as,Fst,PopState])
   
    #print(pop.dvars().gen, D0a, D1a, D2a, D0b, D1b, D2b)
    
    return True


#with open('B.csv', 'w', newline='') as myfile:
#     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
#     wr.writerow(['GAN','D2a','D2as','Fst','State'])

#############################################################################


with open('GangShitFinal.csv', 'a', newline='') as myfile:
     wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
     wr.writerow(['GAN','D0a','D0as','D0asWithFixed','D1a','D1as','D2a','D2as','D0b','D1b','D2b','Fst','Gst', 'ShanDiff','BillGandBratio', 'Fix1', 'Fix2', 'FixT','Jaccard','Sorrensen','BrayCurt', 'State', 'EvnAcr', 'EvnP1','EvnP2', 'PreRate', 'SplitRate', 'PopSize','AlleleProportions','RelativeGens'])


def PrintGen(pop):
    print(pop.dvars().gen)
    
    return True
    
def calcDValues(pop):   
    
    global pratio
    global qratio
    
    if pop.dvars().gen >= BurnInTime:           
        ExpectedratioListPerLoci = []
        stat(pop, alleleFreq=ALL_AVAIL, subPops= ALL_AVAIL)    
        for i in range(LociNumber):
                  
            listofalleleFreq = []
            for A in range(len(alProportions)):
                listofalleleFreq.append(pop.dvars().alleleFreq[i][A])
            
            p = listofalleleFreq[0]
            q = listofalleleFreq[1]
            
            pratio.append(p)
            qratio.append(q)
            
            if pop.dvars().gen < BurnInTime:
                t = 0
            else:
                t = pop.dvars().gen - BurnInTime
    
            
            if p*q == 0:
                x = 0
            else:    
                x = math.sqrt(  (1-(1-(1/(2*sum(PopSize)))**t)) / (4*p*q)   )
            
            ExpectedratioListPerLoci.append(x)

        BillGandBratio = sum(ExpectedratioListPerLoci)/len(ExpectedratioListPerLoci)
        
    elif pop.dvars().gen > BurnInTime:
        for i in range(len(pratio)):
            t = pop.dvars().gen - BurnInTime
            ExpectedratioListPerLoci = []

            if pratio[i]*qratio[i] == 0:
                x = 0
            else:    
                x = math.sqrt(  (1-(1-(1/(2*sum(PopSize)))**t)) / (4*pratio[i]*qratio[i])   )
                
            ExpectedratioListPerLoci.append(x)

        BillGandBratio = sum(ExpectedratioListPerLoci)/len(ExpectedratioListPerLoci)
    else:
        pratio = []
        qratio = []
        BillGandBratio = "Before Split"
    #####'Calculate all D values,Alpha 0a,1a,2a and Beta 0b,1b,2b'
       
    #Alpha Measures and Total measures to be used later for Beta measures
    PerLociTotalPopNumAllelesList = []
    PerLociTotalPopShanDList = []
    PerLociTotalPopHetDList = []
    MITotal = []
    stat(pop, alleleFreq=ALL_AVAIL, subPops= ALL_AVAIL)
    VarPerLociTotal = []
    listofFixedLoci = []
    PerLociTotalPopNumAllelesListWithFixed = []
    
    for i in range(LociNumber):
        
        PerLociTotalPopNumAllelesListWithFixed.append(len({x:y for x,y in pop.dvars().alleleFreq[i].items() if y!=0}))
        
        listofalleleFreq = []
        for A in range(len(alProportions)):
            listofalleleFreq.append(pop.dvars().alleleFreq[i][A])
            
        VariantsAtLociList = []
        for i in range(len(listofalleleFreq)):
            if listofalleleFreq[i] == 0:
                pass
            else:
                VariantsAtLociList.append(i)
        VarPerLociTotal.append(VariantsAtLociList)
            
        ##Q = 0
        #this removes any 0 freq alleles, to then gives a len() of all lesf over alleles (Number of alleles)
        
        ##Q = 1
        if len(listofalleleFreq) == 1:
            listofFixedLoci.append(i)
        else:                
            Slist = []
            for f in listofalleleFreq:
                if f == 0:
                    Slist.append(0)
                else:
                    Slist.append(-(f*math.log(f)))   
            PerLociTotalPopShanDList.append(math.exp(sum(Slist)))
            MITotal.append(sum(Slist))
            
            #this removes any 0 freq alleles, to then gives a len() of all lesf over alleles (Number of alleles)
            #placed here, it only does it if the total is nott fixed, s0 q0>1
            PerLociTotalPopNumAllelesList.append(len({x:y for x,y in pop.dvars().alleleFreq[i].items() if y!=0}))
            ##Q = 2
            Hlist = []
            for f in listofalleleFreq:
                Hlist.append(f**2) 

                
            Het  =  1 - sum(Hlist)
            PerLociTotalPopHetDList.append(1/(1-Het))
            
        
        
        
    #same as before for Subpop 1       
    stat(pop, alleleFreq=ALL_AVAIL, subPops= [0])
    PerLociPop1NumAllelesList = []
    PerLociPop1ShanDList = []
    PerLociPop1HetDList = []
    MIPop1 = []
    VarPerLociPop1 = []
    Allele1PropPop1 = []
    for i in range(LociNumber):
        if i in listofFixedLoci:
            fixd = "yes"
        else:
            listofalleleFreq = []
            Allele1PropPop1.append(pop.dvars().alleleFreq[i][0])
            for A in range(len(alProportions)):
                listofalleleFreq.append(pop.dvars().alleleFreq[i][A])
    
            
            VariantsAtLociList = []
            for i in range(len(listofalleleFreq)):
                if listofalleleFreq[i] == 0:
                    pass
                else:
                    VariantsAtLociList.append(i)                    
            VarPerLociPop1.append(VariantsAtLociList)
                    
                    
            ##Q = 0
            #this removes any 0 freq alleles, to then gives a len() of all lesf over alleles (Number of alleles)
            PerLociPop1NumAllelesList.append(len({x:y for x,y in pop.dvars().alleleFreq[i].items() if y!=0}))
            
            Slist = []
            for f in listofalleleFreq:
                if f == 0:
                    Slist.append(0)
                else:
                    Slist.append(-(f*math.log(f)))   
            PerLociPop1ShanDList.append(math.exp(sum(Slist)))
            MIPop1.append(sum(Slist))
            
            #this removes any 0 freq alleles, to then gives a len() of all lesf over alleles (Number of alleles)
            #placed here, it only does it if the total is nott fixed, s0 q0>1
            ##Q = 2
            Hlist = []
            for f in listofalleleFreq:
                Hlist.append(f**2)   
            Het  =  1 - sum(Hlist)
            PerLociPop1HetDList.append(1/(1-Het))
    
    #same as before for Subpop 2    
    stat(pop, alleleFreq=ALL_AVAIL, subPops= [1])
    PerLociPop2NumAllelesList = []
    PerLociPop2ShanDList = []
    PerLociPop2HetDList = []
    Allele1PropPop2 = []
    MIPop2 = []
    VarPerLociPop2 = []
    for i in range(LociNumber):
        if i in listofFixedLoci:
            fixd = "yes"
        else:
            listofalleleFreq = []
            Allele1PropPop2.append(pop.dvars().alleleFreq[i][0])
            for A in range(len(alProportions)):
                listofalleleFreq.append(pop.dvars().alleleFreq[i][A])
                
            VariantsAtLociList = []
            for i in range(len(listofalleleFreq)):
                if listofalleleFreq[i] == 0:
                    pass
                else:
                    VariantsAtLociList.append(i)
            VarPerLociPop2.append(VariantsAtLociList)
    
            
            ##Q = 0
            #this removes any 0 freq alleles, to then gives a len() of all lesf over alleles (Number of alleles)
            PerLociPop2NumAllelesList.append(len({x:y for x,y in pop.dvars().alleleFreq[i].items() if y!=0}))
            Slist = []
            for f in listofalleleFreq:
                if f == 0:
                    Slist.append(0)
                else:
                    Slist.append(-(f*math.log(f)))   
            PerLociPop2ShanDList.append(math.exp(sum(Slist)))
            MIPop2.append(sum(Slist))
            
            #this removes any 0 freq alleles, to then gives a len() of all lesf over alleles (Number of alleles)
            #placed here, it only does it if the total is nott fixed, s0 q0>1
            ##Q = 2
            Hlist = []
            for f in listofalleleFreq:
                Hlist.append(f**2)   
            Het  =  1 - sum(Hlist)
            PerLociPop2HetDList.append(1/(1-Het))



    T = pop.subPopSizes()[1]+pop.subPopSizes()[0]
    P1 = pop.subPopSizes()[0]
    P2 = pop.subPopSizes()[1]
    K = (-(P1/T)*math.log(P1/T))+(-(P2/T)*math.log(P2/T))
    

    
    AvgBetweenPopD0_AllLoci = []
    AvgBetweenPopD1_AllLoci = []
    AvgBetweenPopD2_AllLoci = []
    AvgBetweenPop_Evenness  = []
    pop1_eveness_allloci = []
    pop2_eveness_allloci = []
    
    for A in range(len(PerLociPop1ShanDList)):
        AvgBetweenPopD0_AllLoci.append((PerLociPop1NumAllelesList[A]*P1 + PerLociPop2NumAllelesList[A]*P2)/(P1+P2))
        AvgBetweenPopD1_AllLoci.append((PerLociPop1ShanDList[A]*P1 + PerLociPop2ShanDList[A]*P2)/(P1+P2))
        AvgBetweenPopD2_AllLoci.append((PerLociPop1HetDList[A]*P1 + PerLociPop2HetDList[A]*P2)/(P1+P2))
        #AvgBetweenPop_Evenness.append(((math.log(PerLociPop1ShanDList[A])/math.log(PerLociPop1NumAllelesList[A]))*P1 + (math.log(PerLociPop2ShanDList[A])/math.log(PerLociPop2NumAllelesList[A])*P2))/(P1+P2))
        if math.log(PerLociPop1NumAllelesList[A]) == 0:
            pop1_eveness_allloci.append(0)
        else:
            pop1_eveness_allloci.append(math.log(PerLociPop1ShanDList[A])/math.log(PerLociPop1NumAllelesList[A]))
        if math.log(PerLociPop2NumAllelesList[A]) == 0:
            pop2_eveness_allloci.append(0)
        else:
            pop2_eveness_allloci.append(math.log(PerLociPop2ShanDList[A])/math.log(PerLociPop2NumAllelesList[A]))
    
    AvgBetweenPopD2_AllLoci = []
    for A in range(len(PerLociPop1HetDList)):
        AvgBetweenPopD2_AllLoci.append((PerLociPop1HetDList[A]*P1 + PerLociPop2HetDList[A]*P2)/(P1+P2))
        
    AvgBetweenPopD0_AllLoci = np.average([PerLociPop1NumAllelesList ,PerLociPop2NumAllelesList], axis=0)
    #AvgBetweenPopD1_AllLoci = np.average([PerLociPop1ShanDList, PerLociPop2ShanDList], axis=0)
    #AvgBetweenPopD2_AllLoci = np.average([PerLociPop1HetDList ,PerLociPop2HetDList], axis=0)
    
    #Average over all loci to get a single Alpha D value for each pop
    D0aWithFix = np.average(PerLociTotalPopNumAllelesListWithFixed)
    D0a = np.average(PerLociTotalPopNumAllelesList)
    D1a = np.average(PerLociTotalPopShanDList)
    D2a = np.average(PerLociTotalPopHetDList)
    
    D0as = np.average(AvgBetweenPopD0_AllLoci)
    D1as = np.average(AvgBetweenPopD1_AllLoci)
    D2as = np.average(AvgBetweenPopD2_AllLoci)
    
    global listOfEvenessTotal
    
    listOfEvenessTotal = []
    for I in range(len(PerLociTotalPopNumAllelesList)):
        if math.log(PerLociTotalPopNumAllelesList[I]) == 0:
            pass
            #listOfEvenessTotal.append(0)
        else:
            listOfEvenessTotal.append(math.log(PerLociTotalPopShanDList[I])/math.log(PerLociTotalPopNumAllelesList[I]))
            #############################
            
    AvgBetweenPop_Evenness = np.average(listOfEvenessTotal)
    
    ###DOING WEIGHTED AVERAGES####
    EvenAvgAcr = np.average(AvgBetweenPop_Evenness)
    Evenpop1   = np.average(pop1_eveness_allloci)
    Evenpop2   = np.average(pop2_eveness_allloci)

     
    Pop1Fix = MIPop1.count(0)
    #print(Pop1Fix)
    Pop2Fix = MIPop2.count(0)
    #print(Pop2Fix)
    PopTFix = MITotal.count(0)
    #print(PopTFix)
    


    
    #divide Totalpop alpha values by averaged between pop alpha value to get Beta values
    #then average over loci to get single values
    
    D1as = np.average(AvgBetweenPopD1_AllLoci)
    
    D0b = np.nanmean(np.divide(PerLociTotalPopNumAllelesList, AvgBetweenPopD0_AllLoci))
    D1b = np.nanmean(np.divide(PerLociTotalPopShanDList, AvgBetweenPopD1_AllLoci))
    D2b = np.nanmean(np.divide(PerLociTotalPopHetDList, AvgBetweenPopD2_AllLoci))
    
    MIPerloci = []
    for I in range(len(AvgBetweenPopD1_AllLoci)):
        MIPerloci.append((PerLociTotalPopShanDList[I] - AvgBetweenPopD1_AllLoci[I])/K)
    MId = np.average(MIPerloci)
    #D0b = np.average(PerLociTotalPopNumAllelesList)/ np.average(AvgBetweenPopD0_AllLoci)
    #D1b = np.average(PerLociTotalPopShanDList)/ np.average(AvgBetweenPopD1_AllLoci)
    #D2b = np.average(PerLociTotalPopHetDList)/ np.average(AvgBetweenPopD2_AllLoci)
    
    JaccardPerloci = []
    SorrensenPerLoci = []
    for i in range(LociNumber):
        R = len(set(VarPerLociPop1[i]) & set(VarPerLociPop2[i]))
        Sy = len(VarPerLociTotal[i])
        Sab = len(VarPerLociPop1[i]) + len(VarPerLociPop2[i])
        
        JaccardPerloci.append(1-(R/Sy))
        SorrensenPerLoci.append(1-((2*R)/Sab))
    
    Jaccard   = np.average(JaccardPerloci)
    Sorrensen = np.average(SorrensenPerLoci)

    stat(pop, structure= ALL_AVAIL, vars=['F_st', 'G_st'])
    
    Fst = pop.dvars().F_st
    Gst = pop.dvars().G_st
    
    if pop.dvars().gen < BurnInTime:
        PopState = 'BEFORE'
    elif pop.dvars().gen >= BurnInTime:
        PopState = 'AFTER'
    elif pop.dvars().gen >= BurnInTime:
        PopState = 'DURING'

    Allele1PropPop1
    Allele1PropPop2
    BrayCurtAllLoci = []
    for i in range(LociNumber):
        if i in listofFixedLoci:
            fixd = "yes"
        else:
            Bray = Allele1PropPop1[i]-Allele1PropPop2[i]
            BrayCurtAllLoci.append(abs(Bray))
    
    BrayCurtAvg = np.average(BrayCurtAllLoci)
            
    with open("GangShitFinal.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow([pop.dvars().gen, D0a, D0as, D0aWithFix, D1a,D1as, D2a,D2as, D0b, D1b, D2b, Fst, Gst,MId,BillGandBratio,Pop1Fix,Pop2Fix,PopTFix,Jaccard,Sorrensen,BrayCurtAvg, PopState,EvenAvgAcr,Evenpop1,Evenpop2, MR, SR, PopSize, alProportions,(pop.dvars().gen-int(BurnInTime/4))])
   
    #print(pop.dvars().gen, D0a, D1a, D2a, D0b, D1b, D2b)

    return True


#Record Allele proportions

def RecordAlleleFreq(pop):
    stat(pop, alleleFreq=ALL_AVAIL, subPops= ALL_AVAIL)
    MaxAlleleFreqPerloci = []
    for i in range(LociNumber):
        ListofalleleFreq = []
        for A in range(len(alProportions)):
                ListofalleleFreq.append(pop.dvars().alleleFreq[i][A])
        MaxAlleleFreqPerloci.append(max(ListofalleleFreq))
    
    WrtieToFile = [pop.dvars().gen , np.average(MaxAlleleFreqPerloci)]
    WrtieToFile.extend(MaxAlleleFreqPerloci)        
    with open("Allele Frequencies over time Total.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(WrtieToFile)
    stat(pop, alleleFreq=ALL_AVAIL, subPops= [0])
    MaxAlleleFreqPerloci = []
    for i in range(LociNumber):
        ListofalleleFreq = []
        for A in range(len(alProportions)):
                ListofalleleFreq.append(pop.dvars().alleleFreq[i][A])
        MaxAlleleFreqPerloci.append(max(ListofalleleFreq))
    
    WrtieToFile = [pop.dvars().gen , np.average(MaxAlleleFreqPerloci)]
    WrtieToFile.extend(MaxAlleleFreqPerloci)        
    with open("Allele Frequencies over time pop1.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(WrtieToFile)
        
    stat(pop, alleleFreq=ALL_AVAIL, subPops= [1])
    MaxAlleleFreqPerloci = []
    for i in range(LociNumber):
        ListofalleleFreq = []
        for A in range(len(alProportions)):
                ListofalleleFreq.append(pop.dvars().alleleFreq[i][A])
        MaxAlleleFreqPerloci.append(max(ListofalleleFreq))
    
    WrtieToFile = [pop.dvars().gen , np.average(MaxAlleleFreqPerloci)]
    WrtieToFile.extend(MaxAlleleFreqPerloci)        
    with open("Allele Frequencies over time pop2.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(WrtieToFile)
    return True 



# -*- coding: utf-8 -*-

def measureFIS(pop):
    import csv
    import numpy
    from simuPOP.utils import  export
    
    export(pop, format='CSV', output='genepopfile.csv')
    file = 'genepopfile.csv'
    
    def func(a):
        if  a[0] == a[-1]:
            return 2
        else:
            return 1
    
     
    #EXPECTED HET PER LOCI
    a = numpy.genfromtxt(file, delimiter=',', skip_header=True)[:, 1:]
    HetE = []
    #ShanE = []
    HetO = []
    #ShanO = []
    #Sis = []
    Fis = []
    for i in range(LociNumber):
        res = numpy.unique(a[:, 1+(i*2):3+(i*2)],return_counts=True)
        y=0

        for i in range(len(res[1])):
            x = (res[1][i] / sum(res[1]))**2
            y = x + y
        HetE.append(1-y)

    #for i in range(LociNumber):
    #    res = numpy.unique(a[:, 1+(i*2):3+(i*2)],return_counts=True)
    #   y=0

    #    for i in range(len(res[1])):
    #        x = (res[1][i] / sum(res[1]))*math.log(res[1][i] / sum(res[1]))
    #        y = x + y
    #    ShanE.append(-y)
    

    #print("Expected Het per Loci", HetE)
    #OBSERVED HET PER LOCI
    for i in range(LociNumber):
        b = numpy.apply_along_axis( func , 1, a[:,1+(i*2):3+(i*2)])
        HetO.append(list(b).count(1)/len(list(b)))
        
    #for i in range(LociNumber):
    #    b = numpy.apply_along_axis( func , 1, a[:,1+(i*2):3+(i*2)])
    #    Prop1 = list(b).count(1)/len(list(b))
    #    Prop2 = list(b).count(2)/len(list(b))
    #    if Prop2 == 0 or Prop1 == 0:
    #        pass
    #    else:
    #        ShanO.append(-sum([(Prop1*math.log(Prop1)),(Prop2*math.log(Prop2))]))
    #print("Observed Het per Loci", HetO)
    
    #FIS per Loci
    for i in range(LociNumber):
        if HetE[i] == 0:
            Loci = "Fixed, so don't worry about it, hey" 
        else:
            Fis.append((HetE[i] - HetO[i])/HetE[i])
    #for i in range(LociNumber):
    #    if ShanE[i] == 0:
    #        Loci = "Fixed, so don't worry about it, hey" 
    #    else:
    #        Sis.append((ShanE[i] - ShanO[i])/ShanE[i])
        
    #print("Fis per Loci", Fis)
    if len(Fis) == 0:
        AvgFis = "Fixed"
    else:
        AvgFis = sum(Fis)/len(Fis)
    #print("Average Fis=", AvgFis)
    #if len(Sis) == 0:
    #    AvgSis = "Fixed"
    #else:
    #    AvgSis = sum(Sis)/len(Sis)
    #print("Average Fis=", AvgFis)
    
    with open("AvgFis.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow([AvgFis])

    return True

def measurePOOLEDFIS(pop):
    import csv
    import numpy
    from simuPOP.utils import  export
    
    export(pop, format='CSV', output='genepopfile.csv')
    file = 'genepopfile.csv'
    
    def func(a):
        if  a[0] == a[-1]:
            return 2
        else:
            return 1
    
     
    #EXPECTED HET PER LOCI
    a = numpy.genfromtxt(file, delimiter=',', skip_header=True)[:, 1:]
    HetE = []
    HetO = []
    
    res = numpy.unique(a[:, 1:],return_counts=True)
    y=0

    for i in range(len(res[1])):
        x = (res[1][i] / sum(res[1]))**2
        y = x + y
    
    HetE = 1-y
    

    #print("Expected Het per Loci", HetE)
    #OBSERVED HET PER LOCI
    HetOList = []
    for i in range(LociNumber):
        b = numpy.apply_along_axis( func , 1, a[:,1+(i*2):3+(i*2)])
        HetOList.append(b)
    flat_HetOList = []
    for sublist in HetOList:
        for item in sublist:
            flat_HetOList.append(item)
    HetO = (flat_HetOList.count(1)/len(flat_HetOList))
    #print("Observed Het per Loci", HetO)
    

    Fis = ((HetE - HetO)/HetE)
        
    with open("Fis.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow([Fis])
    
    return True

def GenoTypeShanPOOLED(pop):
    import csv
    import numpy
    from simuPOP.utils import  export
    import itertools
    
    export(pop, format='CSV', output='genepopfile.csv')
    file = 'genepopfile.csv'
    
    
    def func(b):
        if  b[0] > b[-1]:
            return str(100+b[0]) + str(100+b[-1])
        else:
            return str(100+b[-1]) + str(100+b[0])
    
     
    #EXPECTED HET PER LOCI
    a = numpy.genfromtxt(file, delimiter=',', skip_header=True)[:, 1:]
    Shan = []
    GenoPropList = []
    ShanNoFix = []
    MaxGenoEntropy = []
    EvenNessGeno = []
    
    for i in range(LociNumber):
        b = numpy.apply_along_axis( func , 1, a[:,1+(i*2):3+(i*2)])
        GenoPropList.append(list(b))
    
    #for i in range(LociNumber):
        
    flat_list = []
    for sublist in GenoPropList:
        for item in sublist:
            flat_list.append(item)        
    
    res = numpy.unique(flat_list,return_counts=True)
    
    y=0

    for i in range(len(res[1])):
        p = (res[1][i] / sum(res[1]))
        x = -(p*math.log(p))
        y = x + y
    Shan.append(y)
    
    #for i in range(LociNumber):
    U = numpy.unique(a[:, 1:])
    
    #C is the number of different possible genotypes given len(u) , aka number of variantes
    C = math.factorial(len(U)+2-1)/(math.factorial(2)*math.factorial(len(U)-1))
    MaxGenoEntropy.append(math.log(C))
        
    #Removed Fixed Loci
    #for i in range(LociNumber):
    if Shan == 0:
        Loci = "Fixed, so don't worry about it, hey" 
    else:
        ShanNoFix.append(Shan)
        EvenNessGeno.append(Shan[0]/MaxGenoEntropy[0])

    
    with open("GenoTypeShan.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow([ShanNoFix,EvenNessGeno])

    return True

def GenoTypeShan(pop):
    import csv
    import numpy
    from simuPOP.utils import  export
    import itertools
    
    export(pop, format='CSV', output='genepopfile.csv')
    file = 'genepopfile.csv'
    
    
    def func(b):
        if  b[0] > b[-1]:
            return str(100+b[0]) + str(100+b[-1])
        else:
            return str(100+b[-1]) + str(100+b[0])
    
     
    #EXPECTED HET PER LOCI
    a = numpy.genfromtxt(file, delimiter=',', skip_header=True)[:, 1:]
    Shan = []
    GenoPropList = []
    ShanNoFix = []
    MaxGenoEntropy = []
    EvenNessGeno = []
    
    for i in range(LociNumber):
        b = numpy.apply_along_axis( func , 1, a[:,1+(i*2):3+(i*2)])
        GenoPropList.append(list(b))
    
    for i in range(LociNumber):
        res = numpy.unique(GenoPropList[i],return_counts=True)
        
        y=0

        for i in range(len(res[1])):
            p = (res[1][i] / sum(res[1]))
            x = -(p*math.log(p))
            y = x + y
        Shan.append(y)
    
    for i in range(LociNumber):
        U = numpy.unique(a[:, 1+(i*2):3+(i*2)])
        
        #C is the number of different possible genotypes given len(u) , aka number of variantes
        C = math.factorial(len(U)+2-1)/(math.factorial(2)*math.factorial(len(U)-1))
        MaxGenoEntropy.append(math.log(C))
        
    #Removed Fixed Loci
    for i in range(LociNumber):
        if Shan[i] == 0:
            Loci = "Fixed, so don't worry about it, hey" 
        else:
            ShanNoFix.append(Shan[i])
            EvenNessGeno.append(Shan[i]/MaxGenoEntropy[i])
    if len(ShanNoFix) == 0:
        AvgShan = 0
    else:
        AvgShan = sum(ShanNoFix)/len(ShanNoFix)
    if len(EvenNessGeno) == 0:
        AvgEveness = 0
    else:
        AvgEveness = sum(EvenNessGeno)/len(EvenNessGeno)
    
    FilteredEveness =[]

    for i in range(len(listOfEvenessTotal)):
        if listOfEvenessTotal[i] > 0:
            FilteredEveness.append(0)
            #FilteredEveness.append((listOfEvenessTotal[i]-EvenNessGeno[i])/listOfEvenessTotal[i])
        else:
            pass
    if len(FilteredEveness) == 0:    
        AvgFilteredEveness = 0
    else:
        AvgFilteredEveness = sum(FilteredEveness)/len(FilteredEveness)
    with open("GenoTypeShan.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow([AvgShan,AvgEveness,AvgFilteredEveness])

    return True

def BillExpectedGSToverBrayCurt(pop):
    ExpectedratioListPerLoci = []
    stat(pop, alleleFreq=ALL_AVAIL, subPops= ALL_AVAIL)    
    for i in range(LociNumber):
              
        listofalleleFreq = []
        for A in range(len(alProportions)):
            listofalleleFreq.append(pop.dvars().alleleFreq[i][A])
        
        p = listofalleleFreq[0]
        q = listofalleleFreq[1]
        
        
        t = gen = BurnInTime
        (1-(1-(1/(2*sum(PopSize)))**t)) / (4*p*q)
        
        x = math.sqrt(  (1-(1-(1/(2*sum(PopSize)))**t)) / (4*p*q)   )
        
        ExpectedratioListPerLoci.append(x)
        
        
def WhoIsQ(pop):
    import statistics    
    Qn2Perloc = []
    Qn1Perloc = []
    Q0Perloc = []
    Q1Perloc = []
    Q2Perloc = []
    listofFixedLoci = []
    PerLociTotalPopNumAllelesListWithFixed = []

    stat(pop, alleleFreq=ALL_AVAIL, subPops= ALL_AVAIL)    
    for i in range(LociNumber):
              
        listofalleleFreq = []
        for A in range(len(alProportions)):
            listofalleleFreq.append(pop.dvars().alleleFreq[i][A])
            

        ##Q = 1
        if len(listofalleleFreq) == 1:
            listofFixedLoci.append(i)
        else:                
            Qn2list = []
            Qn1list = []
            Q0list = []
            Q1list = []
            Q2list = []
            
            for f in listofalleleFreq:
                if f == 0:
                    Qn2list.append(0) 
                    Qn1list.append(0) 
                    Q0list.append(0) 
                    Q1list.append(0)
                    Q2list.append(0) 
                else:
                    Qn2list.append(  (f**-2)) 
                    Qn1list.append(  (f**-1)) 
                    Q0list.append(  (f**0)) 
                    Q1list.append(  (f**0.99)) 
                    Q2list.append(  (f**2)) 
                    
            Qn2Perloc.append(sum(Qn2list)**(1/(1--2))) 
            Qn1Perloc.append(sum(Qn1list)**(1/(1--1))) 
            Q0Perloc.append(sum(Q0list)**(1/(1-0)))
            Q1Perloc.append(sum(Q1list)**(1/(1-0.99)))
            Q2Perloc.append(sum(Q2list)**(1/(1-2)))
    
    Qn2 = statistics.mean(Qn2Perloc)
    Qn1 = statistics.mean(Qn1Perloc)    
    Q0  = statistics.mean(Q0Perloc)
    Q1  = statistics.mean(Q1Perloc)
    Q2  = statistics.mean(Q2Perloc)
    
    with open("WhoIsQ.csv", "a", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow([Qn2, Qn1, Q0, Q1, Q2])
    
    return True

#LociNumber = 500
#SplitRate = [0,0.1,0.2,0.3,0.4,0.5]
#MigrationRate =  [0.5, 0.4 ,0.3, 0.2,0.1]
#POPSIZELIST       =  [[500,500],[160,160],[50,50],[16,16]]
#alProportionsLIST = [[0.5,0.5], [0.7,0.3], [0.2,0.2,0.2,0.2,0.2],[0.6,0.1,0.1,0.1,0.1]]


### 36 treatments [3 x 2 x 3 x 2]
LociNumber = 10000
SplitRate = [0]
MigrationRate =  [0.5]
POPSIZELIST       =  [[50,50]]
alProportionsLIST = [[0.5,0.5]]
Progress = 1



AllelPropHeader = ['Gen','Average Max Proportion']
for i in range(LociNumber):
    AllelPropHeader.append('Loci'+str(i))
    
with open("Allele Frequencies over time Total.csv", "w", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(AllelPropHeader)      
with open("Allele Frequencies over time pop1.csv", "w", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(AllelPropHeader)
with open("Allele Frequencies over time pop2.csv", "w", newline='') as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(AllelPropHeader)
        
#calculating Equlibrium Halflife
#HalfLife = math.log(0.5)  /(   math.log( ((1-MigrationRate)**2) * (1-(1/sum(PopSize))) )   )
for SR in SplitRate:
    for MR in MigrationRate:
        if SR > MR:
            pass
        else:
            for PopSize in POPSIZELIST:
                for alProportions in alProportionsLIST:     
                    #ProgressBar
                    Progress = Progress + 1 
                    print(str("Simulation is " + str((Progress/(len(SplitRate)*len(MigrationRate)*len(POPSIZELIST)*len(alProportionsLIST))) *100) + "% done."))
                        
                    for G in range(100):
                    
                        #print(MR)
                        #Calculate gens till fixation. I'll be using a quater of this time till the split, and looking at 50 gen before an after. 
                        BurnInTime = int(-((4 * sum(PopSize) * alProportions[0] * math.log(alProportions[0]))/(1-alProportions[0])))
                        #print(BurnInTime)
                        pop = Population(size=PopSize, loci=LociNumber,infoFields=['migrate_to'])
                        pop.evolve(
                                preOps=[
                                    Migrator(rate=[[MR,MR],[MR,MR]], end=(int(BurnInTime/4))),
                                    Migrator(rate=[[SR,SR],[SR,SR]], begin=(int(BurnInTime/4))),
                                    #Migrator(rate=[[0,MR],[MR,0]], begin=(int(BurnInTime)+20)),
                                ],
                                #setting up the initial perameters, even sex ratio and allele frequencies
                                initOps = [
                                        InitSex(maleProp=0.5),
                                        InitGenotype(freq=alProportions, loci = ALL_AVAIL),
                        
                                        ],
                                  
                                matingScheme= RandomMating(sexMode=(GLOBAL_SEQUENCE_OF_SEX, MALE, FEMALE)),
                                postOps = [
                                        #PyOperator(func=PrintGen, step = 10),
                                        #PyOperator(func=BillsRareAlleleIdea, begin= (int(BurnInTime)-50), end = (int(BurnInTime)+50), step = 1),
                                        #PyOperator(func=WhoIsQ, begin= (int(BurnInTime)-50), end = (int(BurnInTime)+50), step = 1),
                                        #PyOperator(func=measureFIS, begin= (int(BurnInTime)-50), end = (int(BurnInTime)+50), step = 1),
                                        #PyOperator(func=GenoTypeShan, begin= (int(BurnInTime)-50), end = (int(BurnInTime)+50), step = 1),
                                        PyOperator(func=calcDValues, begin= (int(BurnInTime/4)-10), end = (int(BurnInTime/4)+10), step = 1),
                                        
                                        
                                        #PyOperator(func=GenoTypeShan, begin=250, end = 350, step = 1),
                                        #PyOperator(func=measureFIS, begin=250, end = 350, step = 1),
                                        #PyOperator(func=calcDValues, step=1, begin=(int(BurnInTime/2)), end=(int(BurnInTime/2)+int(BurnInTime/10))),
                                        #PyOperator(func=calcDValues, step=20, begin=(int(BurnInTime/2)+int(BurnInTime/10))),
                                        
                                        #PyOperator(func=RecordAlleleFreq, step=1)
                                        
                                                ],
                                gen=(int(BurnInTime)+12)
                                )
                        print(G)
                        
        
