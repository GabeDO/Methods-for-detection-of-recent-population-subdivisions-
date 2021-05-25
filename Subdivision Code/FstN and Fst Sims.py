# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 12:02:45 2021

@author: z5035086
"""
        
        


from simuPOP import *
import csv
import math

with open('FstN and Fst Sim Results.csv', 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["Generations","Fst", "Migration", "Popsize"])
        #writer.writerow(["Generations","Popsize", "Fst", "Gst", "Migration", "WitEst", "GabeEst", "MaxFst", "GabePrediction"])
        
PreRate = [0.5]
SplitRate = [0, 0.001]
PopsizeList = [[25,25],[50,50],[500,500]]
AlleleFreqList = [[0.5,0.5]]
NumLoci = 10000
import time
def makedata(pop):
    import csv
    import statistics
    import math
    
    stat(pop, structure= ALL_AVAIL,alleleFreq= ALL_AVAIL, vars=['F_st', 'G_st', 'F_is'])
    
    Gen = pop.dvars().gen - X
    Fst = pop.dvars().G_st
   
    with open('FistNTestFstDisparityLONGDATApop1000.csv', 'a', newline='') as file:
            writer = csv.writer(file)
            writer.writerow([Gen, Fst, SR, sum(Popsize)])
            #writer.writerow([Gen, sum(Popsize) , Fst, Gst, SR, WitLock, Gabe, MaxFst, GabePrediction])
            
    return True
    

for i in range(1000):
    print(i)
    for AlleleFreq in AlleleFreqList:
        for Popsize in PopsizeList:
            for SR in SplitRate:
                for PR in PreRate:
                    X = int( (-((4 * sum(Popsize) * AlleleFreq[0] * math.log(AlleleFreq[0]))/(1-AlleleFreq[0])))/4)
                    pop = Population(size=Popsize, loci=NumLoci,infoFields=['migrate_to','parent_idx'])
                    pop.evolve(
                        preOps=[
                            Migrator(rate=[[PR,PR],[PR,PR]], end=X),
                            Migrator(rate=[[SR,SR],[SR,SR]], begin=X),        
                            ],
                            
                        initOps = [
                                    InitSex(sex=[MALE, FEMALE], subPops = ALL_AVAIL),
                                    InitGenotype(freq=AlleleFreq, loci = ALL_AVAIL),
                    
                                    ],
                              
                            matingScheme= RandomMating(sexMode=(GLOBAL_SEQUENCE_OF_SEX, MALE, FEMALE)),
                            postOps = [
                                    PyOperator(func=makedata, begin= X-10, step = 1, end = X+10),
                                    PyOperator(func=makedata, begin= X+10, step = 5)
                                    
                                            ],
                            gen=(X+16)
                            )
                    
                