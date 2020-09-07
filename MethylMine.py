# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 18:20:56 2020

@author: s4142554
methylMine8
updated from methylMine and methylMine3 using my - only compares q1s
updated from methylMine and methylMine7 from mDist_024
"""



from operator import itemgetter, is_not
from collections import OrderedDict
import datetime
import time
import sys
from math import log10, sqrt
import re
from functools import partial 
from time import gmtime, strftime
import os
import gc
import shutil
import numpy as np
import pandas as pd
import random
import seaborn as sns
import itertools
from scipy.stats import mannwhitneyu
from statistics import pvariance as pvar

from bioinfokit import analys, visuz


###############################################################################
script = os.path.basename(__file__).split('.',1)[0]
script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
###############################################################################

def filtMatrix(probesList, cancer, tissue, Matrix, num=None):

    
    df = []
    header = []
    
    with open(Matrix, 'r') as file:
        for i, line in enumerate(file):
            line = line.split()
            
            probe = line[0]
            
            if probe == "ProbeID":
                header = ["ProbeID"]
                if num:
                    header.extend(["_".join([cancer, tissue, str(n+1)]) for n in range(num)])
                else:
                    header.extend(["_".join([cancer, tissue, str(n+1)]) for n in range(len(line[1:]))])
                # df.append(rline)
                
            elif probe in probesList:
                # replace NA with mean
                rline = [probe]
                mBeta = np.mean([float(n) for n in line[1:] if n != "NA"])
                beta = [float(n) if n != "NA" else mBeta for n in line[1:]]
                if num:
                    if len(beta) >= num:
                        beta = random.sample(beta, num)
                    else:
                        beta.extend([mBeta for n in range(num-len(beta))])
                    
                rline.extend(beta)
                df.append(rline)
                
    
    df = pd.DataFrame(df, columns=header)
    df = df.set_index("ProbeID")

    
    # print (df.head())
    return df

def heatmap(data, resultFigure):
    
    try:
        g = sns.clustermap(data)
        g.savefig(resultFigure)
    except  MemoryError as m:
        print (m) #break messege
    return g


#matrix1 - ptumor, matrix2 - tnorm
#create dictionary of quantiles of all probes in beta matrix file- reports q1, median(q2), q3
#Matrix1 - pTumor, Matrix2 = tNorm, n1 - pTumor, tNorm - tNorm
def QuantileCalc(InputMatrix1, InputMatrix2, n1=None, n2=None):
    
    quantileDict = OrderedDict()
    # ProbesList = []
    
    with open(InputMatrix1, 'r') as file1, open(InputMatrix2, 'r') as file2:
        for i, (pline, tline) in enumerate(zip(file1, file2)):
            pline = pline.split() #primary tumor
            tline = tline.split() #tissue normal
            
            # process the first i lines
            # if i == 30000:
            #     break


            # if the probes match
            if pline[0] == tline[0] and pline[0] != "ProbeID":
                probe = tline[0]
                
                if i % 10000 == 0:
                    print (i, probe)                
                
                #if all in ptumor line == NA
                if all(n == "NA" for n in pline[1:]):
                    ptumor_q1 = 0.5
                    ptumor_q2 = 0.5
                    ptumor_q3 = 0.5
                else:
                    #convert all beta to floats
                    ptumor_beta = [float(n) for n in pline[1:] if n != "NA"]
                    #get subset if n1 and n2 are specified
                    if n1:
                        if len(ptumor_beta) >= n1:
                            ptumor_beta = random.sample(ptumor_beta, n1)
                        else:
                            mbeta = np.mean(ptumor_beta)
                            ptumor_beta.extend([mbeta for n in range(n1-len(ptumor_beta))])                    
                    
                    ptumor_q1 = np.percentile(ptumor_beta, 25)
                    ptumor_q2 = np.percentile(ptumor_beta, 50)
                    ptumor_q3 = np.percentile(ptumor_beta, 75)
                
                #if  all in tnorm line == NA
                if all(n == "NA" for n in tline[1:]):
                    tnorm_q1 = 0.5
                    tnorm_q2 = 0.5
                    tnorm_q3 = 0.5
                else:
                    #convert all beta to floats
                    tnorm_beta = [float(n) for n in tline[1:] if n != "NA"] 
                    #get subset if n1 and n2 are specified
                    if n2:
                        if len(tnorm_beta) >= n2:
                            tnorm_beta = random.sample(tnorm_beta, n2)
                        else:
                            mbeta = np.mean(tnorm_beta)
                            tnorm_beta.extend([mbeta for n in range(n1-len(tnorm_beta))])
                        
                    tnorm_q1 = np.percentile(tnorm_beta, 25)
                    tnorm_q2 = np.percentile(tnorm_beta, 50)
                    tnorm_q3 = np.percentile(tnorm_beta, 75)                    
                quantileDict[probe] = [ptumor_q1, ptumor_q2, ptumor_q3, tnorm_q1, tnorm_q2, tnorm_q3]
    return quantileDict


#generate Dictionary of Quantiles for each cancer
#folder suffix - cancer_"merged-matrix"
#file suffix - "merged-matrix.txt"
def QuantileLibrary(cancerList, matrixFolder, folderSuffix, file_suffix, n1=None, n2=None):
    
    QuantileDict = OrderedDict()
    Tissues = ["PrimaryTumor", "TissueNormal"]
    
    
    #folder where all matrices are kept
    for i, cancer in enumerate(cancerList):
        cancerMatrixFolder = os.path.join(script, matrixFolder, "_".join([cancer, folderSuffix])) #folder suffix - cancer_merged-matrix
        
        pTumor_Matrix = os.path.join(script, cancerMatrixFolder, "_".join([cancer, Tissues[0], file_suffix]))
        tNorm_Matrix = os.path.join(script, cancerMatrixFolder, "_".join([cancer, Tissues[1], file_suffix]))
        
        print ("Processing...", i+1, 'of', len(cancerList), 'cancers:', cancer)
        
        #check if both files exist before continuing
        if os.path.isfile(pTumor_Matrix) and os.path.isfile(tNorm_Matrix):
            #get quartile (25 and 75) info
            qDict = QuantileCalc(pTumor_Matrix, tNorm_Matrix, n1, n2)
            QuantileDict[cancer] = qDict
            print ('total length:', len(qDict))
            
        
    return QuantileDict


###############################################################################
#main functions here
class DiffAnalysis(object):
    '''
    main differential probe class function
    processes up to 12 cancers simultaneously
    '''
    def __init__(self, cancers,  matrixFolder, folder_suffix, file_suffix, resultFolder, lower_thres_beta, upper_thres_beta, n1=None, n2=None):
        self.cancers = cancers #list of all cancers entered
        self.quantileDict = QuantileLibrary(cancers, matrixFolder, folder_suffix, file_suffix, n1, n2) #[] of all cancers listed
        self.lower_thres_beta = lower_thres_beta #lower beta cutoff e.g. 0.3
        self.upper_thres_beta = upper_thres_beta #upper beta cutoff e.g. 0.7
        self.folder_suffix = folder_suffix
        
        self.matrixFiles = ""
        
        #retrive name of all cancers for calling quantiles
        self.cancer1 = cancers[0]
        self.cancer2 = cancers[1] if len(cancers) >= 2 else ''
        self.cancer3 = cancers[2] if len(cancers) >= 3 else ''    
        self.cancer4 = cancers[3] if len(cancers) >= 4 else ''
        self.cancer5 = cancers[4] if len(cancers) >= 5 else ''
        self.cancer6 = cancers[5] if len(cancers) >= 6 else ''
        self.cancer7 = cancers[6] if len(cancers) >= 7 else ''
        self.cancer8 = cancers[7] if len(cancers) >= 8 else ''
        self.cancer9 = cancers[8] if len(cancers) >= 9 else ''
        self.cancer10 = cancers[9] if len(cancers) >= 10 else ''
        self.cancer11 = cancers[10] if len(cancers) >= 11 else ''
        self.cancer12 = cancers[11] if len(cancers) >= 12 else ''

        #retrive quantiles
        self.cancer1_quantiles = self.quantileDict[self.cancer1] 
        self.cancer2_quantiles = self.quantileDict[self.cancer2] if len(cancers) >= 2 else [{}, {}, {}, {}]
        self.cancer3_quantiles = self.quantileDict[self.cancer3] if len(cancers) >= 3 else [{}, {}, {}, {}]
        self.cancer4_quantiles = self.quantileDict[self.cancer4] if len(cancers) >= 4 else [{}, {}, {}, {}]  
        self.cancer5_quantiles = self.quantileDict[self.cancer5] if len(cancers) >= 5 else [{}, {}, {}, {}]
        self.cancer6_quantiles = self.quantileDict[self.cancer6] if len(cancers) >= 6 else [{}, {}, {}, {}]
        self.cancer7_quantiles = self.quantileDict[self.cancer7] if len(cancers) >= 7 else [{}, {}, {}, {}]
        self.cancer8_quantiles = self.quantileDict[self.cancer8] if len(cancers) >= 8 else [{}, {}, {}, {}]
        self.cancer9_queantiles = self.quantileDict[self.cancer9] if len(cancers) >= 9 else [{}, {}, {}, {}]
        self.cancer10_quantiles = self.quantileDict[self.cancer10] if len(cancers) >= 10 else [{}, {}, {}, {}]
        self.cancer11_quantiles = self.quantileDict[self.cancer11] if len(cancers) >= 11 else [{}, {}, {}, {}]
        self.cancer12_quantiles = self.quantileDict[self.cancer12] if len(cancers) >= 12 else [{}, {}, {}, {}]


        #probes list obtained from cancer1 quantiles
        self.probesList = list(self.cancer1_quantiles.keys())
        
        #emptyList for storing all results
        self.resultsList = []

    def compare(self):
        
        def cancerCombo(num):
            
            compareList = []
            
            cancerList = list(itertools.combinations(self.cancers, num)) #cancerList combos
            
            #for num combination of cancers, add remaining cancers to compare
            #e.g for 3 cancer - [[cancer1, cancer2], cancer3]
            for combo in cancerList:
                comboList = [line for line in self.cancers if line not in combo]
                compareList.append([list(combo), comboList])
                
            return compareList

        #checks for hypermethylation of single or multiple cancers
        def hyperCheck(tnorm_q3s, ptumor_q1s, quantileList):
            
            #quantileList - [list of q1s]
            
            result = False
            #tnorm_13s and ptumor q1s can be either a float or a list
            if isinstance(tnorm_q3s, float) and isinstance(ptumor_q1s, float):
                if ptumor_q1s >= self.upper_thres_beta and tnorm_q3s <= self.lower_thres_beta:
                    if all(n <= self.lower_thres_beta for n in quantileList):
                        result = True
            
            elif isinstance(tnorm_q3s, list) and isinstance(ptumor_q1s, list):
                if all(n >= self.upper_thres_beta for n in ptumor_q1s) and all(n <= self.lower_thres_beta for n in tnorm_q3s):
                    if all(n <= self.lower_thres_beta for n in quantileList):
                        result = True
            
            return result

        def hypoCheck(tnorm_q1s, ptumor_q3s, quantileList):
            
            #quantileList - [list of q1s]
            
            result = False
            
            #tnorm_q1s and ptumor_q3s can be a list of quantiles or a single quantile
            if isinstance(tnorm_q1s, float) and isinstance(ptumor_q3s, float):
                if ptumor_q3s <= self.lower_thres_beta and tnorm_q1s >= self.upper_thres_beta:
                    if all(n >= self.upper_thres_beta for n in quantileList):
                        result = True
            
            elif isinstance(tnorm_q1s, list) and isinstance(ptumor_q3s, list):
                if all(n <= self.lower_thres_beta for n in ptumor_q3s) and all(n >= self.upper_thres_beta for n in tnorm_q1s):
                    if all(n >= self.upper_thres_beta for n in quantileList):
                        result = True

            return result 

        #checks for hypermethylation of tissue (both tnorm and ptumor)
        def tissue_hyperCheck(sampleQuantiles, cancerQunatiles):

            #sampleQuantiles = [ptumor1_q1, tnorm1_q1] 
            #cancerQuanties = [ptumor2_q1, tnorm2_q1, etc]
            
            result = False

            if all(n >= self.upper_thres_beta for n in sampleQuantiles):  #[ptumor_q1, tnorm_q1] 
                if all(n <= self.lower_thres_beta for n in cancerQunatiles):  #[ptumor2_q3, tnorm2_q3, etc]
                    result = True            
            
            return result
        
        
        #checks for hypmethylation of tissue (both tnorm and ptumor)
        def tissue_hypoCheck(sampleQuantiles, cancerQunatiles):
        
            #sampleQuantiles = [ptumor1_q1, tnorm1_q1] 
            #cancerQuanties = [ptumor2_q1, tnorm2_q1, etc]
            
            result = False
        
            if all(n <= self.lower_thres_beta for n  in sampleQuantiles): #[ptumor_q3, tnorm_q3] 
                if all(n >= self.upper_thres_beta for n in cancerQunatiles): #[ptumor2_q1, tnorm2_q1, etc]
                    result = True
                    
            return result

        #check for hypermethylation across all cancers 
        def pancancer_hyperCheck(ptumor_q1s, tnorm_q3s):
            
            #ptumor_q1s - list of q1s from all cancers
            #tnorm_q1s - list of q1s from all cancers
            #in this case, ptumor_q1s and tnorm_q3s are all lists
            
            result = False
            
            if all(n >= self.upper_thres_beta for n in ptumor_q1s) and all(n <= self.lower_thres_beta for n in tnorm_q3s):
                result = True

            return result
        
        #check for hypomethylation across all cancers
        def pancancer_hypoCheck(ptumor_q3s, tnorm_q1s):
            
            #ptumor_q1s - list of q1s from all cancers
            #tnorm_q1s - list of q1s from all cancers
            #in this case, ptumor_q1s and tnorm_q1s are all lists
            
            result = False
            
            if all(n <= self.lower_thres_beta for n in ptumor_q3s) and all(n >= self.upper_thres_beta for n in tnorm_q1s):
                result = True    
                
            return result 
        
        #derived from QuantCompare 
        def compare_probes():
            
            '''
            cancer(s) specific
            + -> hypermethylated
            - -> hypomethylated
            
            tissue(s) specific
            #+ -> hypermethylated
            #- -> hypomethylated
            
            
            '''
            
            #tally
            num1 = 0 # tumor-specific hypermethylation
            num2 = 0 # tumor-specific hypomethylation
            num3 = 0 # tissue-specific hypermethylation
            num4 = 0 # tissue-specific hypomethylation
            num5 = 0 # pan-cancer hypermethylation
            num6 = 0 # pan-cancer hypomethylation
            
            #go through each probe in probes list
            for i, probe in enumerate(self.probesList):

                if i % 10000 == 0:
                    print (i, probe)      
                
                result = []

                ptumor_q1 = {} # dict for storing ptumor q1
                ptumor_q2 = {} # dict for storing ptumor q2 - median
                ptumor_q3 = {} # dict for storing ptumor q3
                tnorm_q1 ={} # dict for storing tnorm q1
                tnorm_q2 = {} # dict for storing tnorm q2 - median
                tnorm_q3 = {} # dict for storing tnorm q3
                
                #pull ptumor and tnorm quartiles into dictionary - [ptumor_q1, ptumor_q2, ptumor_q3, tnorm_q1, tnorm_q2, tnorm_q3
                # q1_pTumorDict[self.cancer1], q3_pTumorDict[self.cancer1], q1_tNormDict[self.cancer1], q3_tNormDict[self.cancer1] = self.cancer1_quantiles[probe]
                
                ptumor_q1[self.cancer1],ptumor_q2[self.cancer1], ptumor_q3[self.cancer1],tnorm_q2[self.cancer1], tnorm_q1[self.cancer1], tnorm_q3[self.cancer1] = self.cancer1_quantiles[probe]
                ptumor_q1[self.cancer2],ptumor_q2[self.cancer2], ptumor_q3[self.cancer2],tnorm_q2[self.cancer2], tnorm_q1[self.cancer2], tnorm_q3[self.cancer2] = self.cancer2_quantiles[probe] if len(self.cancers) >= 2 else [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
                ptumor_q1[self.cancer3],ptumor_q2[self.cancer3], ptumor_q3[self.cancer3],tnorm_q2[self.cancer3], tnorm_q1[self.cancer3], tnorm_q3[self.cancer3] = self.cancer3_quantiles[probe] if len(self.cancers) >= 3 else [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]     
                ptumor_q1[self.cancer4],ptumor_q2[self.cancer4], ptumor_q3[self.cancer4],tnorm_q2[self.cancer4], tnorm_q1[self.cancer4], tnorm_q3[self.cancer4] = self.cancer4_quantiles[probe] if len(self.cancers) >= 4 else [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]     
                ptumor_q1[self.cancer5],ptumor_q2[self.cancer5], ptumor_q3[self.cancer5],tnorm_q2[self.cancer5], tnorm_q1[self.cancer5], tnorm_q3[self.cancer5] = self.cancer5_quantiles[probe] if len(self.cancers) >= 5 else [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]     
                ptumor_q1[self.cancer6],ptumor_q2[self.cancer6], ptumor_q3[self.cancer6],tnorm_q2[self.cancer6], tnorm_q1[self.cancer6], tnorm_q3[self.cancer6] = self.cancer6_quantiles[probe] if len(self.cancers) >= 6 else [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]     
                ptumor_q1[self.cancer7],ptumor_q2[self.cancer7], ptumor_q3[self.cancer7],tnorm_q2[self.cancer7], tnorm_q1[self.cancer7], tnorm_q3[self.cancer7] = self.cancer7_quantiles[probe] if len(self.cancers) >= 7 else [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
                ptumor_q1[self.cancer8],ptumor_q2[self.cancer8], ptumor_q3[self.cancer8],tnorm_q2[self.cancer8], tnorm_q1[self.cancer8], tnorm_q3[self.cancer8] = self.cancer8_quantiles[probe] if len(self.cancers) >= 8 else [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
                ptumor_q1[self.cancer9],ptumor_q2[self.cancer9], ptumor_q3[self.cancer9],tnorm_q2[self.cancer9], tnorm_q1[self.cancer9], tnorm_q3[self.cancer9] = self.cancer8_quantiles[probe] if len(self.cancers) >= 9 else [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
                ptumor_q1[self.cancer10],ptumor_q2[self.cancer10], ptumor_q3[self.cancer10],tnorm_q2[self.cancer10], tnorm_q1[self.cancer10], tnorm_q3[self.cancer10] = self.cancer10_quantiles[probe] if len(self.cancers) >= 10 else [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
                ptumor_q1[self.cancer11],ptumor_q2[self.cancer11], ptumor_q3[self.cancer11],tnorm_q2[self.cancer11], tnorm_q1[self.cancer11], tnorm_q3[self.cancer11] = self.cancer11_quantiles[probe] if len(self.cancers) >= 11 else [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
                ptumor_q1[self.cancer12],ptumor_q2[self.cancer12], ptumor_q3[self.cancer12],tnorm_q2[self.cancer12], tnorm_q1[self.cancer12], tnorm_q3[self.cancer12] = self.cancer12_quantiles[probe] if len(self.cancers) >= 12 else [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
                        
                #read quantiles into q1_pTumorDict[cancer] = ptumor_q1, q1_tNormDict[cancer] = tnorm_q1 etc
                #if only one cancer in list
                if len(self.quantileDict) == 1:
                    pass

                #if there is two cancers in the list
                elif len(self.quantileDict) == 2:
                    #list of quantiles - for cancer hyper and hypo check 
                    pass

                elif len(self.quantileDict) == 3: #ptumor_q3[self.cancer2], tnorm_q1[self.cancer2]
                    #list of quantiles to compare single cancers
                    pass
                
                elif len(self.quantileDict) == 4:
                    #list of quantiles to compare single cancers
                    hyper_quantiles1 = [ptumor_q3[self.cancer2], tnorm_q3[self.cancer2], ptumor_q3[self.cancer3], tnorm_q3[self.cancer3], ptumor_q3[self.cancer4], tnorm_q3[self.cancer4]] #cancer 1 -> cancers 2 and 3, 4
                    hyper_quantiles2 = [ptumor_q3[self.cancer1], tnorm_q3[self.cancer1], ptumor_q3[self.cancer3], tnorm_q3[self.cancer3], ptumor_q3[self.cancer4], tnorm_q3[self.cancer4]] #cancer 2 -> cancers 1 and 3, 4
                    hyper_quantiles3 = [ptumor_q3[self.cancer1], tnorm_q3[self.cancer1], ptumor_q3[self.cancer2], tnorm_q3[self.cancer2], ptumor_q3[self.cancer4], tnorm_q3[self.cancer4]] #cancer 3 -> cancers 1 and 2, 4
                    hyper_quantiles4 = [ptumor_q3[self.cancer1], tnorm_q3[self.cancer1], ptumor_q3[self.cancer2], tnorm_q3[self.cancer2], ptumor_q3[self.cancer3], tnorm_q3[self.cancer3]] #cancer 4 -> cancers 1 and 2, 3
                    
                    hypo_quantiles1 = [ptumor_q1[self.cancer2], tnorm_q1[self.cancer2], ptumor_q1[self.cancer3], tnorm_q1[self.cancer3], ptumor_q1[self.cancer4], tnorm_q1[self.cancer4]] #cancer 1 -> cancers 2 and 3, 4
                    hypo_quantiles2 = [ptumor_q1[self.cancer1], tnorm_q1[self.cancer1], ptumor_q1[self.cancer3], tnorm_q1[self.cancer3], ptumor_q1[self.cancer4], tnorm_q1[self.cancer4]] #cancer 2 -> cancers 1 and 3, 4
                    hypo_quantiles3 = [ptumor_q1[self.cancer1], tnorm_q1[self.cancer1], ptumor_q1[self.cancer2], tnorm_q1[self.cancer2], ptumor_q1[self.cancer4], tnorm_q1[self.cancer4]] #cancer 3 -> cancers 1 and 2, 4
                    hypo_quantiles4 = [ptumor_q1[self.cancer1], tnorm_q1[self.cancer1], ptumor_q1[self.cancer2], tnorm_q1[self.cancer2], ptumor_q1[self.cancer3], tnorm_q1[self.cancer3]] #cancer 4 -> cancers 1 and 2, 3                    
                    

                    #list of quantiles - for tissue-specific hyper and hypo check
                    tis_hyper_quantiles1 = [ptumor_q3[self.cancer2], tnorm_q3[self.cancer2], ptumor_q3[self.cancer3], tnorm_q3[self.cancer3], ptumor_q3[self.cancer4], tnorm_q3[self.cancer4]] # cancer 1 -> cancers 2 and 3, 4
                    tis_hyper_quantiles2 = [ptumor_q3[self.cancer1], tnorm_q3[self.cancer1], ptumor_q3[self.cancer3], tnorm_q3[self.cancer3], ptumor_q3[self.cancer4], tnorm_q3[self.cancer4]] #cancer 2 -> cancers 1 and 3, 4
                    tis_hyper_quantiles3 = [ptumor_q3[self.cancer1], tnorm_q3[self.cancer1], ptumor_q3[self.cancer2], tnorm_q3[self.cancer2], ptumor_q3[self.cancer4], tnorm_q3[self.cancer4]]  #cancer 3 -> cancers 1 and 2, 4
                    tis_hyper_quantiles4 = [ptumor_q3[self.cancer1], tnorm_q3[self.cancer1], ptumor_q3[self.cancer2], tnorm_q3[self.cancer2], ptumor_q3[self.cancer3], tnorm_q3[self.cancer3]] #cancer 4 -> cancers 1 and 2, 3
                    
                    tis_hypo_quantiles1 = [ptumor_q1[self.cancer2], tnorm_q1[self.cancer2], ptumor_q1[self.cancer3], tnorm_q1[self.cancer3], ptumor_q1[self.cancer4], tnorm_q1[self.cancer4]] # cancer 1 -> cancers 2 and 3, 4
                    tis_hypo_quantiles2 = [ptumor_q1[self.cancer1], tnorm_q1[self.cancer1], ptumor_q1[self.cancer3], tnorm_q1[self.cancer3], ptumor_q1[self.cancer4], tnorm_q1[self.cancer4]] #cancer 2 -> cancers 1 and 3, 4
                    tis_hypo_quantiles3 = [ptumor_q1[self.cancer1], tnorm_q1[self.cancer1], ptumor_q1[self.cancer2], tnorm_q1[self.cancer2], ptumor_q1[self.cancer4], tnorm_q1[self.cancer4]]  #cancer 3 -> cancers 1 and 2, 4
                    tis_hypo_quantiles4 = [ptumor_q1[self.cancer1], tnorm_q1[self.cancer1], ptumor_q1[self.cancer2], tnorm_q1[self.cancer2], ptumor_q1[self.cancer3], tnorm_q1[self.cancer3]] #cancer 4 -> cancers 1 and 2, 3
                                        
                    #Cancer 1: check for hypermethylation in cancer1 (compare cancer1 with cancer2, cancer3 and cancer4)
                    if hyperCheck(tnorm_q3[self.cancer1], ptumor_q1[self.cancer1], hyper_quantiles1) == True:
                        result = [probe, self.cancer1, "+"]
                        num1+=1
                        
                    #Cancer 1:  check for hypomethylation in cancer1 (compare cancer1 with cancer2, cancer3 and cancer4)
                    elif hypoCheck(tnorm_q1[self.cancer1], ptumor_q3[self.cancer1], hypo_quantiles1) == True:
                        result = [probe, self.cancer1, "-"]
                        num2+=1
                
                    #Cancer 2: check for hypermethylation in cancer2 (compare cancer2 with cancer1, cancer3 and cancer4)
                    if hyperCheck(tnorm_q3[self.cancer2], ptumor_q1[self.cancer2], hyper_quantiles2) == True:
                        result = [probe, self.cancer2, "+"]
                        num1+=1

                    #Cancer 2:  check for hypomethylation in cancer2 (compare cancer2 with cancer1, cancer3 and cancer4)
                    elif hypoCheck(tnorm_q1[self.cancer2], ptumor_q3[self.cancer2], hypo_quantiles2) == True:
                        result = [probe, self.cancer2, "-"]                                
                        num2+=1
                        
                    #Cancer 3: check for hypermethylation in cancer3 (compare cancer3 with cancer1, cancer2 and cancer4)
                    if hyperCheck(tnorm_q3[self.cancer3], ptumor_q1[self.cancer3], hyper_quantiles3) == True:
                        result = [probe, self.cancer3, "+"]
                        num1+=1

                    #Cancer 3:  check for hypomethylation in cancer3 (compare cancer3 with cancer1, cancer2 and cancer4)
                    elif hypoCheck(tnorm_q1[self.cancer3], ptumor_q3[self.cancer3], hypo_quantiles3) == True:
                        result = [probe, self.cancer3, "-"]                                
                        num2+=1
                        
                    #Cancer 4: check for hypermethylation in cancer4 (compare cancer4 with cancer1, cancer2 and cancer3)
                    if hyperCheck(tnorm_q3[self.cancer4], ptumor_q1[self.cancer4], hyper_quantiles4) == True:
                        result = [probe, self.cancer4, "+"]
                        num1+=1
                        
                    #Cancer 4: check for hypomethylation in cancer4 (compare cancer4 with cancer1, cancer2 and cancer3)
                    elif hypoCheck(tnorm_q1[self.cancer4], ptumor_q3[self.cancer4], hypo_quantiles4) == True:
                        result = [probe, self.cancer4, "-"]                                
                        num2+=1
                    
                    #Cancer 1: check for tissue-specific hypermethylation in cancer 1 (compare cancer1 with cancer2, cancer3 and cancer4)
                    elif tissue_hyperCheck([ptumor_q1[self.cancer1], tnorm_q1[self.cancer1]], tis_hyper_quantiles1):
                        result = [probe, self.cancer1, "#+"]
                        num3+=1
                        
                    #Cancer 1 check for tissue=specific hypomethylation in cancer1 (compare cancer1 with cancer2, cancer3 and cancer4)
                    elif tissue_hypoCheck([ptumor_q3[self.cancer1], tnorm_q3[self.cancer1]], tis_hypo_quantiles1):
                        result = [probe, self.cancer1, "#-"]
                        num4+=1
                        
                    #Cancer 2: check for tissue-specific hypermethylation in cancer 2 (compare cancer2 with cancer1, cancer3 and cancer4)
                    elif tissue_hyperCheck([ptumor_q1[self.cancer2], tnorm_q1[self.cancer2]], tis_hyper_quantiles2):
                        result = [probe, self.cancer2, "#+"]
                        num3+=1
                        
                    #Cancer 2: check for tissue-specific hypomethylation in cancer 2 (compare cancer2 with cancer1, cancer3 and cancer4)
                    elif tissue_hypoCheck([ptumor_q3[self.cancer2], tnorm_q3[self.cancer2]], tis_hypo_quantiles2):
                        result = [probe, self.cancer2, "#-"]                    
                        num4+=1                

                    #Cancer 3: check for tissue-specific hypermethylation in cancer 3 (compare cancer3 with cancer1, cancer2 and cancer4)
                    elif tissue_hyperCheck([ptumor_q1[self.cancer3], tnorm_q1[self.cancer3]], tis_hyper_quantiles3):
                        result = [probe, self.cancer3, "#+"]
                        num3+=1
                        
                    #Cancer 3: check for tissue-specific hypomethylation in cancer 3 (compare cancer3 with cancer1, cancer2 and cancer4)
                    elif tissue_hypoCheck([ptumor_q3[self.cancer3], tnorm_q3[self.cancer3]], tis_hypo_quantiles3):
                        result = [probe, self.cancer3, "#-"]                    
                        num4+=1      

                    #Cancer 4: check for tissue-specific hypermethylation in cancer 4 (compare cancer4 with cancer1, cancer2 and cancer4)
                    elif tissue_hyperCheck([ptumor_q1[self.cancer4], tnorm_q1[self.cancer4]], tis_hyper_quantiles4):
                        result = [probe, self.cancer4, "#+"]
                        num3+=1
                        
                    #Cancer 4: check for tissue-specific hypomethylation in cancer 4 (compare cancer4 with cancer1, cancer2 and cancer4)
                    elif tissue_hypoCheck([ptumor_q3[self.cancer4], tnorm_q3[self.cancer4]], tis_hypo_quantiles4):
                        result = [probe, self.cancer4, "#-"]                    
                        num4+=1    
                        
                    #Pancancer hypermethylation: check for pan-cancer hyper methylation
                    elif pancancer_hyperCheck(list(ptumor_q1.values()), list(tnorm_q3.values())) == True: #ptumor_q1s, tnorm_q1s
                        result = [probe, "4cancers", "*+"]
                        num5+=1
                        
                    #Pancancer hypomethylation: check for pan-cancer hypo methylation
                    elif pancancer_hypoCheck(list(ptumor_q3.values()), list(tnorm_q1.values())) == True:
                        result = [probe, "4cancers", "*-"]
                        num6 +=1   
						
                #only add if result is found and tumor specific 
                if result:
                    self.resultsList.append(result)                        
            
            print ("cancer hypermethylation:", num1)
            print ("cancer hypomethylation:", num2)
            print ("tissue hypermethylation", num3)
            print ("tissue hypomethylation:", num4)
            print ("pan-cancer hypermethylation", num5)
            print ("pan-cancer hypomethylation:", num6)
            
            return ''   
        ####################### s     
        compare_probes() #compares single cancer with other cancers
        # compare_combined_probes() #compares multiple cancers with other cancers
               
        return self.resultsList      					
