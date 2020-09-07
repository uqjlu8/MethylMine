# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 19:10:41 2020

@author: s4142554
#create heatmap with colorbars
repeat mDist_024_Scripts_v1 with MethylMine8
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
# from BioStats import MatrixStatSummary
from methylMine import DiffAnalysis, filtMatrix, heatmap

###############################################################################
script = os.path.basename(__file__).split('.',1)[0]
script_dir = os.path.dirname(__file__) #<-- absolute dir the script is in
###############################################################################
MatrixFolder = "V:/Trau/Jennifer/TCGA_MatrixFiles/" #folder where the matrices for each cancer and tissue type is kept

ResultsFolder = "C:/Users/s4142554/Desktop/mDist_025/mDist_025_Results/" #folder where the result files will be stored.

###############################################################################
def MethyMineProcess(cancerList, MatrixFolder, folder_suffix, matrix_suffix, ResultsFolder, lower_beta, upper_beta):

    Tissues = ["PrimaryTumor", "TissueNormal"]

    ptumor_num = 50
    tnorm_num = 50

    Dict = {
        'cancer':[],
        "+" : [], # cancer-specific
        "-" : [], # cancer-specific
        "#+" : [], # tisue-specific
        "#-" : [], # tisue-specific
        "*+": [], # pan-cancer
        "*-":[] # pan-cancer
        }

    for cancer in cancerList:
        Dict['cancer'].append(cancer)

        print (cancer)
        
        cancers = [cancer]

        MethProbeDict = {"+" : [], # cancer-specific
                         "-" : [], # cancer-specific
                         "#+" : [], # tisue-specific
                         "#-" : [], # tisue-specific
                         "*+": [], # pan-cancer
                         "*-":[] # pan-cancer
                         }

    QuantileResults = DiffAnalysis(cancerList, MatrixFolder, folder_suffix, matrix_suffix, ResultsFolder, lower_beta, upper_beta).compare()#, ptumor_num, tnorm_num).compare()

    for line in QuantileResults:
        probe, cancer, meth = line
        MethProbeDict[meth].append(probe)

    for meth in MethProbeDict:
        Dict[meth].extend(MethProbeDict[meth])

    ProbesList = [n[0] for n in QuantileResults]
    dataframes = []

    # # filter through matrix files
    for i, cancer in enumerate(cancerList):
        
        cancerMatrixFolder = os.path.join(script, MatrixFolder, "_".join([cancer, "Matrix"]))
    
        pTumor_Matrix = os.path.join(script, cancerMatrixFolder, "_".join([cancer, Tissues[0],"matrix.txt"]))
        tNorm_Matrix = os.path.join(script, cancerMatrixFolder, "_".join([cancer, Tissues[1],"matrix.txt"]))
        
        print ("Creating matrix", cancer, Tissues[0])
        
        df = filtMatrix(ProbesList, cancer, Tissues[0], pTumor_Matrix, ptumor_num)

        dataframes.append(df)  #add in Primary Tumor
        
        print ("Creating matrix", cancer, Tissues[1])
        
        df = filtMatrix(ProbesList, cancer, Tissues[1], tNorm_Matrix, tnorm_num)
        
        dataframes.append(df) #add in Tissue normal

    # concat all four df together 
    dataframes = pd.concat(dataframes,axis=1)
    tag = "&".join(cancerList) if len(cancerList) >=2 else cancerList[0]
    resultDF = os.path.join(script, ResultsFolder, "_".join(["methylMine", str(lower_beta), str(upper_beta), tag, str(ptumor_num), str(tnorm_num),"filtMatrix.txt"]))
    result_fig = os.path.join(script, ResultsFolder, "_".join(["methylMine", str(lower_beta), str(upper_beta), tag, str(ptumor_num), str(tnorm_num),"filtMatrix.png"]))
    rowFile = os.path.join(script, ResultsFolder, "_".join(["methylMine", tag, str(ptumor_num), str(tnorm_num),"probeResults.txt"]))
    
    # #write dataframes to text for later heatmap
    print ("Writing df to csv")
    dataframes.to_csv(resultDF, sep="\t")    
    print ("Generating heatmaps")
    
    title = "\t".join(["ProbeID", "cancer", "meth"])+"\n"
    
    with open(rowFile, 'w') as file:
        file.write(title)
        for line in QuantileResults:
            # probe, cancer, meth = line
            line = "\t".join(line)+"\n"
            file.write(line)
    
    
    heatmap(dataframes, result_fig)
    # return Dict
    return MethProbeDict

def MatrixHeatmap(MethProbeDict, cancerList, MatrixFolder, folder_suffix, matrix_suffix, ResultsFolder):
    
    ProbesList = []

    Tissues = Tissues = ["PrimaryTumor", "TissueNormal"]
    
    for meth in MethProbeDict:
        ProbesList.extend(MethProbeDict[meth])
        
    dataframes = [] #for tumor matrix
    

    # # filter through matrix files
    for i, cancer in enumerate(cancerList):
        
        cancerMatrixFolder = os.path.join(script, MatrixFolder, "_".join([cancer, "Matrix"]))
    
        pTumor_Matrix = os.path.join(script, cancerMatrixFolder, "_".join([cancer, Tissues[0],"matrix.txt"]))
        tNorm_Matrix = os.path.join(script, cancerMatrixFolder, "_".join([cancer, Tissues[1],"matrix.txt"]))
        
        print ("Creating matrix", cancer, Tissues[0])
        
        df = filtMatrix(ProbesList, cancer, Tissues[0], pTumor_Matrix, ptumor_num)

        dataframes.append(df)  #add in Primary Tumor
        
        # print ("Creating matrix", cancer, Tissues[1])
        
        # df = filtMatrix(ProbesList, cancer, Tissues[1], tNorm_Matrix, tnorm_num)
        
        # dataframes.append(df) #add in Tissue normal

    # concat all four df together 
    dataframes = pd.concat(dataframes,axis=1)
    tag = "&".join(cancerList) if len(cancerList) >=2 else cancerList[0]
    resultDF = os.path.join(script, ResultsFolder, "_".join(["1.mDist_025_methylMine8", tag, str(ptumor_num), str(tnorm_num),"filtMatrix.txt"]))
    result_fig = os.path.join(script, ResultsFolder, "_".join(["1.mDist_025_methylMine8", tag, str(ptumor_num), str(tnorm_num),"filtMatrix.png"]))
    
    # #write dataframes to text for later heatmap
    print ("Writing df to csv")
    dataframes.to_csv(resultDF, sep="\t")    
    print ("Generating heatmaps")
    
    heatmap(dataframes, result_fig)        
    
    

    
    return  MethProbeDict

###############################################################################

lower_beta = 0.3
upper_beta = 0.6

ptumor_num = 50
tnorm_num = 50

#for full matrix
folder_suffix = "Matrix" #TCGA-ACC_Matrix
matrix_suffix = "matrix.txt"
cancers = [
        "TCGA-BRCA", 
        "TCGA-COAD", 
        "TCGA-LUSC", 
        "TCGA-PRAD", 
        ]



###############################################################################
start = time.time()
probeDict = MethyMineProcess(cancers, MatrixFolder, folder_suffix, matrix_suffix, ResultsFolder, lower_beta, upper_beta)
result = MatrixHeatmap(probeDict, cancers, MatrixFolder, folder_suffix, matrix_suffix, ResultsFolder)
finish = time.time()

print ("Finished", (finish - start)/60, "min")
