import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
import matplotlib.pyplot as plt
import random
from itertools import combinations
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from scipy.stats import qmc
#Danny Ruiz de Castilla 02.12.2025
def strategicRandomSampling(df, partitionStr1, yieldParts, saveStr1 , saveStr2 , targetN):
    #should split up by partition index 
    partitionID = list(df[partitionStr1])

def getCanonical(maskedSMILES):
    canonicalList = []
    for smile in maskedSMILES:
        if isinstance(smile, str):
            mol = Chem.MolFromSmiles(smile)
            canonSMILES = Chem.MolToSmiles(mol, canonical=True)
            canonicalList.append(canonSMILES)
        else:
            continue
    return canonicalList
def getSplits(df , N , split , pop):
    dfList = []
    numSamples = {}
    partitionID = df[split].unique()
    for partition in partitionID:
        dfPart = df[df[split] == partition]
        dfList.append(dfPart)
        if partition in numSamples:
            numSamples[partition ] +=1
        else:
            numSamples[partition ] = 1
    scale = N / pop
    for num in numSamples:
        numSamples[num] = int(numSamples[num] * scale)
    
    return dfList , numSamples
if __name__ == "__main__":
    mainDir = str(sys.argv[1])
    maskedDir = str(sys.argv[2])
    mainDFStr = str(sys.argv[3])
    partitionStr1 = str(sys.argv[4])
    partitionStr2 = str(sys.argv[5])
    targetN = int(sys.argv[6]) 

    yieldRanges =[0 , 15 , 30 , 45 , 60 , 75 , 90]

    maskedDF = pd.read_csv(maskedDir + ".csv")
    inputDF = pd.read_csv(mainDir + "/" + mainDFStr +  ".csv")
    sampleNum = len(inputDF)
    maskedSMILES = maskedDF['SMILES']
    canonMASK = getCanonical(maskedSMILES)

    mainSMILES = inputDF['SMILES']
    canonSMILES = getCanonical(mainSMILES)
    inputDF["Canonical"] = canonSMILES
    dfSplits = getSplits()
