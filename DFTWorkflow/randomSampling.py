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
def strategicRandomSampling(df, targetN, yieldRanges , partitionStr1 ):
    #yieldRanges =[0 , 15 , 30 , 45 , 60 , 75 , 90 , 100]
    sampleCount = 0
    dfSplit = []
    for i in range(len(yieldRanges)-1):
        yieldRange = [yieldRanges[i] , yieldRanges[i + 1]]
        df_ = df[(df[partitionStr1] < yieldRange[1]) & (df[partitionStr1] > yieldRange[0])]
        if len(df_) != 0:
            dfSplit.append(df_)
    dfCounter = len(dfSplit) - 1
    count = 0
    while True:
        sample = dfSplit[count].copy()
        smiles = list(sample['Canonical'])
        if len(smiles) != 0:
            smile = random.choice(smiles)
            if smile not in masterSMILESLIST and smile not in canonMASK:
                masterSMILESLIST.append(smile)
                sampleCount += 1
                dfSplit[count] = sample[sample['Canonical'] != smile]
        count += 1
        if count == dfCounter:
            count = 0
        

        if sampleCount == targetN:
            break

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
    for num in numSamples:
        numSamples[num] = int(numSamples[num]/ pop * N)
    
    return dfList , numSamples
def createCSV(dataFrame , saveDir, saveStr: str):
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    dataFrame.to_csv(saveDir + "/" + saveStr + ".csv", index=True) 
    return dataFrame
def createXLSX(smilesList , saveStr):
    smilesArray = np.array(smilesList)
    dim1 = len(smilesList)
    numbers = []
    for i in range(dim1):  # Generate numbers from 1 to 10
        numb = str(i).zfill(6)
        numbs = "Alkene_" + numb
        numbers.append(numbs)
    columns = ["SMILES", "ID"]
    # Create the DataFrame
    df = pd.DataFrame({columns[0]: smilesArray, columns[1]: numbers})

    df.to_excel(mainDir  + "/" + saveStr +  "DFTInput.xlsx", index=False) 
if __name__ == "__main__":
    mainDir = str(sys.argv[1])
    maskedDir = str(sys.argv[2])
    mainDFStr = str(sys.argv[3])
    partitionStr = str(sys.argv[4])
    yieldStr = str(sys.argv[5])
    targetN = int(sys.argv[6]) 

    yieldRanges =[0 , 15 , 30 , 45 , 60 , 75 , 90]

    maskedDF = pd.read_csv(maskedDir)
    inputDF = pd.read_csv(mainDir + "/" + mainDFStr +  ".csv")
    sampleNum = len(inputDF)
    maskedSMILES = maskedDF['SMILES']
    canonMASK = getCanonical(maskedSMILES) #list of all current SMILES

    mainSMILES = inputDF['SMILES']
    canonSMILES = getCanonical(mainSMILES)
    inputDF["Canonical"] = canonSMILES
    dfSplits , sampleDict = getSplits(inputDF , targetN , partitionStr , sampleNum)
    samples = list(sampleDict.values())
    masterSMILESLIST = []

    for i in range (len(dfSplits)):
        df = dfSplits[i]  
        sample = samples[i]

        strategicRandomSampling(df , sample , yieldRanges, yieldStr )

    createXLSX(masterSMILESLIST ,mainDFStr )

    finalDF = inputDF[inputDF["Canonical"].isin(masterSMILESLIST)]
    df = createCSV(finalDF , mainDir + "/DFTInput/" , mainDFStr)