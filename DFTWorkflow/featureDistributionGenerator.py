import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
#from rdkit.Chem.PandasTools import LoadSDF
import matplotlib.pyplot as plt
import random
import hdbscan
import re
from itertools import combinations
from featureMaping import savePNG , createCSV
#Danny Ruiz de Castilla 02.28.2025
#Imports cluster identifier and feature Dataframe to create 
def getEpsilon(list1):
    list1 = np.array(list1) 
    if np.any(list1 == 0): 
        minNumb = np.min(np.abs(list1)[list1 != 0]) 
        factoredNumb = minNumb * 0.00001
        return factoredNumb
    else:
        return 0
def getKLDiv(dist1, dist2):
    klDiv = np.sum(dist1 * np.log(dist1 / dist2))
    return klDiv
def createHistograms(data1 , data2 , saveDir , saveStr):
    hist1, bins_ = np.histogram(data1, bins="auto", density=True)
    epsilon = getEpsilon(hist1)

    hist2 , _ = np.histogram(data1, bins=bins_, density=True)

    hist1 += epsilon
    hist2 += epsilon

    hist1 /= hist1.sum()
    hist2 /= hist2.sum()

    klDiv = getKLDiv(hist1 , hist2)
    plt.figure(figsize=(8, 5))
    plt.hist(data1, bins=bins_, density=True, alpha=0.6, label="Entire Substrate Scope", color='blue')
    plt.hist(data2, bins=bins_, density=True, alpha=0.6, label="Cluster" + saveStr, color='red')
    plt.title(f"Histogram Comparison\nKL Divergence: {klDiv:.4f}")
    plt.legend()
    plt.savefig(saveDir + "/"  + saveStr + str(int(klDiv)) +  ".png")
    plt.close()

def removeNonLetters(text):
    return re.sub(r'[^a-zA-Z0-9_]', '', text)
def histogramGenerator(df ,  columnList, pngSaveDir , mainDF ):

    if smileStr in columnList:
        columnList.remove(smileStr)
    if pc1Str in columnList:
        columnList.remove(pc1Str)
    if pc2Str in columnList:
        columnList.remove(pc2Str)

    for col in columnList:
        col_ = col
        data1 = list(mainDF[col_])
        if '%' in col:
            col = col.replace('%', '_')
        featureStr = removeNonLetters(col)
        data2 = list(df[col_])
        createHistograms(data1 , data2 , pngSaveDir, featureStr)
def partition(bounds1, bounds2, dfMAST, range1Str, range2Str):
    dim1 = dfMAST[range1Str]
    dim2 = dfMAST[range2Str]
    df = pd.DataFrame(columns=dfMAST.columns)  # Ensures correct structure

    for i, row_ in enumerate(dfMAST.itertuples(index=False, name=None)):  
        if bounds1[0] <= dim1[i] <= bounds1[1] and bounds2[0] <= dim2[i] <= bounds2[1]:
            df = pd.concat([df, pd.DataFrame([row_], columns=dfMAST.columns)], ignore_index=True)

    return df, list(df.columns)

                
if __name__ == "__main__":
    dataframeDir = str(sys.argv[1]) #Must point to a .csv no longer cluster csv 
    pc1Min = float(sys.argv[2])
    pc1Max = float(sys.argv[3])
    pc2Min = float(sys.argv[4])
    pc2Max = float(sys.argv[5])
    pc1Str = str(sys.argv[6])
    pc2Str = str(sys.argv[7])
    smileStr = str(sys.argv[8])
    fileStr = "pc1_" + str(pc1Min) + "_" + str(pc1Max) + "_pc2_" + str(pc2Min) + "_" + str(pc2Max)
    saveDir = dataframeDir + "/" + fileStr
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    featurePCA = glob.glob(dataframeDir + "/*.csv")[0]
    mainDF = pd.read_csv(featurePCA)

    pc1Range = [pc1Min , pc1Max]
    pc2Range = [pc2Min , pc2Max]
    newDF , columnNames = partition(pc1Range , pc2Range , mainDF , pc1Str , pc2Str)
    createCSV(newDF, saveDir , fileStr + "_partitionPCA")



    histogramGenerator(newDF , columnNames, saveDir , mainDF )
