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
    plt.savefig(saveDir + "/Cluster"  + saveStr + str(int(klDiv)) +  ".png")
    plt.close()

def removeNonLetters(text):
    return re.sub(r'[^a-zA-Z0-9_]', '', text)
def histogramGenerator(df , partitionStr , partitionList , columnList, pngSaveDir ):

    for col in columnList:
        if col != partitionStr:
            col_ = col
            data1 = list(df[col_])
            if '%' in col:
                col = col.replace('%', '_')
            featureStr = removeNonLetters(col)
            
            if not os.path.exists(pngSaveDir + "/" + str(featureStr)):
                os.makedirs(pngSaveDir + "/" + featureStr)

            for partition in partitionList:
                partitionType = str(partition)
                data2 = list(df[df[partitionStr] == partition][col_])
                
                createHistograms(data1 , data2 , pngSaveDir + "/" + featureStr , partitionType)
                
if __name__ == "__main__":

    dataframeDir = str(sys.argv[1]) #Must point to a .csv
    pngSaveDir = str(sys.argv[2]) 
    clusterStr = str(sys.argv[3])#what is the name of the column containing all the cluster assingments 
    if not os.path.exists(pngSaveDir):
        os.makedirs(pngSaveDir)
    mainDF = pd.read_csv(dataframeDir)

    if "SMILES" in mainDF.columns: 
        mainDF = mainDF.drop(columns=["SMILES"])
    clusterList = np.unique(mainDF[clusterStr])
    featureList = list(mainDF.columns)
    histogramGenerator(mainDF , clusterStr , clusterList , featureList , pngSaveDir)
