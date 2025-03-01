import sys
import glob
import numpy as np
import os
#from rdkit import Chem
import pandas as pd
#from rdkit.Chem.PandasTools import LoadSDF
import matplotlib.pyplot as plt
import random
import hdbscan
from itertools import combinations
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
#from rdkit.Chem import Descriptors, MACCSkeys
#from rdkit.Chem import AllChem
from hdbscan import HDBSCAN
from expt2_feature_filtering import *
#Danny Ruiz de Castilla 02.28.2025

#pcaDict = {canonicalSmiles : [ yield]} for local Chemistries 
#smilesDict = {canonicalSmiles : [ xPCA , yPCA]}
def dimensionalityReduction(inputDF):
    smiles = list(inputDF["SMILES"])
    X = inputDF.drop(inputDF.columns[:2], axis=1)

    scaler = StandardScaler()
    scaledX = scaler.fit_transform(X)


    pca = PCA(n_components = 2 ,svd_solver="full" )
    xPCA = pca.fit_transform(scaledX)

    dfPCA = pd.DataFrame(xPCA, columns = ["PCA1" , "PCA2"] )
    xAxis = list(dfPCA["PCA1"])
    yAxis = list(dfPCA["PCA2"])
    pcaDict = {}
    for i , smile in enumerate (smiles):
        
        pcaDict[smile] = [xAxis[i] , yAxis[i]]

        
        return pcaDict   
def findHighlights(smilesDict  ,  partition1 , chemistryDict):
    #progressive modification to chemistryDict by feeding the coordinates within the keys of smilesDict 
    substratePool = list(smilesDict.keys())
    coordinates = list(smilesDict.values())
    methodSubstrates = list(chemistryDict.keys())
    methodValues = list(chemistryDict.values())
    highlightDict = {}
    for i , substrate in enumerate(methodSubstrates):
        if substrate in substratePool:
            #Check if yield meets threshold 
            poolIndex = substratePool.index(substrate)
            parameter = methodValues[i]
            if parameter >= partition1:
                xPCA = list(coordinates[poolIndex])[0]
                yPCA = list(coordinates[poolIndex])[1]
                highlightDict[substrate] = [xPCA , yPCA]
                del smilesDict[substrate]

    return highlightDict , smilesDict

def makePlots(pcaDict ,  partitionList , chemistryDicts  ,chemistryStr , colors , partitionStr):
    
    for partition in partitionList:
        for i , chemistryDict in enumerate(chemistryDicts): #this is the highlighted chemistry in the plot
            highlightDict , smilesDict = findHighlights(pcaDict  ,  partition , chemistryDict )
            xHigh , yHigh = zip (*highlightDict.values())
            xBland , yBland = zip(*smilesDict.values())
            chemistryLabel = chemistryStr[i]

            xLabel = "PCA1"
            yLabel = "PCA2"
            dpi = 300
            plt.figure(figsize=(800 / dpi, 600 / dpi), dpi=dpi)  # 800x600 pixels
            plt.scatter(list(xHigh), list(yHigh), c=colors[i], alpha=0.7 , s=10)
            plt.scatter(list(xBland), list(yBland), c="grey", alpha=0.14 , s=10)
            plt.xlabel(xLabel, fontsize=16, fontweight='bold', color='black')
            plt.ylabel(yLabel, fontsize=16, fontweight='bold', color='black')
            plt.title("Substrate Scope for " + str(chemistryLabel) + " at " + str(partition) + "%" + str(partitionStr), fontsize=18, fontweight='bold', color='navy')
            plt.xticks(fontsize=14, color='black')
            plt.yticks(fontsize=14, color='black')
            #plt.grid(True, linestyle='--', alpha=0.5)
            plt.savefig(str(chemistryLabel) + "at" + str(partition) + ".png", dpi=300, bbox_inches='tight')
    
def featureFiltering(outDir , X , feature_labels , featureStr):
    if not os.path.exists(outDir + "/" + str(featureStr) + "featureFiltering.dat"):
        with open(outDir + "/" + str(featureStr) + "featureFiltering.dat", "w") as f:
            f.write(f"Total starting feature count: {len(feature_labels)}")
            f.write("".join([f'\n\t{label}' for label in feature_labels]))
            
            X, feature_labels, dropped_features = remove_by_variance(X, feature_labels)
            text = "\n\n\nFeatures drop due to low variance: " + "".join([f'\n\t{label}' for label in dropped_features])
            f.write(text)

            X, feature_labels, drop_group = correlation_analysis(X, feature_labels, threshold=0.95)
            text = "\n\n\nFeatures drop due correlation: (STILL HAS ISSUES) "
            import json
            text += json.dumps(drop_group, indent=4).replace('\n', '\n\t')
            f.write(text)

            X, feature_labels  = spearmanr_correlation(X, feature_labels, threshold=0.95)
    return X , feature_labels
def transformations():
    
    return 
if __name__ == "__main__":

    chemistryDirs = str(sys.argv[1])
    datasetDir = str(sys.argv[2])

    initdataSets = glob.glob(datasetDir + "/*.csv")
    





