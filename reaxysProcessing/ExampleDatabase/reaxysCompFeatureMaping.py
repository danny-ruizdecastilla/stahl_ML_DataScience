import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
from rdkit.Chem.PandasTools import LoadSDF
import matplotlib.pyplot as plt
import random
import hdbscan
from itertools import combinations
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from rdkit.Chem import Descriptors, MACCSkeys
from rdkit.Chem import AllChem
import datetime
from hdbscan import HDBSCAN
import umap
from expt2_create_features import *
from expt2_feature_filtering import *
#Danny Ruiz de Castilla 02.01.2025

def createCSV(dataFrame , saveDir, saveStr: str):
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    dataFrame.to_csv(saveDir + "/" + saveStr + ".csv", index=True) 
    return dataFrame
def savePNG(df, saveDir:str , saveStr, saveStr_ , cluster:bool ):
    if cluster:
        xLabel = saveStr_ + "1"
        yLabel = saveStr_ + "2"
        clusterLabel = "Clusters"
        xArray = np.array(df[xLabel])
        yArray = np.array(df[yLabel])
        clusterArr = np.array(df[clusterLabel])
        plt.figure(figsize=(8, 6), dpi=100)  # Set figure size and DPI

        # Scatter plot for clustered points
        uniques = np.unique(clusterArr)
        for cluster in uniques:
            if cluster == -1:  
                # Plot noise points in gray
                plt.scatter(xArray[clusterArr == -1], yArray[clusterArr == -1], 
                    color='gray', label="Noise", alpha=0.5, edgecolors='k')
            else:
                # Plot clustered points
                plt.scatter(xArray[clusterArr == cluster], yArray[clusterArr == cluster], color=plt.cm.viridis(cluster / (len(uniques))), 
                label="Cluster " + str(cluster), alpha=0.7, edgecolors='k')

        plt.xlabel(xLabel, fontsize=14, fontweight='bold', color='black')
        plt.ylabel(yLabel, fontsize=14, fontweight='bold', color='black')
        plt.title(saveStr_ + " visualization for " + saveStr, fontsize=18, fontweight='bold', color='navy')
        plt.legend(title='Clusters', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xticks(fontsize=14, color='black')
        plt.yticks(fontsize=14, color='black')
        #plt.grid(True, linestyle='--', alpha=0.5)
        plt.savefig(saveDir + "/saveStr_" + " visualization for" + "saveStr.png" , dpi=300, bbox_inches='tight')


    else:
        xLabel = saveStr_ + "1"
        yLabel = saveStr_ + "2"
        xList = list(df[xLabel])
        yList = list(df[yLabel])
        dpi = 300
        plt.figure(figsize=(800 / dpi, 600 / dpi), dpi=dpi)  # 800x600 pixels
        plt.scatter(xList, yList, c='blue', alpha=0.7, edgecolors='k')
        plt.xlabel(xLabel, fontsize=16, fontweight='bold', color='black')
        plt.ylabel(yLabel, fontsize=16, fontweight='bold', color='black')
        plt.title(saveStr_ + " Cluster visualization for " + saveStr, fontsize=18, fontweight='bold', color='navy')
        plt.xticks(fontsize=14, color='black')
        plt.yticks(fontsize=14, color='black')
        #plt.grid(True, linestyle='--', alpha=0.5)
        plt.savefig(saveDir + "/saveStr_" + " Cluster visualization for" + "saveStr.png" , dpi=300, bbox_inches='tight')

    return saveStr
def DimensionalityReduction(X, Y , sampleScale1, partition1, sampleScale2 , catSplit , saveStr1 , saveDir1 ):
    #PCA then UMAP
    population = len(Y)
    scaler = StandardScaler()
    scaledX = scaler.fit_transform(X)
    pca = PCA(n_components=partition1)  # Keep top 10 principal components
    xPCA= pca.fit_transform(scaledX)

    if catSplit != 1:
        minSamples1 = int(np.sqrt(population))
        #minSamples1 = 0.01*float(len(Y))
        minClusterSize1 =  int(minSamples1*sampleScale1)
    clusterer = HDBSCAN(min_cluster_size=minClusterSize1, min_samples=minSamples1, metric='euclidean')
    clusterLabels = clusterer.fit_predict(xPCA)
    
    #hdbNoise = list(clusterLabels).count(-1)
    #print(hdbNoise)
    numClusters = len(set(clusterLabels)) - (1 if -1 in clusterLabels else 0)

    clusterCounts = uniqueIntCounter(list(clusterLabels))
    umapNeighbors = int(np.mean(clusterCounts.values())* sampleScale2)

    reducer = umap.UMAP(n_components=2, n_neighbors=umapNeighbors, min_dist=0.1, metric='euclidean')
    reducedX = reducer.fit_transform(xPCA)
    dfUMAP = pd.DataFrame(reducedX, columns=['UMAP1', 'UMAP2'])

    xUMAP = list(dfUMAP['UMAP1'])
    yUMAP = list(dfUMAP['UMAP2'])
    saveStr2 = "UMAP"
    savePNG(dfUMAP, saveDir1 , saveStr1, saveStr2 , cluster = False  )

    if numClusters > 6:
        clusterRange = np.arange(numClusters-6, numClusters+6)
    else:
        clusterRange = np.arange(1, 11)
    inertiaList = []

    for k in clusterRange:
        initInertia = []
        count = 0
        while count < 4:
            kmeans = KMeans(n_clusters=k, random_state=42)
            kmeans.fit(reducedX)
            initInertia.append(kmeans.inertia_)
            count +=1 
        inertia = np.mean(initInertia)
        inertiaList.append(inertia)
    specialK = needleAlg(list(clusterRange) , inertiaList)
    kmeansMAST = KMeans(n_clusters = specialK, random_state = 42 )
    kmeansLabels = kmeansMAST.fit(reducedX)

    dfUMAP['Clusters'] = kmeansLabels
    createCSV(dfUMAP , saveDir1, saveStr1 + "_clustered")
    savePNG(xUMAP , yUMAP , saveDir1 , saveStr1 + "_clustered"  )

    return dfUMAP
def slope(m , x, b):
        y = m*x + b
        return y
def needleAlg(xList , yList):
    x0 = xList[0]
    y0 = yList[0]

    x1 = xList[-1]
    y1 = yList[-1]
    m_ = (y1 - y0) / (x1 - x0)
    b_ = y0 - m_ * x0

    diffList = []
    for i in range (xList):
        x_ = xList[i]
        yPred = slope(m_ , x_ , b_) 
        diff = np.abs(yPred - yList[i])
        diffList.append(diff)
    elbowInd = diffList.index(max(diffList))
    elbowX = xList[elbowInd]
    return elbowX
def uniqueIntCounter(list1):
    uniqueDict = {}
    for i in range (len(list1)):
        num = list1[i]
        if num in uniqueDict:
            uniqueDict[num] +=1
        else:
            uniqueDict[num] = 1
    return uniqueDict
if __name__ == "__main__":
#input an xlsx dataframe and depending on what features you want to calculate, export .csv of features, performs feature elimination, maps features, and creates clusters 
    mainDir = str(sys.argv[1]) #Dataframe of Smiles and yields and ID numbers, and citations 
    rdkitAllow = int(sys.argv[2])
    moldrAllow = int(sys.argv[3])
    maccAllow = int(sys.argv[4])
    catSplit = int(sys.argv[5])
    sampleScale1 = float(sys.argv[6])
    chemStr = str(sys.argv[7])
    inputDF = pd.read_excel(mainDir + "MasterDataFrame" + str(chemStr) + ".xlsx")
    
    featureTypes = []   #defines options for features you want 
    if rdkitAllow == 1:
        rdkitSave = "rdkit"
        featureTypes.append(rdkitSave)
    if moldrAllow == 1:
        moldrSave = "moldred"
        featureTypes.append(moldrSave)
    if maccAllow == 1:
        maccSave = "MACCSkeys"
        featureTypes.append(maccSave)
    if catSplit == 1: #Attempts to create initial clusterSize and minSample variables by using the catalyst splitting 
        catalysts = list(inputDF['ChemType'])
        catalystCounts = uniqueIntCounter(catalysts)
        minClusterSize1 =  min(catalystCounts.values()) 
        minSamples1 = int(minClusterSize1/sampleScale1)
    for feature in featureTypes:
        #creates necessary features of interest, and generates plots 
        
    


    

    
    





    
