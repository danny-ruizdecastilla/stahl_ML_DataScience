import sys
#import glob
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
import umap
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
                    color='gray', label="Noise", alpha=0.5, s =10)
            else:
                # Plot clustered points
                plt.scatter(xArray[clusterArr == cluster], yArray[clusterArr == cluster], color=plt.cm.viridis(cluster / (len(uniques))), 
                label="Cluster " + str(cluster), alpha=0.7 , s = 10)

        plt.xlabel(xLabel, fontsize=14, fontweight='bold', color='black')
        plt.ylabel(yLabel, fontsize=14, fontweight='bold', color='black')
        plt.title(saveStr_ + " visualization for " + saveStr, fontsize=18, fontweight='bold', color='navy')
        plt.legend(title='Clusters', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xticks(fontsize=14, color='black')
        plt.yticks(fontsize=14, color='black')
        #plt.grid(True, linestyle='--', alpha=0.5)
        plt.savefig(saveDir  +"/" + saveStr_ + " Cluster visualization for " + saveStr + ".png" , dpi=300, bbox_inches='tight')


    else:
        xLabel = saveStr_ + "1"
        yLabel = saveStr_ + "2"
        xList = list(df[xLabel])
        yList = list(df[yLabel])
        dpi = 300
        plt.figure(figsize=(800 / dpi, 600 / dpi), dpi=dpi)  # 800x600 pixels
        plt.scatter(xList, yList, c='blue', alpha=0.7 , s=10)
        plt.xlabel(xLabel, fontsize=16, fontweight='bold', color='black')
        plt.ylabel(yLabel, fontsize=16, fontweight='bold', color='black')
        plt.title(saveStr_ + " Cluster visualization for " + saveStr, fontsize=18, fontweight='bold', color='navy')
        plt.xticks(fontsize=14, color='black')
        plt.yticks(fontsize=14, color='black')
        #plt.grid(True, linestyle='--', alpha=0.5)
        plt.savefig(saveDir + "/" + saveStr_ + " Cluster visualization for " + saveStr + ".png" , dpi=300, bbox_inches='tight')

    return saveStr
def locateNans(df):
    columns = df.columns.tolist()
    nanDict = {}
    for column in columns:
        col = df[column]
        indCol = col.isna().to_numpy().nonzero()[0]
        if len(indCol) != 0:
            #print("NANS" , indCol)
            nanDict[str(column)] = list(indCol)

    return nanDict
def eliminateNans(df , nanDict):
    #print(df.shape)
    for key in nanDict.keys():
        rows = list(nanDict[key])
        #print(rows)
        try:
            df = df.drop(rows)
        except KeyError:
            continue
    df = df.reset_index(drop=True)
    #print(df.shape)
    return df
def DimensionalityReduction(X,  clusterScale , saveDir , saveStr1,  ):
    scaler = StandardScaler()
    scaledX = scaler.fit_transform(X)
    partMax = int(np.sqrt(len(X)))
    partMin = int(partMax*clusterScale)

    pca = PCA(n_components = 2 ,svd_solver="full" )
    xPCA = pca.fit_transform(scaledX)
    dfPCA = pd.DataFrame(xPCA, columns = ["PCA1" , "PCA2"] )
    savePNG(dfPCA , saveDir , saveStr1 , "PCA" , cluster = False)
    clusterRange = np.arange(partMin, partMax, 2)
    print("clusterRange" , clusterRange)
    inertiaList = []
    for k in clusterRange:
        initInertia = []
        count = 0
        while count < 8:
            kmeans = KMeans(n_clusters = k, random_state = None)
            kmeans.fit(xPCA)
            initInertia.append(kmeans.inertia_)
            count += 1
        inertia = np.mean(initInertia)
        inertiaList.append(inertia)
    print(inertiaList)
    specialK = needleAlg(list(clusterRange) , inertiaList)
    kmeansMAST = KMeans( n_clusters = specialK , random_state = 42)
    kmeansMAST.fit(xPCA)
    dfPCA['Clusters'] = kmeansMAST.labels_
    createCSV(dfPCA, saveDir , saveStr1 + "_clustered")
    savePNG(dfPCA , saveDir , saveStr1 + "_clustered" , "PCA" , cluster = True)

    return list(kmeansMAST.labels_ )
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
    for i in range (len(xList)):
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
if __name__ == "__main__":
#input an xlsx dataframe and depending on what features you want to calculate, export .csv of features, performs feature elimination, maps features, and creates clusters 
    mainDir = str(sys.argv[1]) #Dataframe of Smiles and yields and ID numbers, and citations 
    rdkitAllow = int(sys.argv[2])
    morgAllow = int(sys.argv[3])
    maccAllow = int(sys.argv[4])
    scale1 = float(sys.argv[5]) #scaling for minClusterSize in PCA, take into account X sample size 
    chemStr = str(sys.argv[6])
    inputDF = pd.read_excel(mainDir + "/MasterDataFrame" + str(chemStr) + ".xlsx")
    featureTypes = []   #defines options for features you want 
    if rdkitAllow == 1:
        rdkitSave = "rdkit"
        featureTypes.append(rdkitSave)
    if morgAllow == 1:
        moldrSave = "morgan"
        featureTypes.append(moldrSave)
    if maccAllow == 1:
        maccSave = "MACCSkeys"
        featureTypes.append(maccSave)

    if not os.path.exists(mainDir + "/features"):
        os.makedirs(mainDir + "/features")
    featureMAST = pd.DataFrame()
    featureStr = chemStr + "_"
    for feature in featureTypes:
        #creates necessary features of interest, and generates plots 
        if feature == 'rdkit':
            featureDF = calc_rdkit_desc(smiles = list(inputDF['SMILES']) , multiprocess_ = True)
            createCSV(featureDF , mainDir + "/features", "initialX" + feature + "_" + chemStr  )
            featureMAST = pd.concat([featureMAST, featureDF], axis=1)
            featureStr += feature
        if feature == 'morgan':
            featureDF = calc_morgan_keys(smiles = list(inputDF['SMILES']) , multiprocess_ = True)
            createCSV(featureDF ,mainDir + "/features", "initialX" + feature + "_" + chemStr )
            featureMAST = pd.concat([featureMAST, featureDF], axis=1)
            featureStr += feature
        if feature == 'MACCSkeys':
            featureDF = calc_maccs_keys(smiles = list(inputDF['SMILES']) , multiprocess_ = True)
            createCSV(featureDF ,mainDir + "/features", "initialX" + feature + "_" + chemStr  )
            featureMAST = pd.concat([featureMAST, featureDF], axis=1)
            featureStr += feature

    featureMAST['SMILES'] = list(inputDF['SMILES']) 
    nanDict = locateNans(featureMAST)
    if len(nanDict) != 0:
        featureMAST = eliminateNans(featureMAST , nanDict)
    smilesMain = list(featureMAST['SMILES'])
    featureMAST = featureMAST.drop('SMILES', axis=1)
    featureMAST, featureList = featureFiltering(mainDir + "/features", featureMAST, list(featureMAST.columns) , featureStr)
    featureMAST = createCSV(featureMAST ,mainDir + "/features", featureStr + "_filteredFeatures_" + chemStr )
    kMeanLabels = DimensionalityReduction(featureMAST,  scale1 , mainDir + "/features" , featureStr  )
        
    finalDF = inputDF.loc[inputDF["SMILES"].isin(smilesMain)]
    finalDF['KMeanCluster'] = kMeanLabels

    finaldf = createCSV(finalDF , mainDir + "/features" , "masterDF_KMeansClustered_" + str(featureStr) )

    
