import sys
import glob
import numpy as np
import os
#from rdkit import Chem
import pandas as pd
#from rdkit.Chem.PandasTools import LoadSDF
import matplotlib.pyplot as plt
import random
import umap
from itertools import combinations
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from DFTWorkflow.featureMaping import savePNG , createCSV
from DFTWorkflow.expt2_feature_filtering import *
#Danny Ruiz de Castilla 02.28.2025

#pcaDict = {canonicalSmiles : [ yield]} for local Chemistries 
#smilesDict = {canonicalSmiles : [ xPCA , yPCA]}
def convertCanonical(str):
    from rdkit import Chem
    mol = Chem.MolFromSmiles(str)
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    return canonical
def dimensionalityReduction(X , smiles, reduceStr , saveDir , chemType , specialK):
    if "PCA" in reduceStr:
        np.random.seed(42)
        scaler = StandardScaler()
        scaledX = scaler.fit_transform(X)


        pca = PCA(n_components = X.shape[1] ,svd_solver="full" )
        xPCAFullRank = pca.fit_transform(scaledX)
        if "_" in reduceStr:
            reducerType = reduceStr.split("_")[1]
        else:
            reducerType = "PCA"
        if reducerType == 'umap':
            reducer = umap.UMAP(n_components=2, metric='euclidean', random_state=0)
            X_umap = reducer.fit_transform(xPCAFullRank)


        explainedVar = pca.explained_variance_ratio_
        print(explainedVar)
        if not os.path.exists(saveDir + "/" + "explainedVarr.dat"):
            with open(saveDir + "/" + "explainedVarr.dat", "w") as file:
                for i in range (len(explainedVar)):
                    file.write(f"explained variance ratio: PC {i + 1} {explainedVar[i]:.6f}\n")

        top2 = np.argsort(explainedVar)[-2:][::-1] 
        #print(top2)
        top2PCA = pca.components_[top2]
        xPCA_ = scaledX @ top2PCA.T + 11
        #print(xPCA_)
        loadings = pd.DataFrame(pca.components_.T * np.sqrt(pca.explained_variance_), columns=[f'PC{i+1}' for i in range(pca.n_components_)] , index = X.columns)

        top_features_pc1 = loadings['PC1'].abs().sort_values(ascending=False)
        top_features_pc2 = loadings['PC2'].abs().sort_values(ascending=False)

        topFeatures = pd.DataFrame({
            "Feature": top_features_pc1.index,
            "PC1_Contribution": top_features_pc1.values,
            "PC2_Contribution": top_features_pc2.loc[top_features_pc1.index].values  # Align PC2 with PC1 sorting
        })

        filePath = "topFeatures.csv"
        if not os.path.exists(saveDir + "/" + filePath):
            topFeatures.to_csv(saveDir + "/" + filePath, sep="\t", index=False)
        dfPCA = pd.DataFrame(xPCA_, columns = ["PCA1" , "PCA2"] )
        #print(dfPCA)
        xAxis = list(dfPCA["PCA1"])  #scaling factor to avoid negative numbers in PCA
        yAxis = list(dfPCA["PCA2"])
        dfPCA["SMILES"] = smiles
        kmeansMAST = KMeans( n_clusters = specialK , random_state = 42)
        kmeansMAST.fit(xPCA_)
        dfPCA['Clusters'] = kmeansMAST.labels_
        savePNG(dfPCA , saveDir , chemType + "_clustered" , "PCA" , cluster = True)
        X["SMILES"] = smiles
        X["PCA1"] = xAxis 
        X["PCA2"] = yAxis
        if not os.path.exists(saveDir + "/" + "featureDistribution" ):
            os.makedirs(saveDir + "/" + "featureDistribution")
        createCSV(X, saveDir + "/" + "featureDistribution" , chemType+ "_PCADistributed_Features")

        pcaDict = {}
        for i in range (len(smiles)):
            smile = smiles[i]
            pcaDict[smile] = [xAxis[i] , yAxis[i]]

        try:
            dfUMAP = pd.DataFrame(X_umap, columns=["UMAP1", "UMAP2"])
            umap1 = list(dfUMAP['UMAP1'])
            umap2 = list(dfUMAP['UMAP2'])
            pcaUmapdict = {}
            X["UMAP1"] = umap1 
            X["UMAP2"] = umap2
            createCSV(X, saveDir + "/" + "featureDistribution" , chemType+ "_UMAPDistributed_Features")
            for i in range (len(smiles)):
                smile = smiles[i]
                pcaUmapdict[smile] = [umap1[i] , umap2[i]]
            return pcaUmapdict 
        except NameError:
            return pcaDict  
        
    elif "UMAP" == reduceStr:
        reducer = umap.UMAP(n_components=2, metric='euclidean', random_state=0)
        X_umap = reducer.fit_transform(X)
        dfUMAP = pd.DataFrame(X_umap, columns=["UMAP1", "UMAP2"])
        umap1 = list(dfUMAP['UMAP1'])
        umap2 = list(dfUMAP['UMAP2'])
        umapDict = {}
        X["UMAP1"] = umap1 
        X["UMAP2"] = umap2
        createCSV(X, saveDir + "/" + "featureDistribution" , chemType+ "_UMAPDistributed_Features")
        for i in range (len(smiles)):
            smile = smiles[i]
            umapDict[smile] = [umap1[i] , umap2[i]]
        return umapDict 

def findHighlights(smilesDict, partition1, chemistryDict):
    highlightDict = {}
    hollowDict = {}
    to_remove = []

    for substrate, parameter in chemistryDict.items():
        #print(substrate , parameter)
        if substrate in smilesDict and parameter >= partition1:
            #print("MATCH")
            xPCA, yPCA = smilesDict[substrate]  # Direct unpacking
            highlightDict[substrate] = [xPCA, yPCA , parameter]
            to_remove.append(substrate)
        
        elif substrate in smilesDict and parameter < partition1: #hollow points
            hollowX , hollowY = smilesDict[substrate]
            hollowDict[substrate] = [hollowX , hollowY , parameter]
            to_remove.append(substrate)
    

    # Remove keys after iteration to avoid modifying the dictionary while looping
    for key in to_remove:
        del smilesDict[key]

    return highlightDict, hollowDict , smilesDict


def makePlots(pcaDict ,  partitionList , chemistryDicts  ,chemistryStr , colors , partitionStr , reduceType): 
    xLabel = str(reduceType) + "1"
    yLabel = str(reduceType) + "2"
    for partition in partitionList:
        #print("partition" , partition)
        for i , chemistryDict in enumerate(chemistryDicts): #this is the highlighted chemistry in the plot
            chemistryLabel = chemistryStr[i].split(".")[0]
            if not os.path.exists(reduceDir + "/" + chemistryLabel ):
                os.makedirs(reduceDir  + "/" + chemistryLabel)
            highlightDict , hollowDict , smilesDict = findHighlights(pcaDict.copy()  ,  partition , chemistryDict )
            if not os.path.exists(reduceDir + "/" + chemistryLabel  + "/" +  str(chemistryLabel)  + str(reduceType) + "_Coordinates.dat"):
                with open(reduceDir + "/" + chemistryLabel  + "/" +  str(chemistryLabel)  + str(reduceType) + "_Coordinates.dat" , "w") as file:
                    file.write(f"SMILES,{xLabel},{yLabel},Yield\n")
                with open(reduceDir + "/" + chemistryLabel  + "/" +  str(chemistryLabel)  + str(reduceType) + "_Coordinates.dat" , "a") as file:
                    for smile_ in list(highlightDict.keys()):
                        file.write(f"{smile_},{highlightDict[smile_][0]},{highlightDict[smile_][1]},{highlightDict[smile_][2]}\n")   
                    for smile_ in list(hollowDict.keys()):
                        file.write(f"{smile_},{hollowDict[smile_][0]},{hollowDict[smile_][1]},{hollowDict[smile_][2]}\n")  
            coordinatesHighlighted = [v[:2] for v in list(highlightDict.values())]
            coordinatesHollowed = [v[:2] for v in list(hollowDict.values())]  
            if len(highlightDict) != 0:
                xHigh, yHigh = zip(*coordinatesHighlighted)
                colorScatter = True
            else:
                colorScatter = False

            if len(hollowDict) != 0:
                xHollow , yHollow = zip(*coordinatesHollowed)
                hollowScatter = True
            else:
                hollowScatter = False                

            xBland , yBland = zip(*smilesDict.values())
            if not os.path.exists(reduceDir  + "/" + chemistryLabel  + "/" +  str(chemistryLabel)  + "GreyedOutSubstrates.dat"):
                with open(reduceDir + "/" + chemistryLabel  + "/" +  str(chemistryLabel)  + "GreyedOutSubstrates.dat" , "w") as file:
                    file.write(f"SMILES,{xLabel},{yLabel},Yield\n")
                with open(reduceDir + "/" + chemistryLabel  + "/" +  str(chemistryLabel)  + "GreyedOutSubstrates.dat" , "a") as file:
                    for smile in list(smilesDict.keys()):   
                        file.write(f"{smile},{smilesDict[smile][0]},{smilesDict[smile][1]}\n")   
            dpi = 300
            plt.figure(figsize=(800 / dpi, 600 / dpi), dpi=dpi)  # 800x600 pixels
            plt.scatter(list(xBland), list(yBland), c="grey", alpha=0.14 , s=10)
            if colorScatter:
                plt.scatter(list(xHigh), list(yHigh), facecolors = colors[i], edgecolors = 'none', alpha=0.35 , s=10)
            if hollowScatter:
                plt.scatter(list(xHollow), list(yHollow), facecolors = 'none', edgecolors = colors[i] ,  alpha=0.35 , s=10)
            plt.xlabel(xLabel, fontsize=11,  color='black')
            plt.ylabel(yLabel, fontsize=11,  color='black')
            plt.title("Substrate Scope for " + str(chemistryLabel) + " at " + str(partition) + "%" + str(partitionStr), fontsize=18, fontweight='bold', color='navy')
            plt.xticks(fontsize=7, color='black')
            plt.yticks(fontsize=7, color='black')
            #plt.grid(True, linestyle='--', alpha=0.5)
            plt.savefig(reduceDir + "/" + chemistryLabel + "/" + str(chemistryLabel) + "at" + str(partition) + ".png", dpi=300, bbox_inches='tight')
            plt.close()
    
def featureFiltering(outDir , X , feature_labels , featureStr):
    if not os.path.exists(outDir + "/" + str(featureStr) + "featureFiltering.dat"):
        with open(outDir + "/" + str(featureStr) + "featureFiltering.dat", "w") as f:
            f.write(f"Total starting feature count: {len(feature_labels)}")
            f.write("".join([f'\n\t{label}' for label in feature_labels]))
            
            X, feature_labels, dropped_features = remove_by_variance(X, feature_labels)
            #print("90" , type(X))

            text = "\n\n\nFeatures drop due to low variance: " + "".join([f'\n\t{label}' for label in dropped_features])
            f.write(text)

            X, feature_labels, drop_group = correlation_analysis(X, feature_labels, threshold=0.95)
            text = "\n\n\nFeatures drop due correlation: (STILL HAS ISSUES) "
            #print("97" , type(X))
            import json
            text += json.dumps(drop_group, indent=4).replace('\n', '\n\t')
            f.write(text)

            X, feature_labels  = spearmanr_correlation(X, feature_labels, threshold=0.95)
            #print("103" , type(X))
    return X , feature_labels
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
def eliminateNans(df, nanDict):
    allNanRows = set()
    for key in nanDict.keys():
        rows = nanDict[key]
        allNanRows.update(rows)
    
    df = df.drop(list(allNanRows))
    df = df.reset_index(drop=True)
    return df
def compressData(dataframeDirs , regressionStr , usualSuspects):
    

    dataframeMast = pd.DataFrame()
    
    # Check if input is valid
    if not dataframeDirs:
        return pd.DataFrame(), pd.Series()
    
    try:
        # Process each file
        for dfDir in dataframeDirs:
            
            try:
                df = pd.read_csv(dfDir)
                
                # Handle empty files
                if df.empty:
                    print(f"Warning: Empty file {dfDir}")
                    continue
                    
                if dataframeMast.empty:
                    dataframeMast = df.copy()
                else:
                    commonCols = list(set(df.columns) & set(dataframeMast.columns))
                    if not commonCols:
                        print(f"Warning: No common columns with file {dfDir}")
                        continue
                    dataframeMast = pd.concat([dataframeMast[commonCols], df[commonCols]], ignore_index=True)
            
            except Exception as e:
                print(f"Error processing file {dfDir}: {e}")
        
        # Check if master dataframe is populated
        if dataframeMast.empty:
            return pd.DataFrame(), pd.Series()
        
        # Save SMILES column if it exists
        if 'SMILES' in dataframeMast.columns:
            smileList = dataframeMast['SMILES'].copy()
        # Save Yield column if it exists
        if regressionStr in dataframeMast.columns:
            yieldList = dataframeMast[regressionStr].copy()

        else:
            print("Warning: SMILES column not found in the dataframe")
            df = pd.read_csv(dfDir)
            print(df.columns)
            smileList = pd.Series()
        elimCol = [col for col in dataframeMast.columns if any(frag in col for frag in usualSuspects)]
        dataframeMast = dataframeMast.drop(columns=elimCol)

        return dataframeMast, smileList , yieldList
        
    except Exception as e:
        print(f"Critical error in transformations: {e}")
        return pd.DataFrame(), pd.Series()
def partitionChemistries(substrateSpace , chemistry):
    chemDirList = []
    chemSpaceDicts = []
    substrateSpaces = glob.glob(substrateSpace+ "/*.csv")
    for dir in substrateSpaces:
        #print(dir)
        df = pd.read_csv(dir)
        smileList_ = list(df["SMILES"])
        canonicalSMILES_ = []
        for smile in smileList_:
            canonical = convertCanonical(smile)
            canonicalSMILES_.append(canonical)
        yieldList = list(df["Yield"])
        chemistryDict = dict(zip(canonicalSMILES_, yieldList))
        chemSpaceDicts.append(chemistryDict)
        split1 = dir.split("/")[-1]

        split2 = split1.split(chemistry)[0]
        chemDirList.append(split2)
    return chemSpaceDicts , chemDirList 
def main(substrateSpace , substrateData , chemistry , clusterK ,  elimFile ,outputDir , outputStr ):
    if not os.path.exists(outputDir): 
        os.makedirs(outputDir)
    if os.path.exists(elimFile):
        with open(elimFile, 'r') as file:
            content = file.read()
            eliminatedPhrases = [item.strip() for item in content.split(',') if item.strip()]
    else: 
        eliminatedPhrases = ["SMILES" , "Compound_Name", "Yield", "ChemistryType"  ]
    partitionList = [50 , 70 , 85 , 90]
    colorList = ['red' , 'blue' , 'green' , 'yellow']
    initdataSets = glob.glob(substrateData + "/*.csv")
    initdataSets = sorted(initdataSets)
    Xdataframe , smileList  , yieldList_= compressData(initdataSets , "Yield" , eliminatedPhrases)

    nanDict = locateNans(Xdataframe)
    if len(nanDict) != 0:
        Xdataframe["SMILES"] = smileList
        Xdataframe = eliminateNans(Xdataframe , nanDict)
    smileList = Xdataframe["SMILES"].copy()
    canonicalSMILES = []
    for smile in smileList:
        canonical = convertCanonical(smile)
        canonicalSMILES.append(canonical)
    Xdataframe = Xdataframe.drop("SMILES", axis=1)
    featureLabels = list(Xdataframe.columns)
    
    X , featureLabels  = featureFiltering(reduceDir, Xdataframe ,featureLabels , chemistry)
    dataDict  = dimensionalityReduction(X , canonicalSMILES , outputStr , reduceDir , chemistry , clusterK)
    smilesTot = list(dataDict.keys())
    coordinateTot = list(dataDict.values())
    for i in range(len(smilesTot)):
    
        smiles = smilesTot[i]
        xAxis = coordinateTot[i][0]
        yAxis = coordinateTot[i][1]
        with open(reduceDir + "/" + str(outputStr) + "Coordinates.dat", "a") as file:
            file.write(f"{smiles},{xAxis} ,  {yAxis}\n") 

    chemSpaceDicts , chemDirList = partitionChemistries(substrateSpace , chemistry)
    #print(chemDirList)
    makePlots(dataDict.copy() ,  partitionList , chemSpaceDicts  , chemDirList, colorList , "Yield" , reduceStr)
if __name__ == "__main__":

    chemistryDirs = str(sys.argv[1])
    datasetDir = str(sys.argv[2])
    chemType = str(sys.argv[3])
    specialK = int(sys.argv[4]) #Kmeans cluster
    saveDir = str(sys.argv[5])
    elimFile = str(sys.argv[6])
    reduceStr = str(sys.argv[7])
    reduceDir = saveDir + "/" + reduceStr
    main(chemistryDirs , datasetDir , chemType , specialK , elimFile , reduceDir , reduceStr)
    
