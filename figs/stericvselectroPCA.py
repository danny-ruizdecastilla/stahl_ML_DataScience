import pandas as pd
import os 
import sys
import glob
import re
import chemdraw
import base64
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import plotly
import plotly.io as pio
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from figs.featurePlotter import getFeaturePairs , standardCols
from figs.chemPlotlyV1 import insertIntoDataframe
from figs.chemPlotlyV2 import plotly_template , interactiveFigGenerator
from DFTWorkflow.pitchingATent import compressData , locateNans , eliminateNans , convertCanonical , featureFiltering
#Danny Ruiz de Castilla 05.08.
#Splits features into sterics vs electronics, reduces to PC1 and plots 
def modularPCA(X, numDim: int, outputDir, saveStr):
    np.random.seed(42)
    scaler = StandardScaler()
    scaledX = scaler.fit_transform(X)
    componentFull = min(X.shape[0], X.shape[1])
    pca = PCA(n_components=componentFull, svd_solver="full")
    xPCAFull = pca.fit_transform(scaledX)
    explainedVar = pca.explained_variance_ratio_
    print("Explained variance ratio:", explainedVar)

    explainedVarFile = os.path.join(outputDir, f"{saveStr}_explainedVarr.dat")
    if not os.path.exists(explainedVarFile):
        with open(explainedVarFile, "w") as file:
            for i, ev in enumerate(explainedVar, start=1):
                file.write(f"explained variance ratio: PC{i} {ev:.6f}\n")
    loadings = pd.DataFrame(
        pca.components_.T * np.sqrt(pca.explained_variance_),
        columns=[f'PC{i+1}' for i in range(componentFull)],
        index=X.columns
    )

    dfPCA = pd.DataFrame(
        xPCAFull[:, :numDim],
        columns=[f"PCA{i+1}" for i in range(numDim)]
    )

    featureExpDir = os.path.join(outputDir, "PC_featureExp")
    os.makedirs(featureExpDir, exist_ok=True)

    for j in range(componentFull):
        pcStr = f"PC{j+1}"
        pcFeatures = loadings[pcStr].abs().sort_values(ascending=False)
        filePath = os.path.join(featureExpDir, f"{pcStr}_{saveStr}_loadings.dat")
        with open(filePath, "w") as f:
            for key, value in pcFeatures.items():
                f.write(f"{key}\t{value:.6f}\n")

    return dfPCA

def pcafeatureSplitter(df , axisMotifs:dict , numDim: int , outputDir ):
    for axis in list(axisMotifs.keys()):
        dfMin = pd.DataFrame()
        axisStrings = list(axisMotifs[axis])
        for col in df.columns:
            if any(sub in col for sub in axisStrings):
                dfMin[col] = df[col]   
        dfAxis = modularPCA(dfMin , numDim , outputDir , axis)
        for col in dfAxis.columns:
            newCol = str(axis) + "_" + str(col)
            df[newCol] = dfAxis[col].copy()
        axisMotifs.setdefault(axis, []).append(newCol)
    return df , axisMotifs
def partitionChemistriesDF(chemistriesDir , dfMAST , partition, chemistry):
    substrateSpaces = glob.glob(chemistriesDir+ "/*.csv")
    substrateDict = {}
    for dir in substrateSpaces:
        #print(dir)
        df = pd.read_csv(dir)
        df = standardCols(df)
        smileList_ = list(df["SMILES"])
        canonicalSMILES_ = []
        for smile in smileList_:
            canonical = convertCanonical(smile)
            canonicalSMILES_.append(canonical)
        df["canonicalSMILES"] = canonicalSMILES_
        substrateDF = (dfMAST[dfMAST["canonicalSMILES"].isin(canonicalSMILES_)].drop_duplicates(subset="canonicalSMILES").reset_index(drop=True)) 
        substrateDF = insertIntoDataframe(substrateDF , "canonicalSMILES" , "canonicalSMILES" , df , partition , 2) #  df2 , inputStr , indexInterest  )
        split1 = dir.split("/")[-1]
        split2 = split1.split(chemistry)[0]
        substrateDict[split2] = substrateDF
    return substrateDict

def main(substrateSpace , substrateData , chemistry , outputDir , axisMotifs, elimPhrases , partitionStr ):
    partitionPrompt = "Please Enter a list of " + str(partitionStr) + " values to partition the dataset. (Ex: 40,50,60,70,80,90): "
    reactivityPartitions = listInputs(partitionPrompt)
    reactivityInts = [int(x) for x in reactivityPartitions]
    initdataSets = glob.glob(substrateData + "/*.csv")
    initdataSets = sorted(initdataSets)
    Xdataframe , smileList  , yieldList_= compressData(initdataSets , partitionStr , elimPhrases , outputDir , "rawUnprocessed")
    Xdataframe = standardCols(Xdataframe)
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
    
    X , featureLabels  = featureFiltering(outputDir, Xdataframe ,featureLabels , chemistry)

    axisDF , axisMotifs = pcafeatureSplitter(X , axisMotifs , 1 , outputDir)

    axisLabels = [list(axisMotifs["sterics"] ), list(axisMotifs["electronics"])]
    featureList = list(axisDF.columns)
    featurePairs = getFeaturePairs(featureList , axisLabels)
    axisDF["canonicalSMILES"] = canonicalSMILES

    substrateDict = partitionChemistriesDF(substrateSpace , axisDF , partitionStr, chemistry)
    figList = []
    for pair in featurePairs:
        xStr = str(pair[0])
        yStr = str(pair[1])
        if "PC" in xStr or "PC" in yStr:
            for chemistry_, chemistryDF in substrateDict.items():
                greyDFs = []
                for key, df in substrateDict.items():
                    if key != chemistry_:
                        df = df.copy()
                        greyDFs.append(df)
                greyDF = pd.concat(greyDFs, axis=0).reset_index(drop=True)
                for partition in reactivityInts:
                    titleStr = xStr + " vs " + yStr + " for " + str(chemistry_) + " at " + str(partition) + "% yield" 
                    fig = interactiveFigGenerator(chemistryDF , greyDF , partition , xStr , yStr, str(chemistry_) , "Background Substrates" ,  titleStr)
                    figList.append(fig)
    html_parts = [pio.to_html(fig, include_plotlyjs=False, full_html=False) for fig in figList]
    html_output = """
    <html>
    <head>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    </head>
    <body>
    {}
    </body>
    </html>
    """.format('\n<hr>\n'.join(html_parts))  # <hr> separates plots

    # Save to file
    with open(outputDir + "/combinations.html", "w") as f:
        f.write(html_output)
    
if __name__ == "__main__":
    datasetDir = str(sys.argv[1])
    chemistriesDir = str(sys.argv[2])
    figDir = str(sys.argv[3])
    if not os.path.exists(figDir):
        os.makedirs(figDir)
    axisMotifs = {"sterics" : ["distance" , "Buried" , "angle" , "dihedral" , "%Vbur"] , "electronics" : ["fukui" , "μ" , "ω" , "Dipole" , "NBO" , "polar" , "HOMO" , "NMR" , "η" ]}
    elimFile = str(sys.argv[4])
    chemistry = str(sys.argv[5])
    partitionStr = str(sys.argv[6])
    if os.path.exists(elimFile):
        with open(elimFile, 'r') as file:
            content = file.read()
            eliminatedPhrases = [item.strip() for item in content.split(',') if item.strip()]
    else: 
        eliminatedPhrases = ["SMILES" , "Compound_Name", "Yield", "ChemistryType"  ]
    main(chemistriesDir , datasetDir,chemistry ,  figDir , axisMotifs, eliminatedPhrases , partitionStr)



