import pandas as pd
import glob
import numpy as np
import os 
import sys
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from figs.chemPlotlyV1 import dat2List
from DFTWorkflow.featureMaping import  createCSV
def getDeltaPairs(atom1:str, atom2:str , featureList):
    atom1Features = []
    atom2Features = []
    for feature in featureList:
        if atom1 in feature and atom2 not in feature:
            print("atom1")
            feat = feature.replace(atom1, 'diffAtom')  
            atom1Features.append(feat)
        if atom2 in feature and atom1 not in feature:
            print("atom2")
            feat = feature.replace(atom2, 'diffAtom')  
            atom2Features.append(feat)
    eligibleFeats = set(atom1Features) & set(atom2Features) #only common features prevail
    pairs = []
    for feat in eligibleFeats:
        feat1 = feat.replace("diffAtom" , atom1)
        feat2 = feat.replace("diffAtom" , atom2)
        pairs.append([feat1 , feat2])
    return pairs

def abstractFeats(df , elimSet , atom1, atom2):
    unfilteredColumns = list(df.columns)
    filteredCols = [f for f in unfilteredColumns if not any(sub in f for sub in elimSet)]
    deltaPairs = getDeltaPairs(atom1 , atom2, list(filteredCols))
    
    print("deltaPairs" , deltaPairs)
    for pair in deltaPairs:
        if pair[0] in df.columns and pair[1] in df.columns:
            name = f"delta({pair[0]}_{pair[1]})"
            df[name] = df[pair[0]] - df[pair[1]]
            featMin = df[[pair[0] , pair[1]]].max(axis=1)
            featMax = df[[pair[0] , pair[1]]].min(axis=1)
            if atom1 in pair[0]:
                maxName = pair[0].split(atom1)[0] + "Cmax" + pair[0].split(atom1)[-1]
                minName = pair[0].split(atom1)[0] + "Cmin" + pair[0].split(atom1)[-1]
                df[maxName] = featMax
                df[minName] = featMin
            elif atom1 in pair[1]:
                maxName = pair[1].split(atom1)[0] + "Cmax" + pair[1].split(atom1)[-1]
                minName = pair[1].split(atom1)[0] + "Cmin" + pair[1].split(atom1)[-1]
                df[maxName] = featMax
                df[minName] = featMin               
        else:
            print(f"Columns {pair[0]} or {pair[1]} not found")
    
    return df
if __name__ == "__main__":
    featureDir = str(sys.argv[1])
    atom1 = str(sys.argv[2])
    atom2 = str(sys.argv[3])
    elimFile = str(sys.argv[4])
    outputDir = str(sys.argv[5])
    if os.path.exists(elimFile):
        elimPhrases = dat2List(elimFile , " , ")
    else: 
        elimPhrases = ["SMILES" , "Compound_Name", "Yield", "ChemistryType", "stdev" , "max"]
    elimSet = set(elimPhrases)
    dataframes = glob.glob(featureDir + "/*.csv")
    for dfDir in dataframes:
        df = pd.read_csv(dfDir , encoding='utf-8')
        dfWithDeltas = abstractFeats(df , elimSet , atom1 , atom2)
        createCSV(dfWithDeltas , outputDir , str(dfDir.split("/")[-1].split(".")[0]))