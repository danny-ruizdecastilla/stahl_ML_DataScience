import pandas as pd
import glob
import numpy as np
import os 
import sys
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from figs.chemPlotlyV1 import dat2List
from DFTWorkflow.featureMaping import  createCSV
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from ml_Regression.getDeltas import abstractFeats
def getMinMax(df , colList):
    featMin = df[colList].max(axis=1)
    featMax = df[colList].min(axis=1)
    featMin, featMax
def wildConsolidate(df , elimSet , c1c2Dict):
    unfilteredColumns = list(df.columns)
    filteredCols = [f for f in unfilteredColumns if not any(sub in f for sub in elimSet)]
    wildCols = [col for col in filteredCols if ("Wild" in col)]
    ends = [[col.split("_")[0] , col.split("_")[-1] ] for col in wildCols]
    ends =  [list(t) for t in set(tuple(sublist) for sublist in ends)]
    for end in ends:
        feats = [col for col in wildCols if (end[0] in col) and (end[1] in col)]
        min , max = getMinMax(df, feats)
        df[str(end[0] + "_min_Wild_" + str(end[1]))] = min 

        df[str(end[0] + "_max_Wild_" + str(end[1]))] = max
    return df    
if __name__ == "__main__":
    c1c2 = listInputs(f"Enter the string values for both Alkene atoms. Ex: C1,C2\n")
    c1c2Dict = {}
    for atom in c1c2:
        wildList = listInputs(f"Enter the 2 Wildcard atom string values for Alkene {atom}:\n")
        c1c2Dict[atom] = wildList
    
    featureDir = str(sys.argv[1])
    c1c2list = list(c1c2.keys())
    atom1 = c1c2list[0]
    atom2 = c1c2list[1]
    elimFile = str(sys.argv[2])
    outputDir = str(sys.argv[3])
    if os.path.exists(elimFile):
        elimPhrases = dat2List(elimFile , " , ")
    else: 
        elimPhrases = ["SMILES" , "Compound_Name", "Yield", "ChemistryType", "stdev" , "max"]
    elimSet = set(elimPhrases)
    dataframes = glob.glob(featureDir + "/*.csv")
    for dfDir in dataframes:
        df = pd.read_csv(dfDir , encoding='utf-8')
        dfWithDeltas = abstractFeats(df , elimSet , atom1 , atom2)
        dfSmearedWilds = wildConsolidate(df, elimSet, c1c2Dict)
        createCSV(dfSmearedWilds , outputDir , str(dfDir.split("/")[-1].split(".")[0]))