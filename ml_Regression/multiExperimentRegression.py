#Daniel Ruiz de Castilla
#06.16.2026
#Learns on a reaction data set to classify a substrate based on reaction conditions and predicts the metric of interest 
import os 
import sys
import glob
import pandas as pd 
import numpy as np
import chemdraw
import re
from pathlib import Path
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from DFTWorkflow.pitchingATent import featureFiltering
from dimensionalityReduction.reactivityFeatures import boxGen
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from ml_Regression.ML_RegressionV2 import gradientBoostRegression
def encodeReactionInfo(df):
    cols = list(df.columns)
    boxCol = boxGen(cols)
    encodeLst = listInputs(f"Enter a list of indices for string columns you want to encode as features\n {boxCol}\n")
    dropColsStr = [cols[int(i)] for i in encodeLst]
    for col in dropColsStr:
        types = {colStr : i for i , colStr in enumerate(df[col].unique())}
        df[col] = df[col].map(types)
    return df 
def main(dfMAST , saveDir , saveName , hyperFile):
    cols = list(dfMAST.columns)
    boxCols = boxGen(cols)

    dropIdx = listInputs(
        "Input the indices of all the columns you want to remove before XGB regression\n"
        f"{boxCols}\n"
    )
    dropIdx = [int(idx) for idx in dropIdx]
    dropCols = [cols[idx] for idx in dropIdx]

    dfMAST = dfMAST.drop(columns=dropCols)
    yStr = input("Enter the column for regression modeling: ")
    yCol = dfMAST[yStr]

    dfMAST = dfMAST.drop(columns = [yStr])

    X = encodeReactionInfo(dfMAST)

    reducedX , featureLabels  = featureFiltering(saveDir, X ,list(X.columns), saveName)

    j = 0
    while j < 4:

        newOutputDir = Path(saveDir) / f"trial{j}"
        newOutputDir.mkdir(parents=True, exist_ok=True)
        gradientBoostRegression(reducedX, yCol, hyperFile, newOutputDir , j)
        j +=1
if __name__ == "__main__":
    dfStr = str(sys.argv[1])
    dfMain = pd.read_csv(dfStr)
    saveDir = Path(sys.argv[2])
    saveStr = str(sys.argv[3])

    main(dfMain , saveDir , saveStr , saveDir / "hyperParameterOptimization.dat")




