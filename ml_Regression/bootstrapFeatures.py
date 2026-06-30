#Danny Ruiz de Castilla 06.29.2026
#brute forces the right feature combinations to use for an xgboost model when youre in the low data regime 
import sys
import glob
import numpy as np
import os
import pandas as pd
from pathlib import Path
import json
import random
import plotly.graph_objects as go
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split, GridSearchCV, StratifiedKFold
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from DFTWorkflow.pitchingATent import featureFiltering
from dimensionalityReduction.reactivityFeatures import boxGen
from ml_Regression.multiExperimentRegression import listInputs , encodeReactionInfo
def gbRegressorQuick(X , y , j , nEstimators, maxDepth , learningRate, maxFeatures ):
    X_train_CV, X_test, y_train_CV, y_test = train_test_split(X, y, test_size=0.2, random_state=j)
    gbFinal = GradientBoostingRegressor(n_estimators = nEstimators , max_depth = maxDepth , 
                                        learning_rate = learningRate , max_features = maxFeatures , random_state=j)
    gbFinal.fit(X_train_CV, y_train_CV)
    yTrain = gbFinal.predict(X_train_CV)
    yPred = gbFinal.predict(X_test)
    rmseTest = float(np.sqrt(mean_squared_error(y_test, yPred)))
    r2Test = float(r2_score(y_test, yPred))
    return rmseTest , r2Test
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
        gradientBoostRegression(reducedX, yCol, hyperFile, newOutputDir , j , True)
        j +=1
if __name__ == "__main__":
    dfStr = str(sys.argv[1])
    dfMain = pd.read_csv(dfStr)
    saveDir = Path(sys.argv[2])
    saveStr = str(sys.argv[3])

    main(dfMain , saveDir , saveStr , saveDir / "hyperParameterOptimization.dat")