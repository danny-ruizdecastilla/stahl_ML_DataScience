#Danny Ruiz de Castilla 06.29.2026
#brute forces the right feature combinations to use for an xgboost model when youre in the low data regime 
import sys
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
import shap
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from DFTWorkflow.pitchingATent import featureFiltering
from dimensionalityReduction.reactivityFeatures import boxGen
from ml_Regression.multiExperimentRegression import listInputs , encodeReactionInfo
def greedyBackwardElimination(reducedX, yCol, worstFeats, rmseToBeat, r2ToBeat,
                                nEstimators, maxDepth, learningRate, maxFeatures,
                                nSeeds=5, droppedSoFar=None):

    if droppedSoFar is None:
        droppedSoFar = []

    sortedWorst = sorted(
        ((k, v) for k, v in worstFeats.items() if v >= 3),
        key=lambda kv: kv[1], reverse=True
    )

    candidates = {}

    for feature, _ in sortedWorst:
        newX = reducedX.drop(columns=[feature])

        rmseF = []
        r2F = []
        for j in range(nSeeds):
            rmseTest, r2Test, shapDF = gbRegressorQuick(
                newX, yCol, j, nEstimators, maxDepth, learningRate, maxFeatures
            )
            rmseF.append(rmseTest)
            r2F.append(r2Test)

        rmseComp = np.mean(rmseF)
        r2Comp = np.mean(r2F)

        if rmseComp < rmseToBeat and r2Comp > r2ToBeat:
            candidates[feature] = [
                np.abs(rmseToBeat - rmseComp),
                np.abs(r2ToBeat - r2Comp),
                rmseComp,
                r2Comp
            ]

    # base case: nothing improved this round, stop recursion
    if not candidates:
        return reducedX, droppedSoFar, rmseToBeat, r2ToBeat

    # pick the single best feature to drop this round
    # ranked by rmse improvement magnitude (primary criterion)
    bestFeature = max(candidates.items(), key=lambda kv: kv[1][0])[0]
    bestRmse = candidates[bestFeature][2]
    bestR2 = candidates[bestFeature][3]

    newReducedX = reducedX.drop(columns=[bestFeature])
    droppedSoFar = droppedSoFar + [bestFeature]

    # remove the dropped feature from worstFeats so it's not retested
    newWorstFeats = {k: v for k, v in worstFeats.items() if k != bestFeature}

    print(f"Dropped: {bestFeature} | new RMSE: {bestRmse:.5f} (was {rmseToBeat:.5f}) "
          f"| new R2: {bestR2:.5f} (was {r2ToBeat:.5f})")

    # recurse with updated baseline
    return greedyBackwardElimination(
        newReducedX, yCol, newWorstFeats, bestRmse, bestR2,
        nEstimators, maxDepth, learningRate, maxFeatures,
        nSeeds=nSeeds, droppedSoFar=droppedSoFar
    )
def gbRegressorQuick(X , y , j , nEstimators, maxDepth , learningRate, maxFeatures ):
    X_train_CV, X_test, y_train_CV, y_test = train_test_split(X, y, test_size=0.2, random_state=j)
    gbFinal = GradientBoostingRegressor(n_estimators = nEstimators , max_depth = maxDepth , 
                                        learning_rate = learningRate , max_features = maxFeatures , random_state=j)
    gbFinal.fit(X_train_CV, y_train_CV)
    yTrain = gbFinal.predict(X_train_CV)
    yPred = gbFinal.predict(X_test)
    rmseTest = float(np.sqrt(mean_squared_error(y_test, yPred)))
    r2Test = float(r2_score(y_test, yPred))
    gbExplainer = shap.TreeExplainer(gbFinal , X_test)
    shapValues = gbExplainer(X_test)
    meanAbsShap = np.abs(shapValues.values).mean(axis=0)
    shapDF = pd.DataFrame({"feature" : X_test.columns , "Shap" : meanAbsShap}).sort_values("Shap" , ascending = False).reset_index(drop=False)
    return rmseTest , r2Test , shapDF
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
    targetCols = int(len(list(reducedX.columns))*0.20)

    j = 0
    rmse = []
    r2 = []
    nEstimators = int(input("Enter the number of estimators for this model: "))
    maxDepth = int(input("Enter the max Depth for this model: "))
    learningRate = float(input("Enter the learning rate for this model: "))
    maxFeatures = float(input("Enter the max Features for this model: "))
    worstFeats = {}
    while j < 5:
        rmseTest , r2Test , shapDF = gbRegressorQuick(reducedX , yCol , j, nEstimators , maxDepth, learningRate , maxFeatures )
        rmse.append(rmseTest)
        r2.append(r2Test)
        features = shapDF["feature"]
        nWorst = features[-targetCols:]
        for feature in nWorst:
            if feature in worstFeats.keys():
                worstFeats[feature] += 1 
            else:
                worstFeats[feature] = 1 
        
        j +=1
    
    rmseToBeat = np.mean(rmse)
    r2ToBeat = np.mean(r2)
    greedyBackwardElimination(reducedX, yCol, worstFeats, rmseToBeat, r2ToBeat,
                                    nEstimators, maxDepth, learningRate, maxFeatures)
    


if __name__ == "__main__":
    dfStr = str(sys.argv[1])
    dfMain = pd.read_csv(dfStr)
    saveDir = Path(sys.argv[2])
    saveStr = str(sys.argv[3])

    main(dfMain , saveDir , saveStr , saveDir / "hyperParameterOptimization.dat")