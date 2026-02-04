#Predicts the rank positioning on a 1-d axis based on a feature matrix
#import torch
#import torch.nn as nn
import os 
import sys
import glob
import pandas as pd 
import numpy as np
import chemdraw
import re
from pathlib import Path
from rdkit import Chem
from sklearn.ensemble import RandomForestRegressor , GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.neural_network import MLPRegressor
from sklearn.model_selection import KFold
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from dimensionalityReduction.reactivityFeatures import boxGen
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from ml_Regression.ML_RegressionV2 import GridSearchCV , stratifiedRegressionCV
from DFTWorkflow.pitchingATent import featureFiltering

def canonical_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None
def rFR_2Source(X_trainDF, y_trainDF,X_testDF, yStr , hyperParmFile, outputDir):
    import json 
    X_testDF = X_testDF[X_trainDF.columns]
    if hyperParmFile.exists():
        with open(hyperParmFile, "r") as f:
            savedParms = json.load(f)
        best_params = savedParms["best_params"]
        rfMAST = RandomForestRegressor(**best_params, random_state=1)
        rfFinal = RandomForestRegressor(**best_params, random_state=1)
    else:
        rfGrid = {
            'n_estimators': [100, 200, 300],
            'max_depth': [2, 4, 8, None],
            'max_features': ['sqrt', 0.5, None]
        }
        rf = RandomForestRegressor(random_state=1)
        gridSearch = GridSearchCV(rf, rfGrid, cv=5, scoring="neg_mean_squared_error", n_jobs=-1)
        gridSearch.fit(X_trainDF, y_trainDF)
        with open(hyperParmFile, "w") as f:
            json.dump({"best_score": gridSearch.best_score_, "best_params": gridSearch.best_params_}, f, indent=4)
        rfMAST = RandomForestRegressor(**gridSearch.best_params_, random_state=1)
        rfFinal = RandomForestRegressor(**gridSearch.best_params_, random_state=2)

    mseCV, r2CV = stratifiedRegressionCV(X_trainDF, y_trainDF, 42 , rfMAST, n_splits=5 )

    cvFile = Path(outputDir) / "randomForest" / "crossvalidation" / "scores.dat"
    cvFile.parent.mkdir(parents=True, exist_ok=True)
    with open(cvFile, "w") as f:
        for i in range(len(mseCV)):
            f.write(f"Fold {i} Results: mse {mseCV[i]} | r2 {r2CV[i]}\n")

    rfFinal.fit(X_trainDF, y_trainDF)
    yTrain = rfFinal.predict(X_trainDF)
    yPred = rfFinal.predict(X_testDF)

    rmseTrain = float(np.sqrt(mean_squared_error(y_trainDF , yTrain)))
    r2Train = float(r2_score(y_trainDF , yTrain))

    testFile = Path(outputDir) / "randomForest" / "Eval" / "testScores.dat"
    testFile.parent.mkdir(parents=True, exist_ok=True)
    trainFile = Path(outputDir) / "randomForest" / "Eval" / "trainScores.dat"
    with open(trainFile, "w") as f:
        f.write(f"Train Results: rmse {rmseTrain} | r2 {r2Train}\n")

    importances = rfFinal.feature_importances_
    if isinstance(X_trainDF, pd.DataFrame):
        feature_names = X_trainDF.columns
    else:
        feature_names = [f"f{i}" for i in range(importances.shape[0])]
    featDF = pd.DataFrame({"feature": feature_names, "importance": importances}).sort_values(by="importance", ascending=False)
    featFile = Path(outputDir) / "randomForest" / "Eval" / "testWeights.csv"
    featDF.to_csv(featFile, index=False)
    X_testDF[yStr] = yPred
    return X_testDF
def regressionPipeline(X_trainDF, y_trainDF, X_testDF, yStr, model_class, param_grid=None,
                       hyperParmFile=None, outputDir=None, PipeObj=None, model_name="model"):
    import json

    X_testDF = X_testDF[X_trainDF.columns]

    # Load hyperparameters if they exist
    if hyperParmFile and hyperParmFile.exists():
        with open(hyperParmFile, "r") as f:
            savedParms = json.load(f)
        best_params = savedParms["best_params"]
    else:
        if param_grid is None:
            param_grid = {}
        if PipeObj is None:
            mlObj = model_class(random_state = 42)
            gridSearch = GridSearchCV(mlObj, param_grid, cv=5, scoring="neg_mean_squared_error", n_jobs=-1, refit=True)
        else:
            gridSearch = GridSearchCV(PipeObj, param_grid, cv=5, scoring="neg_mean_squared_error", n_jobs=-1, refit=True)
        gridSearch.fit(X_trainDF, y_trainDF)
        best_params = gridSearch.best_params_
        if hyperParmFile:
            with open(hyperParmFile, "w") as f:
                json.dump({"best_score": gridSearch.best_score_, "best_params": best_params}, f, indent=4)

    modelMAST = model_class(**best_params, random_state=1)
    modelFinal = model_class(**best_params, random_state=2)

    mseCV, r2CV = stratifiedRegressionCV(X_trainDF, y_trainDF, 42, modelMAST, n_splits=5)

    cvFile = Path(outputDir) / model_name / "crossvalidation" / "scores.dat"
    cvFile.parent.mkdir(parents=True, exist_ok=True)
    with open(cvFile, "w") as f:
        for i, (mse, r2) in enumerate(zip(mseCV, r2CV)):
            f.write(f"Fold {i} Results: mse {mse} | r2 {r2}\n")

    modelFinal.fit(X_trainDF, y_trainDF)
    yTrain = modelFinal.predict(X_trainDF)
    yPred = modelFinal.predict(X_testDF)

    rmseTrain = float(np.sqrt(mean_squared_error(y_trainDF, yTrain)))
    r2Train = float(r2_score(y_trainDF, yTrain))

    testFile = Path(outputDir) / model_name / "Eval" / "testScores.dat"
    testFile.parent.mkdir(parents=True, exist_ok=True)
    trainFile = Path(outputDir) / model_name / "Eval" / "trainScores.dat"
    with open(trainFile, "w") as f:
        f.write(f"Train Results: rmse {rmseTrain} | r2 {r2Train}\n")

    importances = getattr(modelFinal, "feature_importances_", None) or getattr(modelFinal, "coef_", None)
    if importances is not None:
        if isinstance(X_trainDF, pd.DataFrame):
            feature_names = X_trainDF.columns
        else:
            feature_names = [f"f{i}" for i in range(importances.shape[0])]
        featDF = pd.DataFrame({"feature": feature_names, "importance": importances}).sort_values(by="importance", ascending=False)
        featFile = Path(outputDir) / model_name / "Eval" / "featureImportance.csv"
        featDF.to_csv(featFile, index=False)

    X_testDF[yStr] = yPred
    return X_testDF

def addToDF(df, colsHash, appendDF):
    for name, info in colsHash.items():
        df[name] = info

    eligibleCols = df.columns.intersection(appendDF.columns)

    newRows = appendDF.loc[:, eligibleCols]

    return pd.concat([df, newRows], ignore_index=True)
def main(training , testing , outputDir  , regressionType , featureStr ):
    trainingDF = pd.read_csv(training, encoding="utf-8" )
    testingDF = pd.read_csv(testing, encoding="utf-8")
    trainingCols = list(trainingDF.columns)
    testingCols = list(testingDF.columns)
    trainingBox = boxGen(trainingCols)
    testingBox = boxGen(testingCols)

    smilesIdxTraining = int(input(f"Enter the index of the smiles column:\n {trainingBox}\n"))
    smilesIdxTesting =  int(input(f"Enter the index of the smiles column:\n {testingBox}\n"))

    yIdxTraining = int(input(f"Enter the index of the Y column:\n {trainingBox}\n"))

    trainingYStr = trainingCols[yIdxTraining]
    yTraining = list(trainingDF[trainingYStr])

    trainingDF["canonical_smiles"] = [canonical_smiles(smi) for smi in trainingDF[trainingCols[smilesIdxTraining]]]
    trainingCanonicalSet = set(trainingDF["canonical_smiles"])

    testingDF["canonical_smiles"] = [canonical_smiles(smi) for smi in testingDF[testingCols[smilesIdxTesting]]]
    mask = testingDF["canonical_smiles"].isin(trainingCanonicalSet)

    rejectedDF = testingDF[mask].reset_index(drop=True)
    rejectedRanks = []
    for smiles in list(rejectedDF["canonical_smiles"]):
        rank = trainingDF.loc[trainingDF["canonical_smiles"] == smiles,trainingYStr]
        rejectedRanks.append(rank)
    rejectedDF[trainingYStr] = rejectedRanks

    testingDF  = testingDF[~mask].reset_index(drop=True)

    trainingDF = trainingDF.drop(columns=["canonical_smiles" , trainingYStr])

    dropCols = listInputs(f"Enter a list of indexes for the columns you want to drop from the test and training data\n {trainingBox}\n")
    dropColsStr = [trainingCols[int(i)] for i in dropCols]

    testingCols = list(testingDF.columns)
    testingBox = boxGen(testingCols)
    keepCols = listInputs(f"Enter a list of indexes for the columns you want to keep for the test data\n {testingBox}\n")
    keepHash = {testingCols[int(idx)] : list(testingDF[str(testingCols[int(idx)])]) for idx in keepCols}
    #print(dropColsStr)
    trainingDF = trainingDF.drop(columns = dropColsStr)
    testingDF = testingDF.drop(columns = dropColsStr)
    #print(list(trainingDF.columns))
    if regressionType == "randomForrest":
        X_training , feature_labels = featureFiltering(outputDir , trainingDF , list(trainingDF.columns), featureStr)
        newCols = set(X_training.columns)
        testCols = set(testingDF.columns)
        removeCols = list(testCols - newCols)
        testingDF = testingDF.drop(columns = removeCols)
        hyperFileRF = Path(outputDir) / "randomForest" / "hyperParameter" / "params.dat"
        if not os.path.exists(str(Path(outputDir) / "randomForest" / "hyperParameter")): 
            os.makedirs(str(Path(outputDir) / "randomForest" / "hyperParameter") ) 
        reducedDimDF = rFR_2Source(X_training, yTraining , testingDF, trainingYStr,  hyperFileRF, outputDir)
        dfFinal = addToDF(reducedDimDF, keepHash , rejectedDF)
        return dfFinal
    elif regressionType == "MLP":
        X_training , feature_labels = featureFiltering(outputDir , trainingDF , list(trainingDF.columns), featureStr)
        newCols = set(X_training.columns)
        testCols = set(testingDF.columns)
        removeCols = list(testCols - newCols)
        testingDF = testingDF.drop(columns = removeCols)
        hyperFileMLP = Path(outputDir) / "MLP" / "hyperParameter" / "params.dat"
        if not os.path.exists(str(Path(outputDir) / "MLP" / "hyperParameter")): 
            os.makedirs(str(Path(outputDir) / "MLP" / "hyperParameter") ) 
        MLPipe = Pipeline([
            ("scaler" , StandardScaler()),
            ("mlp" , MLPRegressor(
                solver = "adam" , 
                max_iter = 2000,
                early_stopping = True,
                validation_fraction = 0.1,
                n_iter_no_change = 20,
                random_state = 0
            ))
        ])
        mlpGrid= { "mlp__hidden_layer_sizes" : [(64,32) , (128,64) , (128,64,32)], 
                  "mlp__activation" : ["relu" , "tanh"], 
                  "mlp__alpha" : [1e-4 , 1e-3 ]}
        reducedDimDF = regressionPipeline(X_training, yTraining, testingDF, trainingYStr, MLPRegressor, param_grid=mlpGrid,
                       hyperParmFile=hyperFileMLP, outputDir=outputDir, PipeObj=MLPipe, model_name="MLP")
        dfFinal = addToDF(reducedDimDF, keepHash , rejectedDF)
        return dfFinal
    elif regressionType == "XGBOOST":
        X_training , feature_labels = featureFiltering(outputDir , trainingDF , list(trainingDF.columns), featureStr)
        newCols = set(X_training.columns)
        testCols = set(testingDF.columns)
        removeCols = list(testCols - newCols)
        testingDF = testingDF.drop(columns = removeCols)
        hyperFileMLP = Path(outputDir) / "GradientBoost" / "hyperParameter" / "params.dat"
        if not os.path.exists(str(Path(outputDir) / "GradientBoost" / "hyperParameter")): 
            os.makedirs(str(Path(outputDir) / "GradientBoost" / "hyperParameter") ) 
        gbGrid = {
            'n_estimators': [100, 200, 300],
            'max_depth': [2, 4, 8, None],
            'learning_rate': [0.01, 0.05, 0.1, 0.2],
            'max_features': ['sqrt', 0.5, None]
        }
        GBPipe = Pipeline([
            ("scaler" , StandardScaler()),
            ("gbR" , GradientBoostingRegressor(
                n_estimators=100,
                random_state = 0
            ))
        ])
        reducedDimDF = regressionPipeline(X_training, yTraining, testingDF, trainingYStr, MLPRegressor, param_grid=gbGrid,
                       hyperParmFile=hyperFileMLP, outputDir=outputDir, PipeObj=GBPipe, model_name="MLP")
        dfFinal = addToDF(reducedDimDF, keepHash , rejectedDF)
        return dfFinal
if __name__ == "__main__":
    trainingCSV = str(sys.argv[1])
    testingCSV = str(sys.argv[2])
    outputDir = str(sys.argv[3])
    if not os.path.exists(outputDir ): 
        os.makedirs(outputDir) 
    regressionStr = str(sys.argv[4])
    featureType = str(sys.argv[5])
    axisDF = main(trainingCSV , testingCSV , outputDir  , regressionStr , featureType) 
    finalDir = Path(outputDir) / "finalDF"
    finalDir.mkdir(parents=True, exist_ok=True)
    finalFile = finalDir / f"{featureType}_ReducedDim.csv"

    axisDF.to_csv(finalFile, index=False)