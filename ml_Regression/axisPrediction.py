#Predicts the rank positioning on a 1-d axis based on a feature matrix
import torch
import torch.nn as nn
import os 
import sys
import glob
import pandas as pd 
import numpy as np
import chemdraw
import re
from pathlib import Path
from rdkit import Chem
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, r2_score
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from dimensionalityReduction.reactivityFeatures import boxGen
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from ml_Regression.ML_RegressionV2 import GridSearchCV , stratifiedRegressionCV
def canonical_smiles(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None
def rFR_2Source(X, y, hyperParmFile, outputDir):
    import json 

    # keep DataFrame so we can get column names for feature importances
    X_train_CV, X_test, y_train_CV, y_test = train_test_split(X, y, test_size=0.2, random_state=j)

    if hyperParmFile.exists():
        with open(hyperParmFile, "r") as f:
            savedParms = json.load(f)
        best_params = savedParms["best_params"]
        rfMAST = RandomForestRegressor(**best_params, random_state=j)
        rfFinal = RandomForestRegressor(**best_params, random_state=j)
    else:
        rfGrid = {
            'n_estimators': [100, 200, 300],
            'max_depth': [2, 4, 8, None],
            'max_features': ['sqrt', 0.5, None]
        }
        rf = RandomForestRegressor(random_state=j)
        gridSearch = GridSearchCV(rf, rfGrid, cv=5, scoring="neg_mean_squared_error", n_jobs=-1)
        gridSearch.fit(X_train_CV, y_train_CV)
        with open(hyperParmFile, "w") as f:
            json.dump({"best_score": gridSearch.best_score_, "best_params": gridSearch.best_params_}, f, indent=4)
        rfMAST = RandomForestRegressor(**gridSearch.best_params_, random_state=j)
        rfFinal = RandomForestRegressor(**gridSearch.best_params_, random_state=j)

    mseCV, r2CV = stratifiedRegressionCV(X_train_CV, y_train_CV, rfMAST, n_splits=5)

    cvFile = Path(outputDir) / "randomForest" / "crossvalidation" / "scores.dat"
    cvFile.parent.mkdir(parents=True, exist_ok=True)
    with open(cvFile, "w") as f:
        for i in range(len(mseCV)):
            f.write(f"Fold {i} Results: mse {mseCV[i]} | r2 {r2CV[i]}\n")

    rfFinal.fit(X_train_CV, y_train_CV)
    yTrain = rfFinal.predict(X_train_CV)
    yPred = rfFinal.predict(X_test)

    rmseTest = float(np.sqrt(mean_squared_error(y_test, yPred)))
    r2Test = float(r2_score(y_test, yPred))

    rmseTrain = float(np.sqrt(mean_squared_error(y_train_CV , yTrain)))
    r2Train = float(r2_score(y_train_CV , yTrain))

    testFile = Path(outputDir) / "randomForest" / "Eval" / "testScores.dat"
    testFile.parent.mkdir(parents=True, exist_ok=True)
    with open(testFile, "w") as f:
        f.write(f"Test Results: rmse {rmseTest} | r2 {r2Test}\n")

    trainFile = Path(outputDir) / "randomForest" / "Eval" / "trainScores.dat"
    with open(trainFile, "w") as f:
        f.write(f"Train Results: rmse {rmseTrain} | r2 {r2Train}\n")

    importances = rfFinal.feature_importances_
    if isinstance(X_train_CV, pd.DataFrame):
        feature_names = X_train_CV.columns
    else:
        feature_names = [f"f{i}" for i in range(importances.shape[0])]
    featDF = pd.DataFrame({"feature": feature_names, "importance": importances}).sort_values(by="importance", ascending=False)
    featFile = Path(outputDir) / "randomForest" / "Eval" / "testWeights.csv"
    featDF.to_csv(featFile, index=False)
def main(training , testing , outputDir  , regressionType ):
    trainingDF = pd.read_csv(training, encoding="utf-8" )
    testingDF = pd.read_csv(testing, encoding="utf-8")
    trainingCols = list(trainingDF.columns)
    testingCols = list(testingDF.columns)
    trainingBox = boxGen(trainingCols)
    testingBox = boxGen(testingCols)

    smilesIdxTraining = int(input(f"Enter the index of the smiles column:\n {trainingBox}\n"))
    smilesIdxTesting =  int(input(f"Enter the index of the smiles column:\n {testingBox}\n"))

    trainingDF["canonical_smiles"] = [canonical_smiles(smi) for smi in trainingDF[trainingCols[smilesIdxTraining]]]
    trainingCanonicalSet = set(trainingDF["canonical_smiles"])

    testingDF["canonical_smiles"] = [canonical_smiles(smi) for smi in testingDF[testingCols[smilesIdxTesting]]]
    testingDF = testingDF[~testingDF["canonical_smiles"].isin(trainingCanonicalSet)].reset_index(drop=True)

    testingDF = testingDF.drop(columns=["canonical_smiles"])
    trainingDF = trainingDF.drop(columns=["canonical_smiles"])


    if regressionType == "randomForrest":

if __name__ == "__main__":
    trainingCSV = str(sys.argv[1])
    testingCSV = str(sys.argv[2])
    outputDir = str(sys.argv[3])
    if not os.path.exists(outputDir ): 
        os.makedirs(outputDir) 
    regressionStr = str(sys.argv[4])
    main(trainingCSV , testingCSV , outputDir  , regressionStr )