#!/usr/bin/env python3
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
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.svm import SVR
from sklearn.metrics import mean_squared_error, r2_score
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
def saveScatter(yPredTest, yTest, yPredTrain , yTrain , outputDir):
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=yPredTrain,
            y=yTrain,
            mode="markers",
            marker=dict(
                size=8,
                color="rgba(001,031,091,0.8)",
                line=dict(color="black", width=1)
            ),
            name="Train"
        )
    )

    fig.add_trace(
        go.Scatter(
            x=yPredTest,
            y=yTest,
            mode="markers",
            marker=dict(
                size=10,
                color="rgba(153,000,000,0.8)",
                line=dict(color="black", width=1)
            ),
            name="Test"
        )
    )

    minVal = min(
        np.min(yPredTrain),
        np.min(yPredTest),
        np.min(yTrain),
        np.min(yTest)
    )

    maxVal = max(
        np.max(yPredTrain),
        np.max(yPredTest),
        np.max(yTrain),
        np.max(yTest)
    )

    fig.add_trace(
        go.Scatter(
            x=[minVal, maxVal],
            y=[minVal, maxVal],
            mode="lines",
            line=dict(color="black", dash="dash"),
            name="Ideal"
        )
    )

    fig.write_html(Path(outputDir) / "gradientBoost" / "Eval" / "parityPlot.html")
class CustomError(Exception):
    pass

def stratifiedRegressionSplit(y, j , n_splits=5, n_bins=10 ):
    y_ser = pd.Series(y).reset_index(drop=True)
    try:
        yBinned = pd.qcut(y_ser, q=n_bins, labels=False, duplicates='drop')
    except Exception:
        # fallback to equal width
        yBinned = pd.cut(y_ser, bins=n_bins, labels=False)

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=j)
    return skf.split(range(len(y_ser)), yBinned)

def stratifiedRegressionCV(X, y, j , model, standardScale ,n_splits=5 ):
    if standardScale:
        scaled = StandardScaler()
        X = scaled.fit_transform(X)
        X_scaled = scaled.fit_transform(X)
        X = pd.DataFrame(X_scaled, columns=X.columns, index=X.index)
    is_df = isinstance(X, pd.DataFrame)
    X_arr = X.values if is_df else np.asarray(X)
    y_arr = np.asarray(y)

    fold_generator = stratifiedRegressionSplit(y_arr, j , n_splits=n_splits )

    mse_scores = []
    r2_scores = []

    for fold_idx, (train_idx, val_idx) in enumerate(fold_generator):
        X_train_fold, X_val_fold = X_arr[train_idx], X_arr[val_idx]
        y_train_fold, y_val_fold = y_arr[train_idx], y_arr[val_idx]

        model.fit(X_train_fold, y_train_fold)
        y_pred = model.predict(X_val_fold)

        mse = mean_squared_error(y_val_fold, y_pred)
        r2 = r2_score(y_val_fold, y_pred)

        mse_scores.append(mse)
        r2_scores.append(r2)

    return mse_scores, r2_scores

def randomForestRegression(X, y, hyperParmFile, outputDir , j ):
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

    mseCV, r2CV = stratifiedRegressionCV(X_train_CV, y_train_CV,  j , rfMAST, True ,  n_splits=5)

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

def supportVectorRegression(kernelStr, X, y, hyperParmFile, outputDir):
    X_train_CV, X_test, y_train_CV, y_test = train_test_split(X, y, test_size=0.2, random_state=j)

    if kernelStr == "linear":
        svrGrid = {'kernel': ['linear'], 'C': [0.1, 1, 10, 100, 1000], 'epsilon': [0.01, 0.1, 0.2, 0.5]}
    elif kernelStr == "rbf":
        svrGrid = {'kernel': ['rbf'], 'C': [0.1, 1, 10, 100, 1000], 'gamma': ['scale', 0.001, 0.01, 0.1, 1], 'epsilon': [0.01, 0.1, 0.2, 0.5]}
    elif kernelStr == "poly":
        svrGrid = {'kernel': ['poly'], 'C': [0.1, 1, 10, 100], 'gamma': ['scale', 0.01, 0.1, 1], 'degree': [2, 3, 4], 'epsilon': [0.01, 0.1, 0.2]}
    else:
        raise ValueError("Unsupported kernel string")

    if hyperParmFile.exists():
        with open(hyperParmFile, "r") as f:
            savedParms = json.load(f)
        best_params = savedParms["best_params"]
        # SVR doesn't accept random_state, so just pass params directly
        svrMAST = SVR(**best_params)
        svrFinal = SVR(**best_params)
    else:
        svr = SVR()
        gridSearch = GridSearchCV(svr, svrGrid, cv=5, scoring="neg_mean_squared_error", n_jobs=-1)
        gridSearch.fit(X_train_CV, y_train_CV)
        with open(hyperParmFile, "w") as f:
            json.dump({"best_score": gridSearch.best_score_, "best_params": gridSearch.best_params_}, f, indent=4)
        svrMAST = SVR(**gridSearch.best_params_)
        svrFinal = SVR(**gridSearch.best_params_)

    mseCV, r2CV = stratifiedRegressionCV(X_train_CV, y_train_CV, svrMAST, True , n_splits=5)

    cvFile = Path(outputDir) / "supportVector" / "crossvalidation" / "scores.dat"
    cvFile.parent.mkdir(parents=True, exist_ok=True)
    with open(cvFile, "w") as f:
        for i in range(len(mseCV)):
            f.write(f"Fold {i} Results: mse {mseCV[i]} | r2 {r2CV[i]}\n")

    svrFinal.fit(X_train_CV, y_train_CV)
    yTrain = svrFinal.predict(X_train_CV)
    yPred = svrFinal.predict(X_test)

    rmseTest = float(np.sqrt(mean_squared_error(y_test, yPred)))
    r2Test = float(r2_score(y_test, yPred))

    rmseTrain = float(np.sqrt(mean_squared_error(y_train_CV , yTrain)))
    r2Train = float(r2_score(y_train_CV , yTrain))

    testFile = Path(outputDir) / "supportVector" / "Eval" / "testScores.dat"
    testFile.parent.mkdir(parents=True, exist_ok=True)
    with open(testFile, "w") as f:
        f.write(f"Test Results: rmse {rmseTest} | r2 {r2Test}\n")

    trainFile = Path(outputDir) / "supportVector" / "Eval" / "trainScores.dat"
    with open(trainFile, "w") as f:
        f.write(f"Train Results: rmse {rmseTrain} | r2 {r2Train}\n")

def gradientBoostRegression(X, y, hyperParmFile, outputDir , j , saveFig : False):
    X_train_CV, X_test, y_train_CV, y_test = train_test_split(X, y, test_size=0.2, random_state=j)

    if hyperParmFile.exists():
        with open(hyperParmFile, "r") as f:
            savedParms = json.load(f)
        best_params = savedParms["best_params"]
        gbMAST = GradientBoostingRegressor(**best_params, random_state=j)
        gbFinal = GradientBoostingRegressor(**best_params, random_state=j)
    else:
        gbGrid = {
            'n_estimators': [100, 200, 300],
            'max_depth': [2, 4, 8, None],
            'learning_rate': [0.01, 0.05, 0.1, 0.2],
            'max_features': ['sqrt', 0.5, None]
        }
        gb = GradientBoostingRegressor(random_state=j)
        gridSearch = GridSearchCV(gb, gbGrid, cv=5, scoring="neg_mean_squared_error", n_jobs=-1)
        gridSearch.fit(X_train_CV, y_train_CV)
        with open(hyperParmFile, "w") as f:
            json.dump({"best_score": gridSearch.best_score_, "best_params": gridSearch.best_params_}, f, indent=4)
        gbMAST = GradientBoostingRegressor(**gridSearch.best_params_, random_state=j)
        gbFinal = GradientBoostingRegressor(**gridSearch.best_params_, random_state=j)

    mseCV, r2CV = stratifiedRegressionCV(X_train_CV, y_train_CV, j , gbMAST, False , n_splits=5) #X, y, j , model, standardScale ,n_splits=5

    cvFile = Path(outputDir) / "gradientBoost" / "crossvalidation" / "scores.dat"
    cvFile.parent.mkdir(parents=True, exist_ok=True)
    with open(cvFile, "w") as f:
        for i in range(len(mseCV)):
            f.write(f"Fold {i} Results: mse {mseCV[i]} | r2 {r2CV[i]}\n")

    gbFinal.fit(X_train_CV, y_train_CV)
    yTrain = gbFinal.predict(X_train_CV)
    yPred = gbFinal.predict(X_test)

    rmseTest = float(np.sqrt(mean_squared_error(y_test, yPred)))
    r2Test = float(r2_score(y_test, yPred))

    rmseTrain = float(np.sqrt(mean_squared_error(y_train_CV , yTrain)))
    r2Train = float(r2_score(y_train_CV , yTrain))

    testFile = Path(outputDir) / "gradientBoost" / "Eval" / "testScores.dat"
    testFile.parent.mkdir(parents=True, exist_ok=True)
    with open(testFile, "w") as f:
        f.write(f"Test Results: rmse {rmseTest} | r2 {r2Test}\n")

    trainFile = Path(outputDir) / "gradientBoost" / "Eval" / "trainScores.dat"
    with open(trainFile, "w") as f:
        f.write(f"Train Results: rmse {rmseTrain} | r2 {r2Train}\n")
    if saveFig:
        saveScatter(yPred, y_test, yTrain, y_train_CV, outputDir )

    importances = gbFinal.feature_importances_
    if isinstance(X_train_CV, pd.DataFrame):
        feature_names = X_train_CV.columns
    else:
        feature_names = [f"f{i}" for i in range(importances.shape[0])]
    featDF = pd.DataFrame({"feature": feature_names, "importance": importances}).sort_values(by="importance", ascending=False)
    featFile = Path(outputDir) / "gradientBoost" / "Eval" / "testWeights.csv"
    featDF.to_csv(featFile, index=False)

def chooseRegressionModels():
    models = {
        0: "Support Vector Machine (SVR)",
        1: "Random Forest Regressor",
        2: "Gradient Boosting Regressor",
    }

    print("Select from the available regression models (separate multiple choices with commas):")
    for i, name in models.items():
        print(f"[{i}] {name}")

    choice_str = input("Enter your choice(s): ")
    try:
        choices = [int(c.strip()) for c in choice_str.split(",")]
    except ValueError:
        print("Invalid input. Please enter integers separated by commas.")
        return []
    valid_choices = [c for c in choices if c in models]
    if not valid_choices:
        print("No valid models selected.")
        return []
    return valid_choices, models

if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise CustomError("Usage: script.py <datasetDir> <outputDir>")

    datasetDir = str(sys.argv[1])
    datasets = glob.glob((datasetDir + "/*.csv"))
    mainOutputDir = str(sys.argv[2])

    if len(datasets) == 1:
        dfMAST = pd.read_csv(datasets[0])
        regressionStr = input(f"Please enter the name of the column in {datasets[0]} that corresponds to the regression metric: ").strip()
        if regressionStr not in dfMAST.columns:
            raise CustomError(f"Column {regressionStr} not found in {datasets[0]}")

        y = dfMAST[regressionStr]
        smilesStr = input(f"Please enter the name of the column in {datasets[0]} that corresponds to SMILES strings (or press Enter if none): ").strip()

        Xdataframe = dfMAST.drop(columns=[smilesStr, regressionStr])
        hyperFileRF = Path(mainOutputDir) / "randomForest" / "hyperParameter" / "params.dat"
        hyperFileRF.parent.mkdir(parents=True, exist_ok=True)
        hyperFileSVM = Path(mainOutputDir) / "supportVector" / "hyperParameter" / "params.dat"
        hyperFileSVM.parent.mkdir(parents=True, exist_ok=True)
        hyperFileGB = Path(mainOutputDir) / "gradientBoost" / "hyperParameter" / "params.dat"
        hyperFileGB.parent.mkdir(parents=True, exist_ok=True)
        print(Xdataframe.shape)
        j = 0
        while j < 4:

            newOutputDir = Path(mainOutputDir) / f"trial{j}"
            newOutputDir.mkdir(parents=True, exist_ok=True)
            gradientBoostRegression(Xdataframe, y, hyperFileGB, newOutputDir)
            randomForestRegression(Xdataframe, y, hyperFileRF, newOutputDir , j )
            supportVectorRegression("rbf", Xdataframe, y, hyperFileSVM, newOutputDir)
            j +=1

    elif len(datasets) == 0:
        raise CustomError(f"{datasetDir} does not have any CSV files")
    else:
        regressionStr = input(f"Please enter the name of the column in the datasets that corresponds to the regression metric: ").strip()
        colDrops = ["SMILES", "Unnamed"]
        Xdataframe, smileList, yieldList = compressData(sorted(datasets), regressionStr, colDrops)
