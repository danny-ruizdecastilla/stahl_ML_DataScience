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
from sklearn.model_selection import train_test_split
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.metrics import mean_squared_error, r2_score
import shap
import sklearn
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from DFTWorkflow.pitchingATent import featureFiltering
from dimensionalityReduction.reactivityFeatures import boxGen
from ml_Regression.multiExperimentRegression import listInputs , encodeReactionInfo
from sklearn.dummy import DummyRegressor 
def plotly_template():  # Credit to Dylan Walsh
    template = go.layout.Template()
    template.layout.font = dict(family="Arial", size=18, color="black")
    template.layout.plot_bgcolor = "white"
    template.layout.width, template.layout.height = 1200, 600
    template.layout.xaxis.tickprefix = "<b>"
    template.layout.xaxis.ticksuffix = "<b>"
    template.layout.xaxis.showline = True
    template.layout.xaxis.linewidth = 5
    template.layout.xaxis.linecolor = "black"
    template.layout.xaxis.ticks = "outside"
    template.layout.xaxis.tickwidth = 4
    template.layout.xaxis.showgrid = False
    template.layout.xaxis.mirror = True
    template.layout.xaxis.zeroline = False
    template.layout.xaxis.tickfont = dict(size=12)
    template.layout.yaxis.tickprefix = "<b>"
    template.layout.yaxis.ticksuffix = "<b>"
    template.layout.yaxis.showline = True
    template.layout.yaxis.linewidth = 5
    template.layout.yaxis.linecolor = "black"
    template.layout.yaxis.ticks = "outside"
    template.layout.yaxis.tickwidth = 4
    template.layout.yaxis.showgrid = False
    template.layout.yaxis.mirror = True
    template.layout.yaxis.zeroline = False
    template.layout.yaxis.tickfont = dict(size=12)

    return template
def greedyBackwardElimination( X , yCol , saveStr , nSeeds , dummyFinal  , hyperFile , targetCols):
    j = 0
    rmse = []
    r2 = []
    worstFeats = {}
    targetCols = int(len(list(X.columns))*0.5)
    while j < 10:
        rmseTest , r2Test , shapDF = gbRegressorQuick(X , yCol , j, hyperFile )
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
    print(rmseToBeat , r2ToBeat)
    sortedWorst = sorted(
        ((k, v) for k, v in worstFeats.items() if v >= 3),
        key=lambda kv: kv[1], reverse=True
    )

    candidates = {}

    for feature, _ in sortedWorst:
        print(feature)
        newX = X.drop(columns=[feature])

        rmseF = []
        r2F = []
        for j in range(nSeeds):
            rmseTest, r2Test, shapDF = gbRegressorQuick(newX, yCol, j, hyperFile)
            rmseF.append(rmseTest)
            r2F.append(r2Test)

        rmseComp = np.mean(rmseF)
        r2Comp = np.mean(r2F)
        print(rmseComp , r2Comp)
        if rmseComp < rmseToBeat and r2Comp > r2ToBeat:
            candidates[feature] = [
                np.abs(rmseToBeat - rmseComp),
                np.abs(r2ToBeat - r2Comp),
                rmseComp,
                r2Comp
            ]
    savePath = hyperFile.parent
    results = savePath  / "dropFeatures.dat"
    if not candidates:
        RMSEMinimizationGraph(results , dummyFinal, nSeeds , saveStr)
        return X, rmseToBeat, r2ToBeat
    else:
        bestFeature = max(candidates.items(), key=lambda kv: kv[1][0])[0]
        bestRmse = candidates[bestFeature][2]
        bestR2 = candidates[bestFeature][3]

        newReducedX = X.drop(columns=[bestFeature])
        newCols = len(newReducedX.columns)
        with open(results , "w") as f:
            f.write(f"Dropped: {bestFeature} | new RMSE: {bestRmse:.5f} (was {rmseToBeat:.5f}) | new R2: {bestR2:.5f} (was {r2ToBeat:.5f}) | totFeats: {newCols}\n")

        return greedyBackwardElimination( newReducedX , yCol , saveStr , nSeeds , dummyFinal  , hyperFile , targetCols)
def gbRegressorQuick(X , y , j , hyperParamFile):
    X_train_CV, X_test, y_train_CV, y_test = train_test_split(X, y, test_size=0.25, random_state=j)

    with open(hyperParamFile, "r") as f:
        savedParms = json.load(f)
    bestParameters = savedParms["best_params"]
    gbFinal = GradientBoostingRegressor(**bestParameters ,random_state=j)
    gbFinal.fit(X_train_CV, y_train_CV)
    yTrain = gbFinal.predict(X_train_CV)
    yPred = gbFinal.predict(X_test)
    rmseTest = float(np.sqrt(mean_squared_error(y_test, yPred)))
    r2Test = float(r2_score(y_test, yPred))
    #print(X_train_CV.shape)
    #print(X_test.shape)
    #print(X_train_CV.columns.equals(X_test.columns))
    gbExplainer = shap.TreeExplainer(gbFinal)
    shapValues = gbExplainer(X_test )
    meanAbsShap = np.abs(shapValues.values).mean(axis=0)
    shapDF = pd.DataFrame({"feature" : X_test.columns , "Shap" : meanAbsShap}).sort_values("Shap" , ascending = False).reset_index(drop=False)
    return rmseTest , r2Test , shapDF
def RMSEMinimizationGraph(resultsPath , dummyRMSE , numSeeds , saveStr):
    resultsHash = {"RMSE" : [] , "numFeats" : [] , "featDropped" : []}
    with open(resultsPath , "r") as f:
        for idx, line in enumerate(f):
            if len(resultsHash["RMSE"]) == 0:
                rmse = float(line.split("|")[1].split("(was")[-1].split(")")[0].strip())
                numCols = int(line.split("| totFeats:")[-1].strip()) + 1
                resultsHash["featDropped"].append("No Features Dropped")
                resultsHash["RMSE"].append(rmse)
                resultsHash["numFeats"].append(numCols)               
            rmse = float(line.split("|")[1].split("RMSE:")[-1].split("(")[0].strip())
            numCols = int(line.split("| totFeats:")[-1].strip())
            featDropped = str(line.split("|")[0].split("Dropped: ")[-1].strip())
            resultsHash["featDropped"].append(featDropped)
            resultsHash["RMSE"].append(rmse)
            resultsHash["numFeats"].append(numCols)
    
    fig = go.Figure(layout=dict(template=plotly_template()))
    names = resultsHash["featDropped"]
    fig.add_trace(go.Scatter(
        x=resultsHash["numFeats"], 
        y=resultsHash["RMSE"],
        mode="lines+markers",
        text=names,                               
        hovertemplate="<b>%{text}</b><br>" +            
                    "Features: %{x}<br>" + 
                    "RMSE: %{y:.4f}<extra></extra>"
    ))
    fig.add_hline(y=dummyRMSE, line_width=2, line_dash="solid", line_color="#990000")
    fig.update_layout(
        xaxis=dict(title='Number of Features Remaining', scaleanchor="y"),  # Keeps x and y scales equal
        yaxis=dict(title=f'Test RMSE of {numSeeds} random splits'),
        plot_bgcolor='rgba(255,255,255,0.1)',  # Light background transparency
        width=600,  
        height=600,  
        margin=dict(l=60, r=60, t=50, b=60),  
        legend=dict(
            x=0.0,  
            y=0.0,  
            bgcolor="rgba(255,255,255,0.7)",  
            bordercolor="black",
            borderwidth=1
        ),
        title=dict(
            text=f"Bootstrap Feature Selection (Greedy Approach) for {saveStr} Feature Class",  # Fixed f-string
            font=dict(size=18, color="black")
        )
    )
    outputDir = resultsPath.parent
    fig.write_html(outputDir / "rmseMinimization.html")
def main(dfMAST , saveDir , saveName , hyperFile , nSeeds):
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
    targetCols = int(len(list(reducedX.columns))*0.5)

    j = 0
    dummyLst = []
    while j < nSeeds:
        X_train_CV, X_test, y_train_CV, y_test = train_test_split(X, yCol, test_size=0.25, random_state=j)
        dummy = DummyRegressor(strategy="median")
        dummy.fit(X_train_CV , y_train_CV)
        yDummy = dummy.predict(X_test)
        dummyRMSE = float(np.sqrt(mean_squared_error(y_test, yDummy)))
        dummyLst.append(dummyRMSE)
        
        j +=1
    dummyRMSEFinal = np.mean(dummyLst)

    finalDF , finalRMSE , finalR2 =  greedyBackwardElimination( reducedX , yCol , saveStr , nSeeds , dummyRMSEFinal , hyperFile , targetCols)
    finalDF[yStr] = yCol 
    savePath = hyperFile.parent
    finalDF.to_csv(savePath / f"{saveStr}_reducedFeatures.csv")
    


if __name__ == "__main__":
    dfStr = str(sys.argv[1])
    dfMain = pd.read_csv(dfStr)
    saveDir = Path(sys.argv[2])
    saveStr = str(sys.argv[3])
    nSeeds = int(sys.argv[4])
    saveDir.parent.mkdir(parents=True, exist_ok=True)
    main(dfMain , saveDir , saveStr , saveDir / "hyperParameterOptimization.dat" , nSeeds)