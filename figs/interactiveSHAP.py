#Danny Ruiz de Castilla 07.30.2026
#Prints out an interactive SHAP plot that can use hover to display chemical structures on a seperate screen 
import pandas as pd
import numpy as np
import re
from pathlib import Path
import shap
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import LinearRegression
import json
import shutil
import html 
import sys 
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from figs.chemPlotlyV2 import createPNGDF,png64
def generateInteractiveSHAP(shapDF, dfMAST ,sortedCols , featureDir):
    shapJSON = shapDF.to_json(orient='records' , force_ascii = False)
    shapJSON = shapJSON.replace('\\/' , '/')

    featJSON = dfMAST.to_json(orient='records' , force_ascii = False)
    featJSON = featJSON.replace('\\/' , '/')

    html = f"""
    <!DOCTYPE html>
    <html>
    <meta charset="UTF-8">
    <head>
        <title>Interactive SHAP for Features </title>
        <script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script>
        <script src="https://unpkg.com/simple-statistics@7.8.3/dist/simple-statistics.min.js"></script>
        <script>
            const jsonSHAP = {json.dumps(shapJSON)};
            const jsonFeats = {json.dumps(featJSON)};
            const figTrace = [];

            funciton plotData(){{
                for (const cols in {sortedCols.values()}){{
                    const x = jsonSHAP.map(p => p[cols]);
                    const featVals = jsonFeats.map(p => p[cols]);
                    cons y = 
                }}
            }}
        </script>
        <style>
            body{{
                font-family: Arial, sans-serif;
                margin: 20px;
                background-color: #f5f5f5;
            }}
            .plot-container {{
                display: flex;
                gap: 20px;
            }}
            #plotContainer {{
                background-color: white;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                flex: 1;
                min-height: 500px;
            }}
            #imageContainer {{
                width: 300px;
                background-color: white;
                padding: 20px;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            #hoverImage {{
                max-width: 100%;
                height: auto;
                border: 1px solid #ddd;
                border-radius: 4px;
            }}
            .image-placeholder {{
                color: #666;
                font-style: italic;
                text-align: center;
                padding: 50px 0;
            }}

    """

def main(df , yCol  , idStr , SMILESStr , outputDir):
    while True:
        modelInt = input(f"Please Select the index for the model you want to fit to:\n\n[1]  XGBOOST\n\n[2]   RandomForest\n\n[3]    Linear Regression\n\n").strip()
        if modelInt == "1":
            model = GradientBoostingRegressor(n_estimators=300 ,max_depth=4 , learning_rate=0.05 )
            break
        elif modelInt == "2" : 
            model = RandomForestRegressor(n_estimators=300,max_depth=4,random_state=42,n_jobs=-1)
            break
        elif modelInt == "3":
            model = LinearRegression()
        else:
            print("Invalid input, enter 1 or 2 only")

    yVals = df[yCol]
    idVals = df[idStr]    
    figDir = Path(outputDir) / "png" 
    dfMAST = createPNGDF(dfMAST ,"SMILES" , str(figDir))
    base64Col = []
    for img in list(dfMAST["pngPath"]):

        base64 = png64(img)
        base64Col.append(base64)
    dfMAST = dfMAST.drop(columns = [yCol , idStr , SMILESStr , "pngPath"])
    model.fit(dfMAST , yVals)
    if modelInt == 3:
        modelExplainer = shap.LinearExplainer(model , dfMAST)
    else:
        modelExplainer = shap.TreeExplainer(model , dfMAST)
    shapValues = modelExplainer.shap_values(dfMAST)
    shapDF = pd.DataFrame(shapValues,index=dfMAST.index,columns=dfMAST.columns)
    shutil.rmtree(figDir)

    featureColsHash = {}
    for idx , feature in enumerate(shapDF.columns):
        shapVal = np.mean([np.abs(val) for val in list(shapDF[feature])])
        featureColsHash[feature] = shapVal
    sortedCols = sorted(
        featureColsHash.items(),
        key=lambda kv: kv[1]
    )
    shapDF["ID"] = idVals
    shapDF["Base64"] = base64Col
    featureDir = Path(outputDir) / "features" 
    featureDir.parent.mkdir(parents=True, exist_ok=True)

    generateInteractiveSHAP(shapDF, dfMAST ,sortedCols , featureDir)




