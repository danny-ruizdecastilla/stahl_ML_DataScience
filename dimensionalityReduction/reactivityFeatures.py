import os 
import sys
import glob
import shutil
import re
import pandas as pd 
import numpy as np
import plotly 
import plotly.graph_objects as go
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from DFTWorkflow.pitchingATent import convertCanonical
from figs.chemPlotlyV1 import plotly_template
def boxGen(list_):
    lines = []
    for i , line_ in enumerate(list_):
        str = f"{line_} == [{i}]"
        lines.append(str)
    maxLine = max(len(line) for line in lines)
    Box = "+" + "-"*(maxLine) + "+\n"
    for line in lines:
        Box += "| " + line.ljust(maxLine) + " |\n"
    Box +=  "+" + "-"*(maxLine) + "+\n" 
    return Box
def selectChemistries(mainDir ):
    csvPool = glob.glob(mainDir + "/*.csv")
    csvOptions = [csv.split("/")[-1] for csv in csvPool]
    csvBox = boxGen(csvOptions)
    while True:
        prompt1 = f"""Here are the dataframe options from this directory:\n{csvBox}
                Please enter the corresponding numbers for each directory you want to generate scatterplots for
                """
        indexList= listInputs(prompt1)
        indexList = [int(index) for index in indexList]
        prompt2 = input("Enter the column name inside the dataframes that you want to use as your Reactivity Metric: ")
        csvList = [csvOptions[index] for index in indexList]
        badDF = []
        for csv in csvList:
            csvDir = mainDir + "/" + csv
            df = pd.read_csv(csvDir)
            if not prompt2 in list(df.columns):
                badDF.append(csv)
        if len(badDF) == 0:
            fileList = csvList
            break
        else:
            rejects = boxGen(badDF)
            prompt3 = f"""These are the dataframes without the required Reactivity Metric: {prompt2}\n{rejects}
            Should we proceed without these dataframes? [1] for Yes and [2] for No
            """
            proceed = input(prompt3)
            while True:
                if proceed.strip() == "2":
                    sys.close()
                elif proceed.strip() == "1":
                    fileList = list(set(csvList) - set(badDF))
                    break
                else:
                    proceed = input("Try again: [1] for Yes and [2] for No")

            break 
    return fileList , prompt2
def plotFeatures(df , outputDir, reactivityStr , featCols):
    figDict = {}
    for feat in featCols:
        xVals = df[feat]
        yVals = df[reactivityStr]
        fig = go.Figure(layout=dict(template=plotly_template()))
        fig.add_trace(go.Scatter(
            x=xVals, 
            y=yVals, 
            mode="markers", 
            name= str(feat) ,
            marker=dict(color='blue' , size = 10), 
            opacity=0.4
        ))
        fig.update_layout(
        xaxis=dict(title=str(feat)),  # Keeps x and y scales equal
        yaxis=dict(title=str(reactivityStr)),
        plot_bgcolor='rgba(255,255,255,0.1)',  # Light background transparency
        title=dict(
            text=("Yield for " + str(feat)),  # Fixed concatenation
            font=dict(size=18, color="black"),  
            x=0.5,  # Center the title
            y=0.95  
        )
    )
        figDict[feat] = fig
    os.makedirs(outputDir, exist_ok=True)

    html_content = "<html><head><script src='https://cdn.plot.ly/plotly-3.0.1.min.js'></script></head><body>\n"
    for feat, fig in figDict.items():
        html_content += fig.to_html(full_html=False, include_plotlyjs=False)
        html_content += "<hr>\n"

    html_content += "</body></html>"

    with open(os.path.join(outputDir, "feature_plots.html"), "w") as f:
        f.write(html_content) 

def featureOverYield(outputDir , reactivityStr, chemistries, dataDir):
    try:
        if dataDir.lower().endswith('.csv'):
            dfMAST = pd.read_csv(dataDir , encoding="utf-8")
        elif dataDir.lower().endswith(('.xlsx', '.xls')):
            dfMAST = pd.read_excel(dataDir)
        else:
            raise ValueError("Master dataset must be CSV or Excel format")
    except Exception as Error:
        print(f"Error loading master dataset: {Error}")
        return None
    chemDFMAST = pd.DataFrame()
    for chemistry in chemistries:
        df = pd.read_csv(chemDir + "/" + chemistry , encoding="utf-8")
        stringCols = list(df.columns)
        stringBox = boxGen(stringCols)
        prompt1 = f"""These are the columns of dtype strings:\n {stringBox}
        Using the corresponding number, select the one with the chemical structures of choice:\n
        """
        while True:
            try:
                smileInd = int(input(prompt1))
                break
            except ValueError:
                print(f"Invalid input. Please enter a single whole number:\n")

        smilesCols = stringCols[smileInd]
        convert = input(f"Do you want to convert to canonical SMILES (recommended)?\n[1] == Yes\n[2] == No")
        while True:
            try:
                convert = int(convert)
                if convert == 1:
                    canonicalList = []
                    smiles = df[smilesCols]
                    for smile in smiles:
                        canonical = convertCanonical(smile)
                        canonicalList.append(canonical)
                    break  
                elif convert == 2:
                    canonicalList = df[smilesCols]
                    break
            except ValueError:
                print(f"Invalid input. Please enter either [1] or [2]\n") 
        tempDF = pd.DataFrame({"canonicalSMILES": canonicalList , "Yield" : list(df[reactivityStr])})
        chemDFMAST = pd.concat([chemDFMAST, tempDF], ignore_index=True)
    stringCols2 = list(dfMAST.select_dtypes(include='object').columns)
    stringBox2 = boxGen(stringCols2)
    prompt = f"""These are the columns of dtype strings:\n {stringBox2}
    Using the corresponding number, select the one with the chemical structures of choice:\n
    """
    while True:
        try:
            smileInd = int(input(prompt))
            break
        except ValueError:
            print(f"Invalid input. Please enter a single whole number:\n") 
    smilesStr = stringCols2[smileInd]
    convert = input(f"Do you want to convert to canonical SMILES (recommended)?\n[1] == Yes\n[2] == No\n")
    while True:
        try:
            convert = int(convert)
            if convert == 1:
                canonicalList = []
                smiles = dfMAST[smilesStr]
                for smile in smiles:
                    canonical = convertCanonical(smile)
                    canonicalList.append(canonical)
                dfMAST["canonicalSMILES"] = canonicalList
                smileStrMAST = "canonicalSMILES"
                break  
            elif convert == 2:
                smileStrMAST = smilesStr
                break
        except ValueError:
            print(f"Invalid input. Please enter either [1] or [2]\n") 
    dropIdx = []
    newCols = [col for col in dfMAST.columns if col not in stringCols2]
    dfMAST = dfMAST.drop_duplicates(subset=[smileStrMAST])
    featureDict = dfMAST.set_index(smileStrMAST).to_dict('index')
    for index, row in chemDFMAST.iterrows():
        canonical = row["canonicalSMILES"]
        if canonical in featureDict:
            matchingFeats = featureDict[canonical]
            for col in newCols:
                chemDFMAST.at[index, col] = matchingFeats[col]
        else:
            dropIdx.append(index)
    chemDFMAST = chemDFMAST.drop(dropIdx)
    plotFeatures(chemDFMAST , outputDir , reactivityStr, newCols)

if __name__ == "__main__":
    chemDir = str(sys.argv[1])
    dataDir = str(sys.argv[2])
    chemList  , Yield= selectChemistries(chemDir)
    outputDir = input("Enter an output directory for the scatterplots: ")
    featureOverYield(outputDir , Yield, chemList , dataDir )

