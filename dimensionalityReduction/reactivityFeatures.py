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
def featureOverYield(outputDir , reactivityStr, chemistries, dataDir):
    if "csv" in dataDir:
        dfMAST = pd.read_csv(dataDir)
    elif "xlsx" in dataDir:
        dfMAST = pd.read_excel(dataDir)
    chemDFMAST = pd.DataFrame()
    for chemistry in chemistries:
        df = pd.read_csv(chemistry)
        stringCols = list(df.select_dtypes(include='string').columns)
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
                
                elif convert == 2:
                    
                break
            except ValueError:
                print(f"Invalid input. Please enter either [1] or [2]\n")  
if __name__ == "__main__":
    chemDir = str(sys.argv[1])
    dataDir = str(sys.argv[2])
    chemList  , Yield= selectChemistries(chemDir)
    outputDir = input("Enter an output directory for the scatterplots: ")
    featureOverYield(outputDir , Yield, chemList , dataDir )

