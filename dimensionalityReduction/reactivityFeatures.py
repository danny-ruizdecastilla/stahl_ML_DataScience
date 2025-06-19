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
def featureOverYield(outputDir , reactivityStr, chemistires, dataDir):
     
if __name__ == "__main__":
    dataDir = str(sys.argv[1])
    chemList  , reactivityParm= selectChemistries(dataDir)
    outputDir = input("Enter an output directory for the scatterplots: ")


