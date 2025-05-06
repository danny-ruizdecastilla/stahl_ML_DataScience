import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
from rdkit.Chem.PandasTools import LoadSDF
import matplotlib.pyplot as plt
import random
import networkx as nx
from networkx import Graph
from itertools import combinations
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from DFTWorkflow.pitchingATent import convertCanonical
def listInputs(prompt:str):
    while True:
        partitionInput = input(prompt)
        partitionList = [part.strip() for part in partitionInput.split(",")]
        if len(partitionList) == 0 or any(part == ' ' for part in partitionList):
            print("Invalid Input. Enter at least 1 comma-seperated string")   
        else:
            print("partition List:", partitionList)
            conf = int(input("Confirm your strings entry by Typing 1:"))
            if conf == 1:
                break
    return partitionList
def getPartitions(chemistryType: str):
    numParts = int(input(f"Enter the number of partitions you want for Chem Type {chemistryType}: "))
    i = 1
    partitionDict = {}
    while i != numParts:
        promptChemistries = "Enter the string representations for partition" + str(i) + "(Ex: Mn,Manganese,Jacobsen) Enter " + str(chemistryType) + " to partition the chemistry of interest without any added catalyst or reagent"
        partitionList = listInputs(promptChemistries)
        if partitionList == chemistryType:
            partitionDict["NoCats"] = ["empty"]
            i += 1
        else:
            key = min(partitionList)
            partitionDict[key] = partitionList
            i += 1

    return partitionDict
def partitionDF(df , groupVar:str):
    dfList = []
    parameters = np.array(df[groupVar])
    parameters = np.unique(parameters)
    
    for i in range (len(parameters)):
        parm = parameters[i]
        smallDF = df[df[groupVar] == parm].copy()
        smallDF.reset_index(drop=True, inplace=True)
        dfList.append(smallDF)
    return dfList
def dataframeDivide(extractingCols , reagentList, reagentSplits, dataframeMAST, chemistry):
    refinedDF = dataframeMAST[dataframeMAST[extractingCols[2]].notna()].copy()
    refinedDF = refinedDF[refinedDF[extractingCols[1]].notna()].copy()
    refinedDF = refinedDF.drop(columns=[col for col in refinedDF.columns if col not in extractingCols])
    dfDict = {}
    substrateDF = pd.DataFrame(columns=refinedDF.columns)

    for index, row in refinedDF.iterrows():
        reagent = str(row[extractingCols[3]])
        if any(sub.lower() in reagent.lower() for sub in reagentList):
            substrateDF = pd.concat([substrateDF, pd.DataFrame([row])], ignore_index=True)
    substrateList = partitionDF(substrateDF , extractingCols[0] )
    finalColumns = [str(extractingCols[0]) , "SMILES" , "Yield" , "Reagent" , str(extractingCols[4]) ,  str(extractingCols[5]), str(extractingCols[6])]
    for reagent in list(reagentSplits.key()):
        reagentDF = pd.DataFrame(columns=finalColumns.columns)
        chemStrs = list(reagentSplits[reagent])
        if reagent == "NoCats":
            chemStr = chemistry + "with NoCats"
            #return only empty catalysts and empty clean reagents
            allChemistries = [s for sublist in reagentSplits.values() for s in sublist]

            for substrateDF in substrateList:
                dfNoCat = substrateDF[substrateDF[extractingCols[4]].isna()]
                if len(dfNoCat) != 0:
                    yieldList = []
                    for ind , row in dfNoCat.iterrows():
                        reagents = row[extractingCols[3]]
                        if not(any(sub.lower() in reagents.lower() for sub in allChemistries)):
                            yieldList.append(float(row[extractingCols[2]]))
                    if len(yieldList) != 0:
                        yieldMast = np.mean(yieldList)
                        smilesString = str(row[extractingCols[1]].split(">>")[0])
                        canonical = convertCanonical(smilesString)
                        newRow = pd.DataFrame([{str(extractingCols[0]): row[extractingCols[0]], 'SMILES': canonical, 'Yield' : yieldMast ,"Reagent" : chemStr , str(extractingCols[5]) : row[extractingCols[5]] ,str(extractingCols[6]) : row[extractingCols[6]]  }])

                        reagentDF= pd.concat([reagentDF, newRow], ignore_index=True)
        else:           
            chemStr = chemistry + " + "  + str(reagent)
            for substrateDF in substrateList:
                
                for index, row in substrateDF.iterrows():
                    reagents = str(row[extractingCols[3]])
                    catalysts = str(row[extractingCols[4]])
                    yieldList = []
                    if any(chem.lower() in reagents.lower() for chem in chemStrs) or any(chem.lower() in catalysts.lower() for chem in chemStrs):
                        yieldList.append(float(row[extractingCols[2]]))

        
    return dfDict
def main(fileDir:str , columnFile: str):
    df = pd.read_excel(fileDir , engine="openpyxl")
    if os.path.isfile(columnFile):
        with open(columnFile , "r") as file:
            headers = [col.strip() for col in file.readline().split(",")] #Reaxys ID number, SMILES of Reaction , Yield (or other numeric reaction assessment), Reagent , Catalyst, Referene , Reaxys Link
        
    chemistryList = listInputs("Enter the string representations for the Reagent type (Ex: NaOCl,bleach,sodium hypochlorite)")
    chemistryType = min(chemistryList)
    reagentPartitions = getPartitions(chemistryType)

    dfList = dataframeDivide(headers , chemistryList , reagentPartitions , masterDF)

    return dfList

if __name__ == "__main__":
    masterDF = sys.argv[1]
    columnFile = sys.argv[2]
    main()