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
from rapidfuzz import process
from rapidfuzz.distance import JaroWinkler
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
        partitionList = listInputs("Enter the string representations for partition" + str(i) + "(Ex: Mn,Manganese,Jacobsen)")
        key = min(partitionList)
        partitionDict[key] = partitionList
        i += 1

    return partitionDict
def dataframeDivide():
    dfList = []


    return dfList
def main(fileDir:str , columnFile: str):
    masterDF = pd.read_csv(fileDir)
    if os.path.isfile(columnFile):
        with open(columnFile , "r") as file:
            headers = file.readline().strip().split(",") #ID , Yield , Reagent , Catalyst, Referene , Reaxys Link
        
    chemistryList = listInputs("Enter the string representations for the Reagent type (Ex: NaOCl,bleach,sodium hypochlorite)")
    chemistryType = min(chemistryList)
    reagentPartitions = getPartitions(chemistryType)

    dfList = dataframeDivide(headers , chemistryList , reagentPartitions , masterDF)

    return dfList

if __name__ == "__main__":
    masterDF = sys.argv[1]
    columnFile = sys.argv[2]
    main()