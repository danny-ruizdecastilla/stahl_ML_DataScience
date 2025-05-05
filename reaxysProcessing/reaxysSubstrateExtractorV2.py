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
def getPartitions(chemistryType: str):
    numParts = int(input(f"Enter the number of partitions you want for Chem Type {chemistryType}: "))
    i = 1
    partitionDict = {}
    while i != numParts:
        while True:
            partitionList = list(input(f"Enter the string representations for partition {i} (Ex: Mn,Manganese,Jacobsen )"))
        key = min(partitionList)
        partitionDict[key] = partitionList
        i += 1

    return partitionDict

def main(fileDir:str):
    masterDF = pd.read_csv(fileDir)
    chemistryType = list(input(f"Enter the string representations for reagent type (Ex: NaOCl,bleach,sodium hypochlorite)").split(","))






if __name__ == "__main__":
    masterDF = sys.argv[1]
    columnFile = sys.argv[2]
    main()