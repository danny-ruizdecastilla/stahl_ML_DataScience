#Danny Ruiz de Castilla | 09.10.2025
import numpy as np 
import os
import sys
import glob
from dscribe.descriptors import SOAP
from ase import Atoms
import sparse
import re
import math
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from breadthFirstSearch.radialBasedCorrelation import getCC
def distance_3d_np(point1: list[float], point2: list[float]) -> float:
    p1 = np.array(point1, dtype=float)
    p2 = np.array(point2, dtype=float)
    if p1.shape != (3,) or p2.shape != (3,):
        raise ValueError("Both points must be 3-dimensional")
    return np.linalg.norm(p1 - p2)
def allUniquePairs(lst): #returns all unique combinations in a list of elements
    pairs = []
    n = len(lst)
    for i in range(n):
        for j in range(i + 1, n):
            print(j)
            pairs.append((lst[i], lst[j]))
    return pairs
def extractAlkenefromCoords(coordHash):
    carbons = {}
    for key , val in coordHash.items():
        atom = val[0]
        coords = val[2:5]
        if atom == "C":
            atomIdx = list(coordHash.keys()).index(key)
            carbons[atomIdx] = coords
    carbonAtoms = list(carbons.keys())
    allcarbonPairs = allUniquePairs(carbonAtoms)
    for pair in allcarbonPairs:
        coords1 = carbons[pair[0]]
        coords2 = carbons[pair[1]]
        euclidDist = distance_3d_np(coords1 , coords2)

if __name__ == "__main__":
    smilesCSV = str(sys.argv[1])
    while True:
        coordType = input(f"Here are the options for importing coordinates\n\n[1]  From .xyz file\n\n[2] From crest_conformer\n\n[3] From DFT log output file\n Select the integer that corresponds to the output file you are importing from").strip()
        optionList = ["1" , "2" , "3"]
        if coordType in optionList:
            break
        else:
            print("You must pick either 1, 2, or 3")
    if coordType == "1":
        
