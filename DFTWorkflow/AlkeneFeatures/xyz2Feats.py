#Danny Ruiz de Castilla
#Generate features from xyz files 
#MM Paper
import numpy as np
import pandas as pd 
from pathlib import Path
import sys
from rdkit import Chem
from itertools import combinations
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
def alkene2xyz(xyzFile):
    xyz = Path(xyzFile)
    atomHash = {}
    atomCoords = []
    carbonList = []
    with open(xyz, "r") as f:
        for lineNum, line in enumerate(f, start=1):
            if lineNum == 1:
                numAtoms = int(line.strip())
                atomHash["totAtoms"] = numAtoms
            elif "." in line: #atom coords
                atomIdx = int(lineNum - 3) 
                parts = line.strip().split()
                atomHash[atomIdx] = [str(parts[0])] + [float(part) for part in parts[1:4]]
                atomCoords.append(np.array([float(part) for part in parts[1:4]]))
                if str(parts[0]) == "C":
                    carbonList.append(atomIdx)
    coords = np.array(atomCoords)
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]

    distanceMatrix = np.linalg.norm(diff, axis=2)
    alkeneMax = 1.430417269
    alkeneMin = 1.308758104
    alkenes = []
    for i, j in combinations(carbonList, 2):
        d = distanceMatrix[i, j]
        if alkeneMin <= d <= alkeneMax:
            alkenes.append((i, j, d))
