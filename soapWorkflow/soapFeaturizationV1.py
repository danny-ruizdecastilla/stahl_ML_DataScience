#Danny Ruiz de Castilla | 09.10.2025
import numpy as np 
import os
import sys
import glob
from rdkit import Chem
from rdkit.Chem import AllChem
from dscribe.descriptors import SOAP
from ase import Atoms
import sparse
import re
import math
import pandas as pd
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from breadthFirstSearch.radialBasedCorrelation import getCC
from DFTWorkflow.conformerSearching.smiles_to_xyz import write_xyz
from dimensionalityReduction.reactivityFeatures import boxGen
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
def trackIdxSMiles_Coords(smile, molec, idxList):
    try:
        m = Chem.AddHs(molec)
    except Exception as e:
        print(f"ERROR: failed to prepare molecule from SMILES {smile}: {e}")
        return {}
    try:
        AllChem.EmbedMolecule(m)
    except Exception as e:
        print(f"ERROR: could not embed 3D coordinates for {smile}: {e}")
        return {}
    conf = m.GetConformer()
    coordHash = {}
    for idx in idxList:
        try:
            pos = conf.GetAtomPosition(idx)
            coordHash[idx] = [pos.x, pos.y, pos.z]
        except Exception as e:
            print(f"Warning: could not extract coordinates for atom {idx}: {e}")

    return coordHash
def soapFromxyz(dfDir , xyzDir):
    if not os.path.exists(xyzDir):
        os.makedirs(xyzDir)
    df = pd.read_csv(dfDir)
    colList = list(df.columns)
    colBox = boxGen(colList)
    while True:
        smileIdx = input(f"Here are the columns provided in the dataframe: {dfDir} Enter the number for the column corresponding to the SMILES strings {colBox}")
        try:
            smilesIdx = int(smileIdx)
            smilesCol = colList[smilesIdx]
            break
        except:
            print("Error: please enter an appropriate integer")
    while True:
        idIdx = input(f"Here are the columns provided in the dataframe: {dfDir} Enter the number for the column corresponding to the ID strings {colBox}")
        try:
            idIdx = int(idIdx)
            idCol = colList[idIdx]
            break
        except:
            print("Error: please enter an appropriate integer")
    alkeneHash = {}
    for _, row in df.iterrows():
        smiles = row[smilesCol]
        cc , mol = getCC(smiles)
        c1Idx = cc[0]
        c2Idx = cc[1]
        id = row[idCol]
        idxList = []
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()  
            idxList.append(idx)
        if smiles not in list(alkeneHash.keys()):
            xyzFile = xyzDir + "/" + str(id) + ".xyz"
            coordHash = trackIdxSMiles_Coords(smiles , mol , idxList)
            atomSymbols = []
            atomCoordinates = []
            for atomIdx , coords in coordHash.items():
                atom = mol.GetAtomWithIdx(atomIdx)
                symbol = atom.GetSymbol()
                atomSymbols.append(str(symbol))
                atomCoordinates.append(coords)
                if atomIdx == c1Idx:
                    c1_x = coords[0]
                    c1_y = coords[1]
                    c1_z = coords[2]
                elif atomIdx == c2Idx:
                    c2_x = coords[0]
                    c2_y = coords[1]
                    c2_z = coords[2]
            alkeneHash[smiles] = {"ID" : id , "file" : xyzFile, "c1Coords" : [c1_x, c1_y, c1_z] , "c2Coords" : [c2_x, c2_y, c2_z] ,}
            write_xyz(atomSymbols, atomCoordinates, xyzFile)
    dfMASt = pd.DataFrame()
    
        
            
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
        xyzDir = input(f"Enter the directory where the .xyz files for {smilesCSV} are contained or are going to be stored")
        alkeneSoapFeatures = soapFromxyz(smilesCSV , xyzDir)
