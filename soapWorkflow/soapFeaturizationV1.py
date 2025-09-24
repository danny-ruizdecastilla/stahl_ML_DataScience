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
import math
import pandas as pd
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from breadthFirstSearch.radialBasedCorrelation import getCC
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
def trackIdxSMiles_Coords(smiles, molec, idxList):
    try:
        m = Chem.AddHs(molec)
    except Exception as e:
        print(f"ERROR: failed to prepare molecule from SMILES {smiles}: {e}")
        return {}
    try:
        AllChem.EmbedMolecule(m)
    except Exception as e:
        print(f"ERROR: could not embed 3D coordinates for {smiles}: {e}")
        return {}
    #print(smiles)
    try:
        AllChem.EmbedMolecule(m, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(m)
        conf = m.GetConformer()
    except Exception as e:
        print(f"ERROR: failed to find a conformer for {smiles}: {e}" )
        return {}
    coordHash = {}
    for idx in idxList:
        try:
            pos = conf.GetAtomPosition(idx)
            coordHash[idx] = [pos.x, pos.y, pos.z]
        except Exception as e:
            print(f"Warning: could not extract coordinates for atom {idx}: {e}")

    return coordHash
def soapFromxyz(dfDir):
    df = pd.read_csv(dfDir)
    colList = list(df.columns)
    colBox = boxGen(colList)
    while True:
        smileIdx = input(f"Here are the columns provided in the dataframe: {dfDir}\n{colBox}\n Enter the number for the column corresponding to the SMILES strings: ")
        try:
            smilesIdx = int(smileIdx)
            smilesCol = colList[smilesIdx]
            break
        except:
            print("Error: please enter an appropriate integer")
    while True:
        yieldIdx = input(f"Here are the columns provided in the dataframe: {dfDir}\n{colBox}\n Enter the number for the column corresponding to the reaction metric: ")
        try:
            idIdx = int(yieldIdx)
            idYield = colList[idIdx]
            break
        except:
            print("Error: please enter an appropriate integer")
    alkeneHash = {}
    masterAtoms = []
    for _ , row in df.iterrows():
        smiles = row[smilesCol]
        metric = row[idYield]
        cc , mol = getCC(smiles)
        if cc[0] == "Error":
            continue
        c1Idx = cc[0]
        c2Idx = cc[1]
        idxList = []
        for atom in mol.GetAtoms():
            idx = atom.GetIdx()  
            idxList.append(idx)
        if smiles not in list(alkeneHash.keys()):
            coordHash = trackIdxSMiles_Coords(smiles , mol , idxList)
            if len(list(coordHash.keys())) == 0:
                continue
            alkeneSymbols = []
            alkeneCoords = []
            aceIdx = 0
            for id , coords in coordHash.items():
                atom = mol.GetAtomWithIdx(id)
                symbol = atom.GetSymbol()
                if symbol not in masterAtoms:
                    masterAtoms.append(symbol)
                alkeneSymbols.append(symbol)
                alkeneCoords.append(np.array(coords))
                if id == c1Idx:
                    aceC1 = aceIdx
                elif id ==c2Idx:
                    aceC2 = aceIdx
                aceIdx += 1
            aceObj = Atoms(symbols=alkeneSymbols, positions=alkeneCoords)
            alkeneHash[smiles] = {"alkeneAce" : aceObj , "alkeneCenters" : [aceC1,aceC2] , "metric" : metric  }
    while True:
        rCut = input(f"Enter the radial cutoff for these SOAP descriptors in Angstroms:")
        try:
            rfloat = float(rCut)
            break
        except:
            print("Error: please enter an appropriate integer/float")
    while True:
        rho = input(f"Enter the number of basis functions for these SOAP descriptors in Angstroms:")
        try:
            rho = int(rho)
            break
        except:
            print("Error: please enter an appropriate integer")
    while True:
        l = input(f"Enter the degree of spherical harmonics for these SOAP descriptors in Angstroms:")
        try:
            l_ = int(l)
            break
        except:
            print("Error: please enter an appropriate integer")
    dfMAST = pd.DataFrame()
    print("allatoms" , len(masterAtoms))
    soap = SOAP(species=masterAtoms,r_cut=rfloat,n_max=rho,l_max=l_,)
    for alkene , info in alkeneHash.items():
        soapHash = {}
        soapHash["SMILES"] = alkene
        soapHash["Y"] = info["metric"]
        mol = info["alkeneAce"]
        alkeneID = info["alkeneCenters"]
        soapParameters = soap.create(mol , centers = alkeneID)
        for i in range(len(soapParameters)):
            alk = "C" + str(i+1)
            for j in range(len(soapParameters[i])):
                soapStr = alk + "_soap_" + str(j)
                soapVal = soapParameters[i][j]
                soapHash[soapStr] = soapVal
        dfMAST = pd.concat([dfMAST, pd.DataFrame([soapHash])], ignore_index=True)
    return dfMAST
            
if __name__ == "__main__":
    smilesCSV = str(sys.argv[1])
    while True:
        coordType = input(f"Here are the options for importing coordinates\n\n[1]  From .xyz file\n\n[2] From crest_conformer\n\n[3] From DFT log output file\n Select the integer that corresponds to the output file you are importing from: ").strip()
        optionList = ["1" , "2" , "3"]
        if coordType in optionList:
            break
        else:
            print("You must pick either 1, 2, or 3")
    if coordType == "1":
        outputDir = input(f"Enter the output directory for {smilesCSV}: ")
        if not os.path.exists(outputDir):
            os.makedirs(outputDir)
        alkeneSoapFeatures = soapFromxyz(smilesCSV)
        print(alkeneSoapFeatures.shape)
        alkeneSoapFeatures.to_csv(outputDir + "/soapFeatures.csv")
