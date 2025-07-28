import os 
import sys
import glob
import shutil
import re
import pandas as pd 
import numpy as np
from rdkit import Chem
from motifExtractorGifGen import motifExtract

parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from DFTWorkflow.fukuiGenerator.addSCF import readDat
from DFTWorkflow.pitchingATent import convertCanonical 
from DFTWorkflow.featureMaping import createCSV
def dict2Dat(dict1, saveDir, saveStr):
    os.makedirs(saveDir, exist_ok=True)
    filePath = os.path.join(saveDir, saveStr + ".dat")
    
    if not os.path.exists(filePath):
        with open(filePath, "w") as file:
            for key, value in dict1.items():
                value_str = ",".join(str(a) for a in value)
                file.write(f"{key}:{value_str}\n")
def getMaxCoordinate( allPaths , coordinate , initList):
    coordinateList = []
    for path in allPaths:
        newPath = [atom for atom in path if atom not in initList]
        if len(newPath) > 0:
            coordinateList.append(int(newPath[coordinate]))
    coordinateList = list(set(coordinateList))
    return coordinateList

def getC1C2(smiles, cutDist):
    ptable = Chem.GetPeriodicTable()
    contactsDict , smilesEdgesDict , allPaths , g, molec = motifExtract(smiles, cutDist)
    weightDict = {}
    for key, contactList in contactsDict.items():
        molWeight = 0
        for idx in contactList:
            atom = molec.GetAtomWithIdx(idx)
            atomicMass = ptable.GetAtomicWeight(atom.GetAtomicNum())
            molWeight += atomicMass
        weightDict[key] = molWeight
    atoms = list(weightDict.keys())
    if weightDict[atoms[0]] >= weightDict[atoms[1]]:
        C1 = atoms[0]
        C2 = atoms[1]
    else:
        C1 = atoms[1]
        C2 = atoms[0]
    wildList = getMaxCoordinate( allPaths,0,[C1 , C2]  )
    return C1 , C2 , wildList

def main(substratesDir:str, contactDist: int  , saveDir , fukuiResults):

    substratePaths = glob.glob(substratesDir + "/*/identification.dat")
    substrateDF = pd.DataFrame(columns = ["SMILES" , "Canonicals" ,"C1_fuk_neg" , "C1_fuk_pos" ,"C1_fuk_neut" ,
                                          "C2_fuk_neg" , "C2_fuk_pos" , "C2_fuk_neut" 
                                             ,"delta_fuk_neg" , "delta_fuk_pos" , "delta_fuk_neut"  , 
                                             "CMax_fuk_neg" , "CMax_fuk_pos" , "CMax_fuk_neut",
                                             "maxRel_fuk_neg", "maxRel_fuk_pos" , "maxRel_fuk_neut"  ])
    for i , path in enumerate(substratePaths):
        identifications = readDat(path)
        smiles = identifications[-1]

        C1 , C2 , wildList = getC1C2(smiles , contactDist)
        C1 += 1 #switch to log file indexing 
        C2 += 1
        wildList = [a + 1 for a in wildList]
        #print(identifications[0])
        print("C1" , C1 , "C2" , C2)
        wildDict = {}
        mainPath = path.split("identification.dat")[0]
        fukuiData = glob.glob(mainPath + f"/*{fukuiResults}*.csv")[0]
        #print(fukuiData)
        fukuiPD = pd.read_csv(fukuiData)
        C1List = C2List = None
        max_fneg = max_fpos = max_fneut = 0

        for _, row in fukuiPD.iterrows():
            atom = row["atom ID"]
            #print(atom)
            f_neg, f_pos, f_neut = row["f_neg"], row["f_pos"], row["f_neut"]
            if f_neg > max_fneg:
                max_fneg = f_neg
            if f_pos > max_fpos:
                max_fpos = f_pos
            if f_neut > max_fneut:
                max_fneut = f_neut
            if atom == C1:
                C1List = [f_neg, f_pos, f_neut]
            elif atom == C2:
                C2List = [f_neg, f_pos, f_neut]
            elif atom in wildList:
                wildDict[atom] = [f_neg, f_pos, f_neut]
        if C1List is None or C2List is None:
            raise ValueError(f"C1 or C2 not found in Fukui data for {path}")
        elif len(wildDict.keys()) == 0:
            raise ValueError(f"Wildcards not found in Fukui data for {path}")
        dict2Dat(wildDict, saveDir , "WildFukuiVals")
        wild0 = []
        wildPlus = []
        wildNeg = []
        for wild , fukuiList in wildDict.items():
            wild0.append(fukuiList[2])
            wildPlus.append(fukuiList[1])
            wildNeg.append(fukuiList[0])
        maxWild_0 = max(wild0)
        maxWild_neg = max(wildNeg)
        maxWild_pos = max(wildPlus)
        maxFukNeg = max([C1List[0] , C2List[0]])
        maxFukPos = max([C1List[1] , C2List[1]])
        maxFukNeut = max([C1List[2] , C2List[2]])
        maxRelNeg = abs(maxFukNeg - max_fneg)
        maxRelPos = abs(maxFukPos - max_fpos)
        maxRelNeut = abs(maxFukNeut - max_fneut)

        deltaFuk = [abs(a - b) for a, b in zip(C1List, C2List)]

        canonical = convertCanonical(smiles)
        newRow = pd.DataFrame([{
            "SMILES": smiles,
            "Canonicals": canonical,
            "C1_fuk_neg": C1List[0],"C1_fuk_pos": C1List[1],"C1_fuk_neut": C1List[2],
            "C2_fuk_neg": C2List[0],"C2_fuk_pos": C2List[1],"C2_fuk_neut": C2List[2],
            "delta_fuk_neg": deltaFuk[0],"delta_fuk_pos": deltaFuk[1],"delta_fuk_neut": deltaFuk[2], 
            "CMax_fuk_neg": maxFukNeg , "CMax_fuk_pos":maxFukPos , "CMax_fuk_neut": maxFukNeut,
            "maxRel_fuk_neg" : maxRelNeg, "maxRel_fuk_pos": maxRelPos , "maxRel_fuk_neut": maxRelNeut  
        }])
        substrateDF = pd.concat([substrateDF, newRow], ignore_index=True)
    return substrateDF
if __name__ == "__main__":
    fukuiOutputs = str(sys.argv[1])
    contactDist = int(sys.argv[2])
    saveDir = str(sys.argv[3])
    saveStr = str(sys.argv[4])
    fukuiResults = str(sys.argv[5])
    substrateDF = main(fukuiOutputs , contactDist , saveDir , fukuiResults)
    charge_theory = str(sys.argv[6])
    columns = substrateDF.columns.tolist()
    for i, col in enumerate(columns):
        if "fuk" in col:
            columns[i] = charge_theory + "_" + col
    substrateDF.columns = columns
    createCSV(substrateDF , saveDir , saveStr)
