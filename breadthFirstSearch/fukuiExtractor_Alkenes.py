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
    return C1 , C2

def main(substratesDir:str, contactDist: int  ):

    substratePaths = glob.glob(substratesDir + "/*/identification.dat")
    substrateDF = pd.DataFrame(columns = ["SMILES" , "Canonicals" ,"C1_fuk_neg" , "C1_fuk_pos" ,
                                           "C1_fuk_neut" ,"C2_fuk_neg" , "C2_fuk_pos" , "C2_fuk_neut" 
                                             ,"delta_fuk_neg" , "delta_fuk_pos" , "delta_fuk_neut"  ])
    for i , path in enumerate(substratePaths):
        identifications = readDat(path)
        smiles = identifications[-1]

        C1 , C2 = getC1C2(smiles , contactDist)
        C1 += 1 #switch to log file indexing 
        C2 += 1
        mainPath = path.split("identification.dat")[0]
        fukuiData = glob.glob(mainPath + "/*.csv")[0]
        fukuiPD = pd.read_csv(fukuiData)
        C1List = C2List = None
        for _, row in fukuiPD.iterrows():
            atom = int(row["atoms"])
            f_neg, f_pos, f_neut = row["f_neg"], row["f_pos"], row["f_neut"]

            if atom == C1:
                C1List = [f_neg, f_pos, f_neut]
            elif atom == C2:
                C2List = [f_neg, f_pos, f_neut]
        if C1List is None or C2List is None:
            raise ValueError("C1 or C2 not found in Fukui data")

        deltaFuk = [abs(a - b) for a, b in zip(C1List, C2List)]
        canonical = convertCanonical(smiles)
        newRow = pd.DataFrame([{
            "SMILES": smiles,
            "Canonicals": canonical,
            "C1_fuk_neg": C1List[0],
            "C1_fuk_pos": C1List[1],
            "C1_fuk_neut": C1List[2],
            "C2_fuk_neg": C2List[0],
            "C2_fuk_pos": C2List[1],
            "C2_fuk_neut": C2List[2],
            "delta_fuk_neg": deltaFuk[0],
            "delta_fuk_pos": deltaFuk[1],
            "delta_fuk_neut": deltaFuk[2]
        }])
        substrateDF = pd.concat([substrateDF, newRow], ignore_index=True)
    return substrateDF
if __name__ == "__main__":
    fukuiOutputs = str(sys.argv[1])
    contactDist = int(sys.argv[2])
    saveDir = str(sys.argv[3])
    saveStr = str(sys.argv[4])
    substrateDF = main(fukuiOutputs , contactDist)
    createCSV(substrateDF , saveDir , saveStr)
