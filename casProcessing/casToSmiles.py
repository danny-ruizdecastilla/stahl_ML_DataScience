import pandas as pd
import os 
import sys
import glob
import json
import numpy as np
import time
from pathlib import Path
import cirpy
#import requests
#from bs4 import BeautifulSoup
from rdkit import Chem
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from dimensionalityReduction.reactivityFeatures import boxGen
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
#Go from CAS to SMILES lists using beautifl soup 
#takes in csv of CAS numbers
#Danny Ruiz de Castilla | 11.10.2024

'''def pubChemScrapper():
    First cas result will be a direct search inside sigma, need to click on first result to access smiles webpage 
    '''
def testIsAlkene(smiles , alkCount):
    #Checks if a smiles string contains alkCount amount of c=c bonds 
    try:
        molec = Chem.MolFromSmiles(smiles)
        molecH = Chem.AddHs(molec)
        doubleCount = 0
        for bond in molecH.GetBonds():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6 and bond.GetBondType() == Chem.BondType.DOUBLE:
                c1 = bond.GetBeginAtomIdx()
                c2 = bond.GetEndAtomIdx()
                doubleCount +=1
        if doubleCount != alkCount:
            return False , ["Error" , "Error"] , molec , smiles
        else:
            return True , [c1 , c2] , molec , smiles


    except:
        #smiles string is wrong or corrupted or cant be kekulized 
        return False , ["Fatal" , "Fatal"] , 0 , smiles

def getSMILES(casList):
    casHash = {}
    for cas in casList:
        if "/" in cas:
            casParts = cas.split("/")
            newCas = ""
            for i in range(len(casParts)):
                part = casParts[i]
                if i == len(casParts)-1:
                    newCas += str(part)
                else:
                    newCas += (str(part) + "_")
                
                               
def main(chemistryFile):
    dfMAST = pd.read_csv(chemistryFile)
    colList = list(dfMAST.columns)
    colBox = boxGen(colList)
    casInt = int(input(f"Here are all the columns from the dataframe presented {colBox} select the index that matches with the columnn name for CAS numbers"))
    casStr = colList[casInt]
    dfReduced = dfMAST[dfMAST[casStr].notna()]
    casListRaw = set(list(dfReduced[casStr]))
    casHash = getSMILES(list(casListRaw))



    

    

if __name__ == "__main__":
    chemistryDir = str(sys.argv[1])
    main(chemistryDir)