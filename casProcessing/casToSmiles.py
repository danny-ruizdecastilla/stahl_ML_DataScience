import pandas as pd
#import os 
import sys
import glob
#import json
import time
from pathlib import Path
import cirpy
#import requests
#from bs4 import BeautifulSoup
from rdkit import Chem
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from dimensionalityReduction.reactivityFeatures import boxGen
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
            return False , smiles
        else:
            return True ,  smiles


    except:
        #smiles string is wrong or corrupted or cant be kekulized 
        return False , smiles

def getAlkeneSMILES(casList):
    casHash = {}
    for cas in casList:
        if "/" in cas:
            casParts = cas.split("/")
            newCas = "_".join(casParts[:-1]) + ("_" + casParts[-1])
        else:
            newCas = cas
        smilesStr = cirpy.resolve(newCas , 'smiles')
        if smilesStr == "None":
            '''
            smilesStr = pubChemScrapper(newCas)
            isAlkene , smiles = testIsAlkene(smilesStr , 1)
            if isAlkene:
                casHash[cas] = [smiles]
            else:
                pass
            '''
            pass
        else:
            isAlkene , smiles = testIsAlkene(smilesStr , 1)
            if isAlkene:
                casHash[cas] = [smiles]
            else:
                pass
    return casHash

def main(chemistryFile):
    chemPath = Path(chemistryFile)
    chemName = str(chemPath.root.split(".")[0])
    dfMAST = pd.read_csv(chemistryFile)
    colList = list(dfMAST.columns)
    colBox = boxGen(colList)
    casInt = int(input(f"Here are all the columns from the dataframe presented {colBox} select the index that matches with the columnn name for CAS numbers"))
    casStr = colList[casInt]
    dfReduced = dfMAST[dfMAST[casStr].notna()]
    casListRaw = set(list(dfReduced[casStr]))
    smilesHash = getAlkeneSMILES(list(casListRaw))
    alkeneCASNums = list(smilesHash.keys())
    alkeneDF = pd.DataFrame()
    for idx , row in dfReduced.iterrows():
        casNumber = row[casStr]
        if casNumber in alkeneCASNums:
            rowDict = row.to_dict() 
            rowDict["SMILES"] = str(smilesHash[casNumber])
            dfRow = pd.DataFrame([rowDict])
            alkeneDF = pd.concat([dfRow,alkeneDF] , ignore_index=True)
        else:
            continue
    return alkeneDF , chemName

if __name__ == "__main__":
    chemistryFile = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    alkeneDF , mainName = main(chemistryFile)
    alkeneDF.to_csv(outputDir + f"/{mainName}_withSMILES.csv")