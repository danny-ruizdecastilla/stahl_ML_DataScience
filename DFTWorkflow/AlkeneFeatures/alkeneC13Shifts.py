import os 
import sys
import glob
import shutil
import re
import pandas as pd 
import numpy as np
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(parentDir)
from breadthFirstSearch.radialBasedCorrelation import getCC
from DFTWorkflow.featureMaping import createCSV
from DFTWorkflow.pitchingATent import convertCanonical 
from DFTWorkflow.fukuiGenerator.addSCF import readDat
def main(substratesDir:str,  csvStr , theoryStr):
    substratePaths = glob.glob(substratesDir + "/*/identification.dat")
    substrateDF = pd.DataFrame()
    for i , path in enumerate(substratePaths):
        identifications = readDat(path)
        smiles = identifications[-1]
        CC , molec = getCC(smiles)
        C1 = CC[0] + 1
        C2 = CC[1] + 1
        #print(path , CC)
        mainPath = path.split("identification.dat")[0]
        nmrData = glob.glob(mainPath + f"/*{csvStr}*.csv")[0]
        nmrDF = pd.read_csv(nmrData)
        for _ , row in nmrDF.iterrows():
            atom = row["atom ID"]
            if int(atom) == int(C1):
                c1NMR = row["NMR Shift"]
            elif int(atom) == int(C2):
                c2NMR = row["NMR Shift"]
        cMax = max([c1NMR , c2NMR])
        cMin = min([c1NMR , c2NMR]) 
        cMean = np.mean([c1NMR , c2NMR])
        cDelta = abs(cMax - cMin)
        canonical = convertCanonical(smiles)
        newRow = pd.DataFrame([{
            "SMILES": smiles, "Canonicals": canonical, "Cmax_NMR_" + str(theoryStr) : cMax , 
            "Cmin_NMR_" + str(theoryStr) : cMin , theoryStr + "_NMR_Mean" : cMean , 
            theoryStr + "_Delta_CAlk" : cDelta
        }])
        substrateDF = pd.concat([substrateDF , newRow] , ignore_index=True)
    return substrateDF
if __name__ == "__main__":
    resultsDir = str(sys.argv[1])
    saveStr = str(sys.argv[2])
    csvStr = str(sys.argv[3])
    charge_theory = str(sys.argv[4])
    substrateDF = main(resultsDir, csvStr , charge_theory)
    createCSV(substrateDF , resultsDir , saveStr)