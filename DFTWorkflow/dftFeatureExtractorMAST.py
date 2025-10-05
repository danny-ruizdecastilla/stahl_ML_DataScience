#Master Feature Extractor
import os 
import sys
import glob
import tkinter
from pathlib import Path
import pandas as pd 
import numpy as np
import chemdraw
import base64
from rdkit import Chem
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from DFTWorkflow.cleanLogs import basicTerm
from DFTWorkflow.ionComGenerator import locateinLog
from DFTWorkflow.fukuiGenerator.fukuiExtractorV1 import  getBoltzmannWeightsGauss
from DFTWorkflow.AlkeneFeatures.alkeneNBOExtract import alkeneNBOExtractor
from dimensionalityReduction.reactivityFeatures import boxGen
from breadthFirstSearch.radialBasedCorrelation import getCC
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
def extractShiftsByIdx(logFile: str, extract1:str, extract2:str,location1:str ,location2, idxList  , m , b):
    lowerInd = locateinLog(logFile , extract1, location1)
    upperInd = locateinLog(logFile , extract2, location2 )
    if upperInd == "Poison" or lowerInd == "Poison":
        raise ValueError(f"The Log file {logFile} did not terminate properly")
    atomShifts = {}
    with open(logFile , "r") as f:
        for idx, line in enumerate(f):
            if idx > int(lowerInd)  and idx < int(upperInd):
                
                if "   Isotropic =  " in line:
                    atomInd = line.split("   Isotropic =  ")[0].strip()
                    atomIdx = int(atomInd.split(" ")[0].strip())
                    if atomIdx in idxList:
                        shiftsTots = line.split("   Isotropic =  ")[-1].strip()
                        shifts = shiftsTots.split("Anisotropy")[0]
                        atomName = str(atomInd.split(" ")[-1].strip())
                        atomShifts[atomIdx] = [ atomName,(b - float(shifts.strip())) / m]
    return atomShifts  
def getAlkenes(substratesHash , smilesHash , featureList):

    for id , smilesStr in smilesHash.items():
        substrateHash = {}
        cc , molec = getCC(smilesStr)
        conformerFiles = substratesHash[id]
        boltzmannDF = getBoltzmannWeightsGauss(conformerFiles, 298, "electronic")
        if "C13_shift" in featureList:
            C13_C1 = []
            C13_C2 = []
            
            for name in list(boltzmannDF["logID"]):
                fileStr = f"{name}.log"
                conformer  = next((f for f in conformerFiles if fileStr in f.name), None)
                alkeneHash = extractShiftsByIdx(conformer, "SCF GIAO Magnetic shielding tensor (ppm):", "Eigenvalues:  ", "earliest", "latest", cc , 1.0081 , 195.6683)
                if len(C13_C1) == 0:
                    #C1 and C2 not yet defined 
                    atom0 = alkeneHash[cc[0]][1]
                    atom1 = alkeneHash[cc[1]][1]
                    if atom0 >= atom1:#C1 is smallest shift, C2 is largest shift
                        C1 = cc[1]
                        C2 = cc[0]
                    else:
                        C1 = cc[0]
                        C2 = cc[1]
                C13_C1.append(alkeneHash[C1][1])
                C13_C2.append(alkeneHash[C2][1])
            weights = boltzmannDF["boltzWeights"]
            C13_Cmx = ( (C13_C2 * weights).sum()    / weights.sum()  )
            C13_Cmn = ( (C13_C1 * weights).sum()    / weights.sum()  )
            C13_delta = C13_Cmx - C13_Cmn
            C13_Mean = np.mean([C13_Cmx,C13_Cmn])
            featureList = featureList.remove("C13_shift")
        for feature in featureList:
            if feature == "NBO7":
                nboHash = alkeneNBOExtractor(conformerFiles , C1 , C2 , "electronic", logNameMAST , smilesStr)

def compartmentalization(logDir , outputDir , substrateFile):
    substrateDF = pd.read_csv(substrateFile)
    cols = list(substrateDF.columns())
    boxCols = boxGen(cols)
    while True:
        colId = input(f"Here are the columns for the file {substrateFile}\n {boxCols}\n Enter the ID number for the column corresponding to the SMILES strings").strip()
        try:
            smilesId = int(colId)
            smilesCol = cols[smilesId]
            smilesMAST = substrateDF[smilesCol]
            break
        except:
            print("Try again, enter an integer")
    while True:
        colId = input(f"Here are the columns for the file {substrateFile}\n {boxCols}\n Enter the ID number for the column corresponding to the substrate ID strings").strip()
        try:
            substrateId = int(colId)
            subId = cols[substrateId]
            idMAST = substrateDF[subId]
            break
        except:
            print("Try again, enter an integer")
    


    logPaths = Path(logDir) 
    logFiles = logPaths.glob('*.log')
    fileSplit = input(f"{logFiles[0]} Enter the string iteral that seperates the common name with the conf. type : ")

    substrateHash = {}
    smilesHash = {}
    for log in logFiles:
        termError = basicTerm(log , "Error termination" , "Normal termination")
        if not termError:
            substrates = list(substrateHash.keys())
            fileID = log.name.split(fileSplit)[0]
            if fileID in substrates:
                substrateHash[fileID].append(log)
            else:
                substrateHash[fileID] = [log]
                smilesHash[fileID] = smilesMAST[idMAST.index(fileID)]
    extractNum = input(f"Enter the number corresponding to which substructre you want to extract information from:\n [0] Alkenes\n")
    if int(extractNum) == 0:
        localStrs = ["C13_shift" , "NBO7" , "fukuiParameters" , "%Vbur" , "EvsZ" , "Dist." ]
        localDescriptorsInput = boxGen(localStrs)
        featureList = listInputs(f"Enter the indexes corresponding to the features you would like to extract\n{localDescriptorsInput}")
        featuresMAST = []
        for idx in featureList:
            feature = featureList[idx]
            featuresMAST.append(feature)
        substratesMAST = getAlkenes(substrateHash , smilesHash , featuresMAST)




    

if __name__ == "__main__":
    logDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    substrateCSV = str(sys.argv[3])
    compartmentalization(logDir , outputDir, substrateCSV)