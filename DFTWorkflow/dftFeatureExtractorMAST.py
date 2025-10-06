#Master Feature Extractor
import os 
import sys
import glob
import tkinter
import pandas as pd 
import numpy as np
import chemdraw
import base64
from pathlib import Path
from morfeus import BuriedVolume
from rdkit import Chem
from networkx import Graph
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from DFTWorkflow.cleanLogs import basicTerm
from DFTWorkflow.ionComGenerator import locateinLog , getAtomCoords
from DFTWorkflow.fukuiGenerator.fukuiExtractorV1 import  getBoltzmannWeightsGauss
from DFTWorkflow.AlkeneFeatures.alkeneNBOExtract import alkeneNBOExtractor
from DFTWorkflow.AlkeneFeatures.alkeneFukuiExtractor import getAlkeneFukuiFunctions
from DFTWorkflow.AlkeneFeatures.alkeneSubstitution import eVszAlkenes
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
def hCount(molec, wildCards):
    hydrogens = 0
    for wildcard in wildCards:
        atom = molec.GetAtomWithIdx(wildcard)
        atom1 = str(atom.GetSymbol())
        if atom1 == "H":
            hydrogens += 1
    return hydrogens
def getAlkenes(substratesHash , smilesHash , featureList , **kwargs):
    if kwargs.get("anions") is None:
        pass
    else:
        anionFiles = str(kwargs["anions"])
    if kwargs.get("cations") is None:
        pass
    else:
        cationFiles = str(kwargs["cations"])
    featuresMASTDF = pd.DataFrame()
    for id , smilesStr in smilesHash.items():
        hashList = []
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
                        C1 = cc[1] + 1
                        C2 = cc[0] + 1
                    else:
                        C1 = cc[0] + 1
                        C2 = cc[1] + 1
                C13_C1.append(alkeneHash[C1][1])
                C13_C2.append(alkeneHash[C2][1])
            weights = boltzmannDF["boltzWeights"]
            C13_Cmx = ( (C13_C2 * weights).sum()    / weights.sum()  )
            C13_Cmn = ( (C13_C1 * weights).sum()    / weights.sum()  )
            C13_delta = C13_Cmx - C13_Cmn
            C13_Mean = np.mean([C13_Cmx,C13_Cmn])
            featureList = featureList.remove("C13_shift")
            c13Hash = {"C13_Cmx" : C13_Cmx , "C13_Cmn" : C13_Cmn , "C13_delta" : C13_delta , "C13_Mean" : C13_Mean}
            hashList.append(c13Hash)
        for feature in featureList:
            if feature == "NBO7":
                nboHash = alkeneNBOExtractor(conformerFiles , C1 , C2 , "electronic", id , smilesStr)
                hashList.append(nboHash)
            elif feature == "fukuiParameters":
                
                fukuiHash = getAlkeneFukuiFunctions(conformerFiles , cationFiles , anionFiles , boltzmannDF , C1 , C2 , "NBO7" , id , smilesStr )
                hashList.append(fukuiHash)
            elif feature == "%Vbur":
                Vbur_Cmn_2 = []
                Vbur_Cmn_3 = []
                Vbur_Cmx_2 = []
                Vbur_Cmx_3 = []
                
                for name in list(boltzmannDF["logID"]):
                    fileStr = f"{name}.log"
                    conformer  = next((f for f in conformerFiles if fileStr in f.name), None)

                    coordHash = getAtomCoords(str(conformer) , "GINC-COMPUTE" , 5 )
                    elements = []
                    coordinates = []
                    for _, coords in coordHash.items:
                        elements.append(str(coords[0]))
                        coordinates.append(np.array(coords[2:5]))
                    vburC1_2 = BuriedVolume(elements, coordinates, int(C1-1), include_hs=True, radius=2.0)
                    vburC1_3 = BuriedVolume(elements, coordinates, int(C1-1), include_hs=True, radius=3.5)
                    vburC2_2 = BuriedVolume(elements, coordinates, int(C2-1), include_hs=True, radius=2.0)
                    vburC2_3 = BuriedVolume(elements, coordinates, int(C2-1), include_hs=True, radius=3.5)
                    Vbur_Cmn_2.append(vburC1_2)
                    Vbur_Cmx_2.append(vburC2_2)
                    Vbur_Cmn_3.append(vburC1_3)
                    Vbur_Cmx_3.append(vburC2_3)

                weights = boltzmannDF["boltzWeights"]
                Vbur_Cmx_2Ang = ( (Vbur_Cmx_2 * weights).sum()    / weights.sum()  )
                Vbur_Cmn_2Ang= ( (Vbur_Cmn_2 * weights).sum()    / weights.sum()  )    
                Vbur_Cmx_3Ang = ( (Vbur_Cmx_3 * weights).sum()    / weights.sum()  )
                Vbur_Cmn_3Ang= ( (Vbur_Cmn_3 * weights).sum()    / weights.sum()  )               

                Vbur_2Ang_delta = Vbur_Cmx_2Ang - Vbur_Cmn_2Ang
                Vbur_3Ang_delta = Vbur_Cmx_3Ang - Vbur_Cmn_3Ang
                Vbur_2Ang_mean =  np.mean([Vbur_Cmx_2Ang, Vbur_Cmn_2Ang])
                Vbur_3Ang_mean = np.mean([Vbur_Cmx_3Ang, Vbur_Cmn_3Ang])
                BurVolHash = {"2Ang_Cmx" :Vbur_Cmx_2Ang , "2Ang_Cmn" :Vbur_Cmn_2Ang ,
                              "3Ang_Cmx" :Vbur_Cmx_3Ang , "3Ang_Cmn" :Vbur_Cmn_3Ang ,
                               "2Ang_delta" :Vbur_2Ang_delta , "3Ang_delta" :Vbur_3Ang_delta ,
                              "3Ang_mean" :Vbur_3Ang_mean , "2Ang_mean" :Vbur_2Ang_mean }
                hashList.append(BurVolHash)
            elif feature == "EvsZ":
                conformer = conformerFiles[0]
                coordHash = getAtomCoords(str(conformer) , "GINC-COMPUTE" , 5 )
                g = Graph()
                #print("CC" , CC)
                for bond in molec.GetBonds():
                    start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                    g.add_edge(start, end)
                CminNeighbors = list(g.neighbors(int(C1-1))).remove(int(C2-1))
                CminHcount = hCount(molec, CminNeighbors)
                CmaxNeighbors = list(g.neighbors(int(C2-1))).remove(int(C1-1))
                CmaxHcount = hCount(molec, CmaxNeighbors)

                hCount = CmaxHcount + CminHcount
                if hCount <=2: #Easy check for E vs Z 
                    alkeneType = eVszAlkenes(molec , g ,[int(C1-1) , int(C2-1)] ,  coordHash)# 0 if neither E or Z, -1 if Z, 1 if E
                else:
                    alkeneType = 0
                alkeneHash = {"HCount" : hCount , "EvsZ" : alkeneType}
                hashList.append(alkeneHash)
            elif feature == "Dist.":
                distList = []
                for name in list(boltzmannDF["logID"]):
                    fileStr = f"{name}.log"
                    conformer  = next((f for f in conformerFiles if fileStr in f.name), None)

                    coordHash = getAtomCoords(str(conformer) , "GINC-COMPUTE" , 5 )
                    idx = 0
                    for _ , coords in coordHash.items():
                        idx += 1
                        if idx == C1:
                            c1_coords = coords[2:5]
                        elif idx == C2:
                            c2_coords = coords[2:5]
                    dist1 = np.linalg.norm(c1_coords - c2_coords)
                    distList.append(dist1)
                distMAST = ( (distList * weights).sum()    / weights.sum()  )  
                distHash = {"C1C2Dist" : distMAST} 
                hashList.append(distHash)
        masterDict = {}
        for d in hashList:
            masterDict.update(d)
        df_row = pd.DataFrame([masterDict]) 
        featuresMASTDF = pd.concat([featuresMASTDF, df_row], ignore_index=True)
    return featuresMASTDF
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
        if "fukuiParameters" in featuresMAST:
            cationDir = input(f"Enter the directory corresponding to the cation files for {logDir}: ")
            cationPaths = Path(cationDir)
            anionDir = input(f"Enter the directory corresponding to the cation files for {logDir}: ")
            anionPaths = Path(anionDir)
            anionFiles = anionPaths.glob("*.log")
            cationFiles = cationPaths.glob("*.log")
            substratesMAST = getAlkenes(substrateHash , smilesHash , featuresMAST , anions = anionFiles, cations = cationFiles)
            outputFile = Path(outputDir) / "alkeneFeaturesMAST.csv"
            substratesMAST.to_csv(outputFile , index=False )
        else:

            substratesMAST = getAlkenes(substrateHash , smilesHash , featuresMAST)
            outputFile = Path(outputDir) / "alkeneFeaturesMAST.csv"
            substratesMAST.to_csv(outputFile , index=False )




    

if __name__ == "__main__":
    logDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    substrateCSV = str(sys.argv[3])
    compartmentalization(logDir , outputDir, substrateCSV)