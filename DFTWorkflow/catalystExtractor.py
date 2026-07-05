#Danny Ruiz de Castilla 
#07.04.2026
import os 
import sys
import pandas as pd 
import numpy as np
import chemdraw
import base64
import re
from pathlib import Path
from morfeus import SASA , Sterimol, BuriedVolume
from rdkit import Chem
from networkx import Graph
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
from DFTWorkflow.cleanLogs import basicTerm
from dimensionalityReduction.reactivityFeatures import boxGen
from DFTWorkflow.dftFeatureExtractorMAST import time_to_seconds , extractNBOOccupancies , getBoltzmannWeightsGauss , getAtomCoordsRobust , extractShiftsByIdx , getGlobalGreeks

def locateCatRows(logDir , splitStr , catalystFile):
    catalystDF = pd.read_csv(catalystFile)
    cols = list(catalystDF.columns)
    boxCols = boxGen(cols)
    while True:
        colId = input(f"Here are the columns for the file {catalystFile}\n {boxCols}\n Enter the ID number for the column corresponding to the SMILES strings: ").strip()
        try:
            smilesId = int(colId)
            smilesCol = cols[smilesId]
            smilesMAST = catalystDF[smilesCol]
            break
        except:
            print("Try again, enter an integer")
    while True:
        colId = input(f"Here are the columns for the file {catalystFile}\n {boxCols}\n Enter the ID number for the column corresponding to the catalyst ID strings: ").strip()
        try:
            catalystId = int(colId)
            subId = cols[catalystId]
            idMAST = catalystDF[subId]
            break
        except:
            print("Try again, enter an integer")
    #returns hash table of log paths that lead to unique catalysts with solvents 
    path = Path(logDir)
    logFiles = list(path.glob("/*/comFiles/logs/*.log"))
    catalystHash = {}
    smilesHash = {}
    for log in logFiles:
        termError = basicTerm(log , "Error termination" , "Normal termination")
        if not termError:
            catalysts = list(catalystHash.keys())
            fileID = log.name.split(splitStr)[0]
            solventID = Path(str(log).split("comFiles")[0]).name
            try:
                pos = idMAST[idMAST == fileID].index[0]
                if f"{fileID}_{solventID}" in catalysts:
                    catalystHash[f"{fileID}_{solventID}"].append(log)
                else:
                    catalystHash[f"{fileID}_{solventID}"] = [log]
                    pos = idMAST[idMAST == fileID].index[0]
                    smilesHash[fileID] = smilesMAST[pos]
            except:
                continue
    return catalystHash , smilesHash
def getChiralKetoneFeatures(catalystHash , smilesHash):
