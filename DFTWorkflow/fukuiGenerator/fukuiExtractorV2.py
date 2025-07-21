import os 
import sys
import glob
import shutil
import re
import pandas as pd 
import numpy as np

parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(parentDir)
from DFTWorkflow.featureMaping import createCSV
from DFTWorkflow.fukuiGenerator.fukuiExtractorV1 import  getBoltzmannWeightsGauss , extractChargesCHELPG_Mull_NBO_Hirsh
#Danny Ruiz de Castilla
#Generates Fukui maps for each molecule 
def extractChargesByDensity(log, densityStr):
    if densityStr == "NBO7":
        return extractChargesCHELPG_Mull_NBO_Hirsh(
            log, "Atom No    Charge        Core ", "* Total *",
            "earliest", "earliest", 1, -1, 1, 0)
    elif densityStr == "mull":
        return extractChargesCHELPG_Mull_NBO_Hirsh(
            log, " Mulliken charges", " Sum of Mulliken charges =",
            "earliest", "earliest", 1, 0, 0, 1)
    elif densityStr == "chelpg":
        return extractChargesCHELPG_Mull_NBO_Hirsh(
            log, " ESP charges:", " Sum of ESP charges =",
            "latest", "latest", 1, 0, 0, 1)
    elif densityStr == "hirsh":
        return extractChargesCHELPG_Mull_NBO_Hirsh(
            log, "              Q-H", "       Tot ",
            "earliest", "latest", 0, 0, 0, 1)
    else:
        raise ValueError(f"Unknown density string: {densityStr}")
def processIonLogs(logs, keyword, densityStr):
    result = {}
    for log in logs:
        try:
            name = os.path.basename(log).split(keyword)[0]
        except Exception:
            name = os.path.basename(log).split(".")[0]
        result[name] = extractChargesByDensity(log, densityStr)
    return result
def consolidatePaths(pathList):
    split1 = input(f"This is the root of all common files: {pathList[-1]} Enter the front split to find the common suffix: ")
    paths = []
    for log in pathList:
        name = log.split(split1)[0]
        #print(log)
        paths.append(name)
    substrates = list(set(paths))
    return sorted(substrates)
def extractFromKeys(dict1 , name1):
    return {key: value.copy() for key, value in dict1.items() if name1 in key}
def atomicBoltzmannConstruction(masterLogs , boltzWeights,  substratesDict , colList  ):
    substrates = substratesDict.keys()
    df = pd.DataFrame(columns = colList)
    for i , log in enumerate(masterLogs):
        weight = boltzWeights[i]
        substrate = [substrate for substrate in substrates if substrate in log][0]
        substrateVals = substratesDict[substrate]
        for key , val in substrateVals.items():
            atomID = int(key)
            atomType = str(val[0])
            scaledValue = weight* val[-1]
            if atomID not in df[colList[0]].values:
                new_row = {colList[0]: atomID, colList[1]: atomType, colList[-1]: scaledValue}
                df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)
            else:
                idx = df[df[colList[0]] == atomID].index[0]
                df.at[idx, colList[-1]] += scaledValue
    return df 

def main(logDir, cationDir, anionDir, substrateCSV, outputDir, densityStr , bltzmannStr):
    substrateScope = pd.read_csv(substrateCSV, encoding='utf-8')

    neutralLogs = glob.glob(os.path.join(logDir, "*.log"))
    cationLogs = glob.glob(os.path.join(cationDir, "*.log"))
    anionLogs = glob.glob(os.path.join(anionDir, "*.log"))

    neutrals = processIonLogs(neutralLogs, "neutral", densityStr)
    anions = processIonLogs(anionLogs, "anion", densityStr)
    cations = processIonLogs(cationLogs, "cation", densityStr)
    set1 = list(set(list(neutrals.keys())) & set(list(anions.keys())))
    commonPaths = list(set(set1) & set(list(cations.keys())))         
    
    substrateList = consolidatePaths(commonPaths)

    for substrate in substrateList:
        if not os.path.exists(outputDir + "/" + str(substrate)): 
            os.makedirs(outputDir + "/" + str(substrate))
        #print(substrate)
        substrateSMILES = substrateScope.loc[substrateScope['ID'] == substrate, "SMILES"].values[0]

        neutralCharges = extractFromKeys(neutrals , substrate) #atomID : [atomType , number]
        cationCharges = extractFromKeys(cations , substrate)
        anionCharges = extractFromKeys(anions, substrate)
        logPaths = sorted([log for log in neutralLogs if substrate in log])

        neutralBoltzmannWeights = getBoltzmannWeightsGauss(logPaths, 298, bltzmannStr)
        createCSV(neutralBoltzmannWeights ,outputDir + "/" + str(substrate) , "boltzmannWeights" )
        logList = neutralBoltzmannWeights["logID"]
        boltzWeights = neutralBoltzmannWeights["boltzWeights"]

        dfNeut = atomicBoltzmannConstruction(logList , boltzWeights,  neutralCharges, ["atom ID" , "atom Type" , "pop_Neut"]  )
        dfCat = atomicBoltzmannConstruction(logList , boltzWeights,  cationCharges, ["atom ID" , "atom Type" , "pop_Cat"]  )
        dfAnion = atomicBoltzmannConstruction(logList , boltzWeights,  anionCharges, ["atom ID" , "atom Type" , "pop_An"]  )
        dfMAST = pd.DataFrame(columns = ["atom ID" , "atom Type" , "pop_Neut" , "pop_Cat" , "pop_An"])
        for index, row in dfNeut.iterrows():
            atomID = row["atom ID"]
            catIdx = dfCat[dfNeut["atom ID"] == atomID].index[0]
            popCat = dfCat.at[catIdx, "pop_Cat"]
            anIdx = dfAnion[dfNeut["atom ID"] == atomID].index[0]
            popAn = dfAnion.at[anIdx, "pop_An"]
            new_row = {"atom ID": atomID, "atom Type": row["atom Type"], "pop_Neut": row["pop_Neut"] , "pop_Cat" : popCat , "pop_An" : popAn }
            dfMAST = pd.concat([dfMAST, pd.DataFrame([new_row])], ignore_index=True)
        # f- (nucleophilic attack): q(N) - q(N+1) 
        # f+ (electrophilic attack): q(N-1) - q(N)
        # f0 (radical attack): 0.5 * [q(N-1) - q(N+1)]
        f_Neg = dfMAST["pop_Neut"] - dfMAST["pop_An"]
        dfMAST["f_neg"] = f_Neg
        f_Pos = dfMAST["pop_Cat"] - dfMAST["pop_Neut"]
        dfMAST["f_pos"] = f_Pos
        f_Neut = 0.5 * (dfMAST["pop_Cat"] - dfMAST["pop_An"])
        dfMAST["f_neut"] = f_Neut
        with open(outputDir + "/" + str(substrate) + "/identification.dat", "w") as f:
            f.write(f"{substrate},{substrateSMILES}")
        createCSV(dfMAST ,outputDir + "/" + str(substrate), substrate)

if __name__ == "__main__":
    neutralDir = str(sys.argv[1])
    cationDir = str(sys.argv[2])
    anionDir = str(sys.argv[3])
    substrateCSV = str(sys.argv[4])
    outputDir = str(sys.argv[5])
    chargeStr = str(sys.argv[6])
    bltzmannStr = str(sys.argv[7])
    main(neutralDir , cationDir , anionDir , substrateCSV , outputDir , chargeStr , bltzmannStr)