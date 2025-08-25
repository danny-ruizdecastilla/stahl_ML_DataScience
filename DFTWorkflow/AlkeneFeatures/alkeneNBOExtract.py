import os 
import sys
import glob
import pandas as pd
import numpy as np
from pathlib import Path
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__) , "../.."))
sys.path.append(parentDir)
from breadthFirstSearch.radialBasedCorrelation import getCC
from DFTWorkflow.fukuiGenerator.fukuiExtractorV2 import extractChargesByDensity
from DFTWorkflow.fukuiGenerator.fukuiExtractorV1 import getBoltzmannWeightsGauss

def selectColumns(df , promptStr):
    print("\nAvailable columns:\n")
    for i, col in enumerate(df.columns):
        print(f"{i}: {col}")
    
    while True:
        choice = input(f"\n{promptStr}").strip()
        if choice.isdigit():
            idx = int(choice)
            if 0 <= idx < len(df.columns):
                selected_col = df.columns[idx]
                print(f"\n✅ You selected column: {selected_col}")
                return selected_col
            else:
                print("❌ Invalid index. Try again.")
        elif choice in df.columns:
            print(f"\n✅ You selected column: {choice}")
            return choice
        
        else:
            print("❌ Invalid input. Please try again.")
def extractPiBond(logPath , C1 , C2):
    matchStr = f"BD(2)C{C1}-C{C2}"
    with open(logPath, "r") as f:
        for line in f:
            newLine = line.replace(" ", "")
            if matchStr in newLine:
                matchLine = line.strip()
                occupancy = float(matchLine.split("BD")[0].split("(")[1].split(")")[0].strip())
                return occupancy

def alkeneNBOExtractor(logList , C1 , C2 , energyStr , logNameMAST , smiles):
    weightsDF = getBoltzmannWeightsGauss(logList , 298 , energyStr)
    Cmin_NBO = []
    Cmax_NBO = []
    piBond = []
    for idx , row in weightsDF.iterrows():
        logFile = row["logID"]
        logPath = [log for log in logList if logFile in log][0]
        chargeHash = extractChargesByDensity(logPath , "NBO7")
        alkeneHash = {key: chargeHash[key] for key in chargeHash if key == C1 or key == C2}
        c1NBO = float(alkeneHash[C1][1])
        c2NBO = float(alkeneHash[C2][1])
        if c2NBO >= c1NBO:
            Cmax_NBO.append(c2NBO)
            Cmin_NBO.append(c1NBO)
        else:
            Cmax_NBO.append(c1NBO)
            Cmin_NBO.append(c2NBO)

        piOccupancy = extractPiBond(logPath , C1 , C2)
        piBond.append(piOccupancy)
    weightsDF["NBO_Cmin"] = Cmin_NBO
    weightsDF["NBO_Cmax"] = Cmax_NBO
    weightsDF["piBond"] = piBond
    finalHash = {}
    finalHash["ID"] = logNameMAST
    finalHash["SMILES"] = smiles
    finalHash["NBO_mxAlk"] =((weightsDF["NBO_Cmax"] * weightsDF["boltzWeights"]).sum() / weightsDF["boltzWeights"].sum())
    finalHash["NBO_mnAlk"] =((weightsDF["NBO_Cmin"] * weightsDF["boltzWeights"]).sum() / weightsDF["boltzWeights"].sum())
    finalHash["piBond"] =((weightsDF["piBond"] * weightsDF["boltzWeights"]).sum() / weightsDF["boltzWeights"].sum())
    return finalHash
def main(nboDir , substrateCSV , outputDir):
    dfMAST = pd.DataFrame()
    substrateScopeDF = pd.read_csv(substrateCSV , encoding='utf-8')
    smilesCol = selectColumns(substrateScopeDF , "Select the colunn name for SMILES strings: ")
    idCol = selectColumns(substrateScopeDF , "Select the colunn name for molecule ID: ")
    logFiles = glob.glob(nboDir + "/*.log")
    log1 = logFiles[0]
    fileSplit = input(f"{log1} Enter the string iteral that seperates the common name with the conf. type : ")
    
    for idx , row in substrateScopeDF.iterrows():
        smiles = row[smilesCol]
        CC , molec = getCC(smiles)
        C1 = CC[0] + 1
        C2 = CC[1] + 1
        C1 = min([C1, C2])
        C2 = max([C1 , C2])

        molecStr = row[idCol] 
        fileBase = str(molecStr) + str(fileSplit)
        #print(fileBase)
        eligibleLogs = [log for log in logFiles if fileBase in log]
        if len(eligibleLogs) != 0:
            #print(eligibleLogs)
            nboHash = alkeneNBOExtractor(eligibleLogs , C1 , C2 , 'electronic' , molecStr, smiles )
            dfMAST = pd.concat([dfMAST , pd.DataFrame([nboHash])] , ignore_index=True)
    dfMAST["NBO_mean"]  = (dfMAST["NBO_Cmax"] + dfMAST["NBO_Cmin"]) / 2
    dfMAST["NBO_delta"]  = (dfMAST["NBO_Cmax"] - dfMAST["NBO_Cmin"]) 
    dfMAST.to_csv(outputDir  +"/AlkeneNBO.csv")

if __name__ == "__main__":
    nboDir = str(sys.argv[1])
    substrateCSV = str(sys.argv[2])
    outputDir = str(sys.argv[3])
    if not os.path.exists(outputDir ): 
        os.makedirs(outputDir) 
    main(nboDir , substrateCSV , outputDir)
