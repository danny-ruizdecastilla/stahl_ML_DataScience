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
from DFTWorkflow.AlkeneFeatures.alkeneNBOExtract import selectColumns
# f- (nucleophilic attack): q(N) - q(N+1) 
# f+ (electrophilic attack): q(N-1) - q(N)
# f0 (radical attack): 0.5 * [q(N-1) - q(N+1)]
def getConfNum(numStr):
    num = ""
    for s in numStr:
        try: 
            int(s)
            num += s
        except:
            pass
    num = int(num)
    return num
def getAlkeneFukuiFunctions(neutralFiles , cationFiles , anionFiles , weightsDF , C1 , C2 , chargeStr , logNameMAST , smiles ):
    #print(logNameMAST)
    f_pos_Cmin = []
    f_neg_Cmin = []
    f_neut_Cmin = []

    f_pos_Cmax = []
    f_neg_Cmax = []
    f_neut_Cmax = []

    for name in list(weightsDF["logID"]):
        numStr = name.split("conf")[-1]
        confNum = getConfNum(numStr)

        neutralStr = str(name) + ".log"
        cationStr = "_" + str(confNum) + "_cation"
        anionStr = "_" + str(confNum) + "_anion"

        neutralFile = [f for f in neutralFiles if neutralStr in f.name][0]
        cationFile  = [f for f in cationFiles  if cationStr  in f.name][0]
        anionFile   = [f for f in anionFiles   if anionStr   in f.name][0]

        neutralCharges = extractChargesByDensity(neutralFile , chargeStr)
        cationCharges = extractChargesByDensity(cationFile , chargeStr)
        anionCharges = extractChargesByDensity(anionFile , chargeStr)

        fukuiPlus = []
        fukuiNeg = []
        fukuiNeut = []

        for atomIdx in list(neutralCharges.keys()):
            atomN = str(neutralCharges[atomIdx][0])
            atomC = str(cationCharges[atomIdx][0])
            atomA = str(anionCharges[atomIdx][0])

            if atomN == atomC == atomA:
                f_minus = float(neutralCharges[atomIdx][-1]) - float(anionCharges[atomIdx][-1])
                f_plus = float(cationCharges[atomIdx][-1]) - float(neutralCharges[atomIdx][-1])
                f_neut = 0.5 * (float(cationCharges[atomIdx][-1]) - float(anionCharges[atomIdx][-1]))
                fukuiPlus.append(f_plus)
                fukuiNeg.append(f_minus)
                fukuiNeut.append(f_neut)
                if atomIdx == C1:
                    f_min_C1 = f_minus
                    f_neg_Cmin.append(f_min_C1)
                    f_plus_C1 = f_plus
                    f_pos_Cmin.append(f_plus_C1)
                    f_neut_C1 = f_neut
                    f_neut_Cmin.append(f_neut_C1)
                elif atomIdx == C2:
                    f_min_C2 = f_minus
                    f_neg_Cmax.append(f_min_C2)
                    f_plus_C2 = f_plus
                    f_pos_Cmax.append(f_plus_C2)
                    f_neut_C2 = f_neut
                    f_neut_Cmax.append(f_neut_C2)
    weightsDF["f_neg_mxAlk"] = f_neg_Cmax
    weightsDF["f_pos_mxAlk"] = f_pos_Cmax
    weightsDF["f_neut_mxAlk"] = f_neut_Cmax
    weightsDF["f_neg_mnAlk"] = f_neg_Cmin
    weightsDF["f_pos_mnAlk"] = f_pos_Cmin
    weightsDF["f_neut_mnAlk"] = f_neut_Cmin

    finalHash = {}
    finalHash["ID"] = logNameMAST
    finalHash["SMILES"] = smiles
    finalHash["f_neg_mxAlk"] =((weightsDF["f_neg_mxAlk"] * weightsDF["boltzWeights"]).sum() / weightsDF["boltzWeights"].sum())
    finalHash["f_pos_mxAlk"] =((weightsDF["f_pos_mxAlk"] * weightsDF["boltzWeights"]).sum() / weightsDF["boltzWeights"].sum())
    finalHash["f_neut_mxAlk"] =((weightsDF["f_neut_mxAlk"] * weightsDF["boltzWeights"]).sum() / weightsDF["boltzWeights"].sum())

    finalHash["f_neg_mnAlk"] =((weightsDF["f_neg_mnAlk"] * weightsDF["boltzWeights"]).sum() / weightsDF["boltzWeights"].sum())
    finalHash["f_pos_mnAlk"] =((weightsDF["f_pos_mnAlk"] * weightsDF["boltzWeights"]).sum() / weightsDF["boltzWeights"].sum())
    finalHash["f_neut_mnAlk"] =((weightsDF["f_neut_mnAlk"] * weightsDF["boltzWeights"]).sum() / weightsDF["boltzWeights"].sum())
    
    finalHash["f_neg_Delta"] = finalHash["f_neg_mxAlk"] - finalHash["f_neg_mnAlk"]
    finalHash["f_pos_Delta"] = finalHash["f_pos_mxAlk"] - finalHash["f_neut_mnAlk"]
    finalHash["f_neut_Delta"] = finalHash["f_neut_mxAlk"] - finalHash["f_neut_mnAlk"]

    finalHash["f_neg_mean"]  = (finalHash["f_neg_mxAlk"]  + finalHash["f_neg_mnAlk"])  / 2
    finalHash["f_pos_mean"]  = (finalHash["f_pos_mxAlk"]  + finalHash["f_pos_mnAlk"])  / 2
    finalHash["f_neut_mean"] = (finalHash["f_neut_mxAlk"] + finalHash["f_neut_mnAlk"]) / 2
    return finalHash

def main(logDir , cationDir , anionDir , substrateCSV , outputDir , chargeStr , bltzmannStr):
    dfMAST = pd.DataFrame()
    substrateScopeDF = pd.read_csv(substrateCSV , encoding='utf-8')
    smilesCol = selectColumns(substrateScopeDF , "Select the colunn name for SMILES strings: ")
    idCol = selectColumns(substrateScopeDF , "Select the colunn name for molecule ID: ")

    neutralLogs = glob.glob(os.path.join(logDir, "*.log"))
    #print("neutral Logs" , neutralLogs)
    cationLogs = glob.glob(os.path.join(cationDir, "*.log"))
    anionLogs = glob.glob(os.path.join(anionDir, "*.log"))
    log1 = neutralLogs[0]
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
        #print("fileBase" , fileBase)
        eligibleNeutrals = sorted([log for log in neutralLogs if fileBase in os.path.basename(log)])
        eligibleCations = sorted([log for log in cationLogs if fileBase in os.path.basename(log)])
        eligibleAnions = sorted([log for log in anionLogs if fileBase in os.path.basename(log)])
        if len(eligibleNeutrals) == len(eligibleCations) == len(eligibleAnions) and len(eligibleNeutrals) != 0:
            weightsDF = getBoltzmannWeightsGauss(eligibleNeutrals , 298 , bltzmannStr)
            mainHash = getAlkeneFukuiFunctions(eligibleNeutrals , eligibleCations , eligibleAnions, weightsDF , C1 , C2 , chargeStr , molecStr , smiles )
            dfMAST = pd.concat([dfMAST , pd.DataFrame([mainHash])] , ignore_index=True)
    dfMAST.to_csv(outputDir  +"/AlkenefukuiMAST.csv")
if __name__ == "__main__":
    neutralDir = str(sys.argv[1])
    cationDir = str(sys.argv[2])
    anionDir = str(sys.argv[3])
    substrateCSV = str(sys.argv[4])
    outputDir = str(sys.argv[5])
    chargeStr = str(sys.argv[6])
    bltzmannStr = str(sys.argv[7])
    main(neutralDir , cationDir , anionDir , substrateCSV , outputDir , chargeStr , bltzmannStr)