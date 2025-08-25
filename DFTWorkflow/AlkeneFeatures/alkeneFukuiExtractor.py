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
def getAlkeneFukuiFunctions(neutralFiles , cationFiles , anionFiles , weightsDF , C1 , C2 , chargeStr ):
    f_pos_Cmin = []
    f_neg_Cmin = []
    f_neut_Cmin = []

    f_pos_Cmax = []
    f_neg_Cmax = []
    f_neut_Cmax = []

    f_pos_Scaled = []
    f_neg_Scaled = []
    f_neut_Scaled = []
    for name in list(weightsDF["logID"]):
        numStr = name.split("conf")[-1]
        confNum = getConfNum(numStr)

        neutralStr = str(name) + ".log"
        cationStr = "_" + confNum + "_cation"
        anionStr = "_" + confNum + "_anion"

        neutralFile = [file for file in neutralFiles if neutralStr in file][0]
        cationFile = [file for file in cationFiles if cationStr in file][0]
        anionFile = [file for file in anionFiles if anionStr in file][0]

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
                    f_plus_C1 = f_plus
                    f_neut_C1 = f_neut
                elif atomIdx == C2:
                    f_min_C2 = f_minus
                    f_plus_C2 = f_plus
                    f_neut_C2 = f_neut 
        f_neg_Scaled.append(np.mean([f_min_C1 , f_min_C2]) - np.mean(fukuiNeg) )/ np.max(fukuiNeg)
        f_pos_Scaled.append(np.mean([f_plus_C1 , f_plus_C2]) - np.mean(fukuiPlus) )/ np.max(fukuiPlus)
        f_neut_Scaled.append(np.mean([f_neut_C1 , f_neut_C2]) - np.mean(fukuiNeut) )/ np.max(fukuiNeut)



def main(logDir , cationDir , anionDir , substrateCSV , outputDir , chargeStr , bltzmannStr):
    dfMAST = pd.DataFrame()
    substrateScopeDF = pd.read_csv(substrateCSV , encoding='utf-8')
    smilesCol = selectColumns(substrateScopeDF , "Select the colunn name for SMILES strings: ")
    idCol = selectColumns(substrateScopeDF , "Select the colunn name for molecule ID: ")

    neutralLogs = glob.glob(os.path.join(logDir, "*.log"))
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
        #print(fileBase)
        eligibleNeutrals = [log for log in neutralLogs if fileBase in log].sort()
        eligibleCations = [log for log in cationLogs if fileBase in log].sort()
        eligibleAnions = [log for log in anionLogs if fileBase in log].sort()
        if len(eligibleNeutrals) == len(eligibleCations) == len(eligibleAnions):
            weightsDF = getBoltzmannWeightsGauss(eligibleNeutrals , 298 , bltzmannStr)
            weightsDF = getAlkeneFukuiFunctions(eligibleNeutrals , eligibleCations , eligibleAnions, weightsDF , C1 , C2 , chargeStr )


if __name__ == "__main__":
    neutralDir = str(sys.argv[1])
    cationDir = str(sys.argv[2])
    anionDir = str(sys.argv[3])
    substrateCSV = str(sys.argv[4])
    outputDir = str(sys.argv[5])
    chargeStr = str(sys.argv[6])
    bltzmannStr = str(sys.argv[7])
    main(neutralDir , cationDir , anionDir , substrateCSV , outputDir , chargeStr , bltzmannStr)