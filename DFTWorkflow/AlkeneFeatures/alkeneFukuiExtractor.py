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
            



if __name__ == "__main__":
    neutralDir = str(sys.argv[1])
    cationDir = str(sys.argv[2])
    anionDir = str(sys.argv[3])
    substrateCSV = str(sys.argv[4])
    outputDir = str(sys.argv[5])
    chargeStr = str(sys.argv[6])
    bltzmannStr = str(sys.argv[7])
    main(neutralDir , cationDir , anionDir , substrateCSV , outputDir , chargeStr , bltzmannStr)