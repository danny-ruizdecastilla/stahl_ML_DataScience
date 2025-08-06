import os 
import sys
import glob
import shutil
import re
import pandas as pd 
import numpy as np
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../"))
sys.path.append(parentDir)
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from DFTWorkflow.featureMaping import createCSV
from DFTWorkflow.fukuiGenerator.fukuiExtractorV1 import  getBoltzmannWeightsGauss
from DFTWorkflow.fukuiGenerator.fukuiExtractorV2 import consolidatePaths , extractFromKeys , atomicBoltzmannConstruction
from DFTWorkflow.ionComGenerator import locateinLog

def extractShifts(logFile: str, extract1:str, extract2:str,location1:str ,location2, atomList ):
    lowerInd = locateinLog(logFile , extract1, location1)
    upperInd = locateinLog(logFile , extract2, location2 )
    #print(lowerInd)
    #print(upperInd)
    if upperInd == "Poison" or lowerInd == "Poison":
        raise ValueError(f"The Log file {logFile} did not terminate properly")
    atomShifts = {}
    with open(logFile , "r") as f:
        for idx, line in enumerate(f):
            if idx > int(lowerInd)  and idx < int(upperInd):
                
                if "   Isotropic =  " in line:
                    atomInd = line.split("   Isotropic =  ")[0].strip()
                    if any(atom in atomInd for atom in atomList):
                        shiftsTots = line.split("   Isotropic =  ")[-1].strip()
                        shifts = shiftsTots.split("Anisotropy")[0]
                        atomIdx = int(atomInd.split(" ")[0].strip())
                        atomName = str(atomInd.split(" ")[-1].strip())
                        atomShifts[atomIdx] = [ atomName, float(shifts.strip())]
    return atomShifts  
def main(nmrDir, substrateCSV , outputDir):
    substrateScope = pd.read_csv(substrateCSV, encoding='utf-8' )
    nmrPaths = glob.glob(nmrDir + "/*.log")
    

    atomList = listInputs(f"Enter the list of all atom types you wnat to extract shifts for Ex: [C,H,O,N] ")
    nmrShifts = {}
    for atom in atomList:
        yIntercept = input(f"Enter the y intercept for the NMR shift scaling for {atom}: ")
        slope = input(f"Enter the slope for the NMR shift scaling for {atom}: ")
        nmrShifts[atom] = [float(yIntercept) , float(slope)]

    molecularShifts = {}
    extract1 = "SCF GIAO Magnetic shielding tensor (ppm):"
    extract2 = "Eigenvalues:  "
    for path in nmrPaths:
        shifts = extractShifts(path, extract1, extract2, "earliest", "latest", list(nmrShifts.keys()))
        name = os.path.basename(path).split(".")[0]
        molecularShifts[name] = shifts
    substrateList = consolidatePaths(list(molecularShifts.keys()))

    for substrate in substrateList:
        substrateSMILES = substrateScope.loc[substrateScope['ID'] == substrate, "SMILES"].values[0]
        if not os.path.exists(outputDir + "/" + str(substrate)): 
            os.makedirs(outputDir + "/" + str(substrate))
        substratePaths = sorted([path for path in nmrPaths if substrate in path])
        boltzmannDF = getBoltzmannWeightsGauss(substratePaths, 298, "electronic")
        createCSV(boltzmannDF ,outputDir + "/" + str(substrate) , "boltzmannWeights" )
        substrateShifts = extractFromKeys(molecularShifts , substrate)
        logList = boltzmannDF["logID"]
        boltzWeights = boltzmannDF["boltzWeights"]
        dfMAST = atomicBoltzmannConstruction(logList , boltzWeights,  substrateShifts, ["atom ID" , "atom Type" , "NMR Shift"]  )
        for index, row in dfMAST.iterrows():
            atomType = str(row["atom Type"])
            if atomType in nmrShifts:
                y_intercept, slope = nmrShifts[atomType]
                scaledShift = (y_intercept - row["NMR Shift"]) / (-1 * slope)
                dfMAST.loc[index, "NMR Shift"] = scaledShift
        with open(outputDir + "/" + str(substrate) + "/identification.dat", "w") as f:
            f.write(f"{substrate},{substrateSMILES}")
        createCSV(dfMAST ,outputDir + "/" + str(substrate), substrate)
if __name__ == "__main__":
    nmrDir = str(sys.argv[1])
    substrateCSV = str(sys.argv[2])
    outputDir = str(sys.argv[3])
    main(nmrDir , substrateCSV, outputDir )