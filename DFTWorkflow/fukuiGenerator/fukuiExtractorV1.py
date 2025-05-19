import os 
import sys
import glob
import shutil
import re
import pandas as pd 
import numpy as np
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(parentDir)
from DFTWorkflow.ionComGenerator import locateinLog

#Danny Ruiz de Castilla
#Generates Fukui maps for each molecule 
def fukuiMap(substrate:dict):
    substrateDF = pd.DataFrame()
    for atom , presences in enumerate(substrate):
        #presences = [neutral , cation , anion]
        f_0 = 0.5*(float(presences[2]) - float(presences[1]) )
        f_plus = float(presences[2]) - float(presences[0])
        f_minus = float(presences[0]) - float(presences[1])

        substrateDF[atom+ "_f_0"] = f_0
        substrateDF[atom+ "_f_plus"] = f_plus
        substrateDF[atom+ "_f_minus"] = f_minus
    return substrateDF
def extractOccupancies(logFile:str, extract1:str, extract2:str , location:str):

    lowerInd = locateinLog(logFile , extract1, location)
    upperInd = locateinLog(logFile , extract2, location )
    atomOccupancies = {}
    with open(logFile , "r") as f:
        for idx, line in enumerate(f):
            if idx > lowerInd + 1  and idx < upperInd-1:
                occupancies = line.strip()
                occupancyLines = occupancies.split("    ")
                #print(occupancyLines)
                atomNums = occupancyLines[0].split("   ")
                atomInd = ''.join(atomNums[0].split())
                #print(atomInd)
                #print(occupancyLines[1:])
                nums = []
                for num in occupancyLines[1:]:
                    num = float(''.join(num.split()))
                    nums.append(num)
                
                #print(np.max(nums))
                atomOccupancies[atomInd] = np.max(nums)
    #print(atomOccupancies)
    return atomOccupancies


def main(logDir , cationDir, anionDir ):
    neutralLogs = glob.glob(logDir + "/*.log")
    substrateDicts = {}

    for log in neutralLogs:
        logList = []
        logName = log.split("/")[-1].split(".")[0]
        logList.append(logName)
        

if __name__ == "__main__":
    neutralDir = str(sys.argv[1])
    cationDir = str(sys.argv[2])
    anionDir = str(sys.argv[3])
    main(neutralDir , cationDir , anionDir)

