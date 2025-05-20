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
def fukuiFunction(substrateDF):
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
        logList.append(log)
        if os.path.exists(cationDir + "/" + logName + "_cation.log"):
            logList.append(cationDir + "/" + logName + "_cation.log")
        else:
            print("Missing the cation .log File for " + logName)
            continue
        if os.path.exists(anionDir + "/" + logName + "_anion.log"):
            logList.append(anionDir + "/" + logName + "_anion.log")     
         
        else:
            print("Missing the anion .log File for " + logName)
            continue
        if len(logList) == 3:
            substrateDicts[logName] = logList
    substrateDictMAST = {}
    for substrate, dirs in substrateDicts.items():
        substrateDF = pd.DataFrame()
        for dir in dirs:
            name = dir.split("/")[-1]
            #print(name)
            if "ion" in name:
                name1 = name.split("_")[-1].split(".")[0]
            else:
                name1 = "neutral"
            occupanciesDict = extractOccupancies(dir,"Core      Valence    Rydberg" , "* Total *" , "earliest")
            if not "atoms" in substrateDF.columns:
                substrateDF["atoms"] = occupanciesDict.keys()
            substrateDF[name1 + "_occupancies"] = list(occupanciesDict.values())

        
        substrateDictMAST[substrate] = substrateDF



        

if __name__ == "__main__":
    neutralDir = str(sys.argv[1])
    cationDir = str(sys.argv[2])
    anionDir = str(sys.argv[3])
    main(neutralDir , cationDir , anionDir)

