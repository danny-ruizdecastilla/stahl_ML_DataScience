#clean up log files before extracting information
#Danny Ruiz de Castilla 09.19.2025
import sys
import glob
import numpy as np
import os
from pathlib import Path

def basicTerm(logFile, errorPhrase , termPhrase):
    termCount = 0
    with open(logFile, 'r') as f:
        for line in f:
            if errorPhrase in line:
                return True
            elif termPhrase in line:
                termCount +=1
    if termCount == 0:
        return True
    else:
        return False
def emptyLog(logFile , minLines):
    lineCount = 0
    with open(logFile , 'r') as f:
        for line in f:
            lineCount +=1
    if lineCount <= minLines:
        return True
    else:
        return False
def resub(logFiles , comFiles ):
    resubDict = {}
    com1 = comFiles[0]
    directory = Path(com1).parent
    resub = directory / "resub"
    resub.mkdir(exist_ok=True) 

    for log in logFiles:
        empty = emptyLog(log , 25)
        if empty:
            file1 = Path(log)
            fileName = str(file1.name)
            mainName = fileName.split(".")[0]
            comName = mainName + ".com"
            comFile = [file for file in comFiles if comName in file][0]
            resubDict[mainName] = [str(log) , str(comFile) ]
    
    for name , paths in resubDict.items():
        logFile = paths[0]
        comFile = paths[1]
        logPath = Path(logFile)

        if logPath.exists():
            logPath.unlink()
            print(f"{logPath} deleted.")
        comPath = Path(comFile)
        newCom = resub / comPath.name
        comPath.rename(newCom)

if __name__ == "__main__":
    logDir = sys.argv[1]
    functionStr = sys.argv[2]

    p = Path(logDir)
    logFiles = list(glob.glob(logDir + "/*.log"))

    if functionStr == "resub":
        while True:
            comDir = input("Enter the directory where the com files are located: ")
            comP = Path(comDir)
            if comP.exists() and comP.is_dir():
                break
            else:
                print(f"The chosen path {comDir} does not exist")

        comFiles = list(glob.glob(comDir + "/*.com"))
        resub(logFiles, comFiles)