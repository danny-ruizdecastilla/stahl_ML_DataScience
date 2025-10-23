#creates com files from any Gaussian Log outputs that exist in a directory
#Uses optimized geometries 
#Danny Ruiz de Castilla 08.20.2025
import sys
import glob
import numpy as np
import os
from pathlib import Path
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from DFTWorkflow.ionComGenerator import copyChks ,  getAtomCoords
from DFTWorkflow.conformerSearching.crestComWriter import binaryInput , addLink
def comWriter(comFile:str, fromChk , **kwargs ):
    try:
        nprocs = int(kwargs["nprocs"])
        mem = int(kwargs["mem"])
        chk = str(kwargs["chk"])
        InputGeomLine = str(kwargs['geomLine'])
        netCharge = int(kwargs["netCharge"])
        spin = int(kwargs["spin"])
    except KeyError as e:
        raise ValueError(f"Missing required parameter: {e}")
    with open(comFile, "w") as f:
        print("writing")
        f.write(f"%nprocs={nprocs}\n")
        f.write(f"%mem={mem}GB\n")
        f.write(f"%chk={chk}\n")
        #fileStr = chk.split("/")[-1].split(".")[0]
        fileStr = chk.root.split(".")[0]
        if fromChk:
            f.write(f"{InputGeomLine}\n\n")
            f.write(" Post Opt\n\n")
            f.write(f"{netCharge} {spin}\n")

        else:
            f.write(f"{InputGeomLine}\n\n")
            f.write(f"{fileStr}.xyz\n\n")
            f.write(f"{netCharge} {spin}\n")
            coordDict = kwargs["coordinates"]
            for atom in coordDict.keys():
                row = coordDict[atom]
                f.write(f"{','.join(str(i) for i in row)}\n")
            f.write("\n\n")

    return fileStr
def main(logDir  , outputDir ):
    logs = list(Path(logDir).glob("*.log"))
    if not os.path.exists(outputDir): 
        os.makedirs(outputDir)
    geomOpt = str(input("Please enter the geometry optimization line: "))
    netCharge = int(input("Please enter the net charge for these jobs: "))
    spin = int(input("Please enter the spin for these jobs: (2s+1): "))
    if netCharge < 0:
        addendum = "_anion.com"
    elif netCharge > 0:
        addendum = "_cation.com"
    else:
        addendum = ".com"
    for log in logs:
        fileName = str(log.name.split(".")[0])
        chkFile = fileName + ".chk"
        #print(chkFile)
        fileName = str(fileName) + str(addendum)
        if not os.path.exists(chkFile):
            print("Missing .chk file for " + fileName)
            atomsDict = getAtomCoords(log , "GINC-COMPUTE" , 5)
            if atomsDict == "Poison":
                print(f"Failed to find coordinates for {log}")
                continue
            outputPath = Path(outputDir) / fileName
            fileName = comWriter(outputPath, nprocs = int(16) , mem = int(48) ,  chk =Path(chkFile) ,
                             netCharge =  netCharge , spin = spin , coordinates = atomsDict ,geomLine = geomOpt , fromChk = False )
        else:
            chkFile = copyChks(chkFile , outputDir)
if __name__ == "__main__":
    logDir = str(sys.argv[1])
    comDir = str(sys.argv[2])
    chkDir = str(sys.argv[3])
    main(logDir , comDir )
    while True:
        link = binaryInput(f"Select if you want to add a --link--\n[0]   No Link\n[1]    Add Link\n")
        if link:
            linkStr = input(f"Write out the input line for your --link--\n")
            linkName = input(f"Enter the title for the --link--: ")
            netCharge = int(input("Please enter the net charge for these jobs: "))
            spin = int(input("Please enter the spin for these jobs: (2s+1): "))
            comFiles = list(Path(comDir).glob("*.com"))
            for com in comFiles:
                addLink(com , linkStr , linkName , netCharge , spin)
        else:
            break