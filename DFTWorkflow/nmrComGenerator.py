import os
import sys
import glob
import shutil
import re
#Danny Ruiz de Castilla
#Edits com files
def editComs(comFile , comDict , comType , chkName):
    machineStr = f"%nprocs=16\n%mem=48GB\n%chk={chkName}.chk\n"
    geomTheory = comDict["geom"]
    solventType = comDict["solvent"]
    HNMR = comDict["HNMR"]
    CNMR = comDict["CNMR"]
    
    if comType == 1:
        geomTheory = f"# {geomTheory} empiricaldispersion=GD3BJ int=(grid=fine) opt freq=noraman\n"
        spEnergy = comDict["energy"]
        linkList = [f"--Link1--\n%nprocs=16\n%mem=48GB\n%chk={chkName}.chk\n# {spEnergy} int=(grid=fine) SP guess=read geom=check {solventType}\n\n SPE \n\n0 1\n" , 
                    f"\n--Link1--\n%nprocs=16\n%mem=48GB\n%chk={chkName}.chk\n# {HNMR} int=(grid=fine) nmr=giao guess=read geom=check {solventType}\n\n proton NMR\n\n0 1\n" , 
                    f"\n--Link1--\n%nprocs=16\n%mem=48GB\n%chk={chkName}.chk\n# {CNMR} int=(grid=fine) nmr=giao guess=read geom=check {solventType}\n\n c13 NMR\n\n0 1\n\n"]

    elif comType == 2:
        geomTheory = f"# {geomTheory} int=(grid=fine) opt freq=noraman {solventType}\n"
        linkList = [f"--Link1--\n%nprocs=16\n%mem=48GB\n%chk={chkName}.chk\n# {HNMR} int=(grid=fine) nmr=giao guess=read geom=check {solventType}\n\n proton NMR\n\n0 1\n" , 
                    f"\n--Link1--\n%nprocs=16\n%mem=48GB\n%chk={chkName}.chk\n# {CNMR} int=(grid=fine) nmr=giao guess=read geom=check {solventType}\n\n c13 NMR\n\n0 1\n\n"]

    modifiedLines = []
    with open(comFile, 'r') as f:
        for line in f:
            if "!" in line:
                modifiedLines.append(machineStr)
            elif "#" in line:
                modifiedLines.append(geomTheory)
            else:
                modifiedLines.append(line)
    with open(comFile, 'w') as f:
        f.writelines(modifiedLines)
        f.writelines(linkList)

def readComs(mainDir , comDict  , comType):
    allComs = glob.glob(mainDir + "/*.com")
    for com in allComs:
        chkName = com.split("/")[-1].split(".")[0]
        editComs(com , comDict , comType , chkName)

if __name__ == "__main__":
    comDir = str(sys.argv[1])
    while True:
        comType = input("Do you want method [1] for Speed + Efficiency or method [2] for High Accuracy? ")
        if comType in ("1", "2"):
            comType = int(comType)
            break
        else:
            print("Invalid input. Please enter 1 or 2.")
    if comType == 1:
        inputDict = {"geom" : "B3LYP/6-31G(d)" , "energy" : "B3LYP-D3/6-31G(d)" , "solvent": "SCRF=(PCM,Solvent=Generic)" , "HNMR": "WPO4/jul-cc-pVDZ" , "CNMR": "wB97X-D/def2-SVP",  }
    elif comType == 2:
        inputDict = {"geom" : "B3LYP-D3/6-31G(d)" , "energy" : "B3LYP-D3/6-31G(d)" , "solvent": "SCRF=(PCM,Solvent=Generic)" , "HNMR": "WPO4/6-311++G(2d,p)" , "CNMR": "wB97X-D/def2-SVP",  }
    readComs(comDir , inputDict , comType)