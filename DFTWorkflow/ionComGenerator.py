import os 
import sys
import glob
import shutil
import re
#Danny Ruiz de Castilla
#Prepares anion and cation .com files
def copyChks(chkFile , outputDir):
    if not chkFile.endswith(".chk"):
        raise ValueError("input is not a .chk File")
    if not os.path.isfile(chkFile):
        raise FileNotFoundError(f"Source File does not exist: {chkFile}")
    if not os.path.isdir(outputDir):
        os.makedirs(outputDir)
    outputFile = os.path.join(outputDir , os.path.basename(chkFile))
    with open (chkFile , 'rb') as srcFile:
        with open(outputFile , 'wb') as destFile:
            shutil.copyfileobj(srcFile,destFile)
    return str(chkFile.split("/")[-1])
def theoryLevel(option):
    if option == 1:
        functional = "M062X"
        basisSet = "def2tzvp"
    if option == 2:
        functional = "B3LYP"
        basisSet = "6-31G(d,p)"
    return functional + "/" + basisSet
def comWriter(comFile:str, **kwargs ):
    if kwargs.get("coordinates") is None:
        fromChk = True
    else:
        fromChk = False
    try:
        nprocs = int(kwargs["nprocs"])
        mem = int(kwargs["mem"])
        chk = str(kwargs["chk"])
        theory = str(kwargs["theory"])
        netCharge = int(kwargs["electron"]) 
        spin = int(kwargs["spin"])
    except KeyError as e:
        raise ValueError(f"Missing required parameter: {e}")
    if kwargs.get("geom") is None:
        geom =  "checkpoint"
    else:
        geom = str(kwargs["geom"])
    if kwargs.get("guess") is None:
        guess = "read"
    else:
        guess = str(kwargs["guess"])
    if kwargs.get("symmetry") is None:
        symmetry = "loose"
    else:
        symmetry = str(kwargs["symmetry"])
    if kwargs.get("prop") is None:
        prop = "efg"
    else:
        prop = str(kwargs["prop"])

    with open(comFile, "w") as f:
        print("writing")
        f.write(f"%nprocs={nprocs}\n")
        f.write(f"%mem={mem}GB\n")
        f.write(f"%chk={chk}\n")
        fileStr = chk.split("/")[-1].split(".")[0]
        if fromChk:
            f.write(f"#{theory} geom={geom} pop=nbo7 guess={guess} symmetry={symmetry} "
                    "int=(grid=ultrafine) SP\n\n")
            f.write("Commentline\n\n")
            f.write(f"{netCharge} {spin}\n\n")

        else:
            f.write(f"#{theory} symmetry={symmetry} pop=nbo7 "
                    "scf=qc int=(grid=ultrafine) SP\n\n") 
            f.write(f"{fileStr}.xyz\n\n")
            f.write(f"{netCharge} {spin}\n")
            coordDict = kwargs["coordinates"]
            for atom in coordDict.keys():
                row = coordDict[atom]
                f.write(f"{','.join(str(i) for i in row)}\n")
            f.write("\n\n")
        
    return fileStr
def locateinLog(logFile, textStr, returnType: str):
    matchingLines = []
    
    with open(logFile, "r") as f:
        for idx, line in enumerate(f):
            if textStr in line:
                matchingLines.append(idx)
    if returnType == "earliest":
        return matchingLines[0]

    if returnType == "latest":
        return matchingLines[-1]
def getAtomCoords(logFile , xyzStr , commaSplit:int , ):
    #Extracts atom coordinates into a dict from a log file 
    atomCoords = {}
    lowerIdx = locateinLog(logFile , xyzStr, "latest" )
    upperIdx = locateinLog(logFile, "The archive entry for this job was punched." , "latest")
    masterStr = ""
    with open(logFile , "r") as f:
        for idx, line in enumerate(f):
            if idx >= lowerIdx and idx < upperIdx:
                cleaned = re.sub(r'\s+', '' , line)
                masterStr += cleaned

    masterList = masterStr.split("\\")
    for i ,  phrase in enumerate(masterList):
        atomStr = phrase.split(",")
        #print(atomStr)
        if len(atomStr) == commaSplit:
            atomCoords[i] = atomStr[:commaSplit]
    return atomCoords
    
def main(logDir , deltaE , comDir , chkDir):

    if deltaE < 0:
        ion = "anion"
    if deltaE > 0:
        ion = "cation"
    theory = theoryLevel(1)
    logs = glob.glob(logDir + "/*.log")
    print(logs)
    ionComDir = comDir + "/" + ion + "s"
    if not os.path.exists(ionComDir): 
        os.makedirs(ionComDir)
    for log in logs:
        file = log.split("/")[-1]
        fileName = file.split(".")[0]

        chkFile = chkDir + "/"  + fileName + ".chk"
        #print(chkFile)
        fileName = str(fileName) + "_" + str(ion) + ".com"
        if not os.path.exists(chkFile):
            print("Missing .chk file for " + fileName)
            atomsDict = getAtomCoords(log , "GINC-COMPUTE" , 5)
            fileName = comWriter(ionComDir + "/" + fileName , nprocs = int(16) , mem = int(48) , theory = str(theory) , chk =chkFile ,
                             electron =  deltaE , spin = int(2) , coordinates = atomsDict )          
        else:
            chkFile = copyChks(chkFile , ionComDir)

            fileName = comWriter(ionComDir + "/" + fileName , nprocs = int(16) , mem = int(48) , theory = str(theory) , chk =chkFile ,
                              geom = "checkpoint" , electron =  deltaE , spin = int(2)  )


if __name__ == "__main__":
    logDir = str(sys.argv[1])
    deltaE = int(sys.argv[2])
    comDir = str(sys.argv[3])
    chkDir = str(sys.argv[4])
    main(logDir , deltaE , comDir , chkDir)
