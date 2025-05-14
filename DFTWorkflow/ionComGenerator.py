import os 
import sys
import glob
#Danny Ruiz de Castilla
#Prepares anion and cation .com files
def theoryLevel(option):
    if option == 1:
        functional = "B3LYP"
        basisSet = "6-31G(d,p)"
    
    return functional + "/" + basisSet
def comWriter(comFile:str, **kwargs ):
    if kwargs["coordinates"] is None:
        fromChk = True
    else:
        fromChk = False
    try:
        nprocs = int(kwargs["nprocs"])
        mem = int(kwargs["mem"])
        chk = str(kwargs["chk"])
        theory = str(kwargs["theory"])
        geom = str(kwargs["geom"])
        netCharge = int(kwargs["electron"]) 
        spin = int(kwargs["spin"])
    except KeyError as e:
        raise ValueError(f"Missing required parameter: {e}")
    if kwargs["guess"] is None:
        guess = "read"
    if kwargs["symmetry"] is None:
        symmetry = "loose"    
    with open(comFile, "w") as f:
        f.write(f"%nprocs={nprocs}\n")
        f.write(f"%mem={mem}GB\n")
        f.write(f"%chk={chk}\n")
        f.write(f"#{theory} geom={geom} guess={guess} symmetry={symmetry} "
                "empiricaldispersion=GD3BJ int=(grid=ultrafine) SP\n\n")
        fileStr = chk.split("/")[-1].split(".")[0]
        if fromChk:
            f.write(f"{netCharge} {spin}\n")
            f.write(f"{netCharge} {spin}")
            f.write("\n")
        else:
            f.write(f"{fileStr}.xyz\n\n")
            f.write(f"{netCharge} {spin}\n")
            coordDict = kwargs["coordinates"]
            for atom in coordDict.keys():
                row = coordDict[atom]
                f.write(f"{' '.join(str(i) for i in row)}\n")
            f.write("\n\n")
    return fileStr

def main(logDir , deltaE , comDir , chkDir):

    if deltaE < 0:
        ion = "anion"
    if deltaE > 0:
        ion = "cation"
    theory = theoryLevel(1)
    logs = glob.glob(logDir + "/*.log")
    if not os.path.exists(comDir):
        os.makedirs(comDir)
    for log in logs:
        file = log.split("/")[-1]
        fileName = file.split(".")[0]

        chkFile = chkDir + fileName + ".chk"
        if os.file.exists(chkFile):
            chkPath = chkFile.copy()
        else:
            raise ValueError(f"Missing chk File {chkPath}")

        fileName = str(fileName) + "_" + str(ion) + ".com"
        fileName = comWriter(comDir + "/" + fileName , nprocs = int(16) , mem = int(48) , theory = str(theory) , chk =chkPath ,
                              geom = "checkpoint" , electron =  deltaE , spin = int(2)  )


if __name__ == "__main__":
    logDir = str(sys.argv[1])
    deltaE = int(sys.argv[2])
    comDir = str(sys.argv[3])
    chkDir = str(sys.argv[4])
    main(logDir , deltaE , comDir , chkDir)