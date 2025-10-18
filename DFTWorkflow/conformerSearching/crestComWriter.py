import os
import sys
import glob
#Danny Ruiz de Castilla
#writes geom optimization com files from crest 
def energyCutoff(energiesFile):
    energiesDict = {}
    with open(energiesFile , 'r') as file:
        for idx, line in enumerate(file):
            inputs = line.split("        ")
            energiesDict[idx] = float(inputs[-1].strip())
    if len(list(energiesDict.keys())) >= 20:
        energyCutoff = 0.75
        finalDict = {key: val for key, val in energiesDict.items() if val <= energyCutoff}
        cutoffKey = list(finalDict.keys())[-1]
    else:
        cutoffKey = len(list(energiesDict.keys()))
    return cutoffKey
def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
def groupThresh(values, threshold):
    values = sorted(values)
    groups = [[values[0]]]
    
    for v in values[1:]:
        if abs(v - groups[-1][-1]) <= threshold:
            groups[-1].append(v)
        else:
            groups.append([v])
    
    return groups
def xyzExtractor(coordsFile , pathNameMAST ,numComs, numAtoms ):
    if numComs == 0:
        numComs +=1
    comCount = 0
    confHash = {}
    termMiddle = False
    with open(coordsFile , 'r') as file:
        coordHash = None
        for idx , line in enumerate(file):
            #print(line)
            if comCount ==numComs and coordHash is not None:
                print(coordHash)
                confHash[pathNameMAST + "_conf_" + str(comCount)] = {"EnergyLevel" : energyLevel , "coordinates" : coordHash}
                termMiddle = True
                break
            elif line.strip() == str(numAtoms):#Begin
                if coordHash is None:#intialize
                    coordHash = {}
                else: #new 3d conformer
                    confHash[pathNameMAST + "_conf_" + str(comCount)] = {"EnergyLevel" : energyLevel , "coordinates" : coordHash}
                    coordHash = {}
                    comCount +=1
            elif line.split("         ")[0].strip()[:1].isalpha():
                lineOptions = line.strip().split("    ")
                atom = str(lineOptions[0])
                coords = [num.strip() for num in lineOptions if is_float(num.strip())]
                coordHash[str(idx) + "," + atom] = coords
                #print("atom List:" , len(coordHash.keys()))
            elif is_float(line.strip()):
                #Energy level
                energyLevel = float(line.strip())
        if not termMiddle:
            confHash[pathNameMAST + "_conf_" + str(comCount)] = {"EnergyLevel" : energyLevel , "coordinates" : coordHash}
    confHash = conformerDownsize(confHash , 15)
    return confHash
def conformerDownsize(xyzHash , popThresh):
    population = list(xyzHash.keys())
    if len(population) >= popThresh:
        #We must downsize by reducing the number of keys
        print(len(population))
        energyList = []
        for key , val in xyzHash.items():
            try:
                energy = val["EnergyLevel"]
                energyList.append(energy)
            except:
                print(key , "no energy")
        energyGroups = groupThresh(energyList ,0.001)
        allowedEnergies = []
        for group in energyGroups:
            energy = min(group)
            allowedEnergies.append(float(energy))
        newHash = {key: conf for key, conf in xyzHash.items() if float(conf["EnergyLevel"]) in allowedEnergies}
        return newHash
    else:
        return xyzHash
def geomComWriter(saveStr , coords , output , **kwargs):
    comFile = str(output) + "/" + str(saveStr) + ".com"
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
        f.write(f"{InputGeomLine}\n\n")
        f.write(" Energy Minimization\n\n")
        f.write(f"{netCharge} {spin}\n")
        for atom, coordinates in coords.items():
            f.write(f"{atom.split(',')[-1]} {coordinates[0]} {coordinates[1]} {coordinates[2]}\n")
        f.write("\n")
    
def binaryInput(inputStr):
    binaryVal = input(f"{inputStr}").strip()
    if binaryVal not in ("0", "1"):
        raise ValueError("Input must be '0' or '1'")
    return bool(int(binaryVal))
def addLink(comFile , linkStr , linkName , charge , spin):
    with open(comFile, 'r') as f:
        for idx , line in enumerate(f):
            if "nprocs" in line:
                nproc = int(line.split("=")[-1].strip())
            elif "mem=" in line:
                mem = str(line.split("=")[-1].strip())
    lnkStr = comFile.split("/")[-1].split(".")[0]
    chkFile = lnkStr + ".chk"
    newName = lnkStr + linkName
    with open(comFile, 'a') as f:
        f.write(f"--Link1--\n")
        f.write(f"%nprocs={nproc}\n")
        f.write(f"%mem={mem}\n")
        if charge > 0:
            f.write(f"%oldchk={chkFile}\n")
            newChk = lnkStr + "_cation.chk"
            f.write(f"%chk={newChk}\n")
        elif charge < 0:
            f.write(f"%oldchk={chkFile}\n")
            newChk = lnkStr + "_anion.chk"
            f.write(f"%chk={newChk}\n")    
        elif charge == 0:
            f.write(f"%chk={chkFile}\n")  
        f.write(f"{linkStr}\n\n")
        f.write(f" {newName}\n\n")
        f.write(f"{charge} {spin}\n\n")
            
def main(masterDir , outputDir):
    pathAvailables = glob.glob(masterDir + "/*/crest.energies")
    geomOpt = str(input("Please enter the geometry optimization line: "))
    netCharge = int(input("Please enter the net charge for these jobs: "))
    spin = int(input("Please enter the spin for these jobs: (2s+1): "))
    for path in pathAvailables:
        pathNameMAST = str(path.split("/")[-2].strip())
        cutoffKey = energyCutoff(path)
        pathDir = path.split("crest.energies")[0]
        pathXYZ = pathDir + "/" + str(pathNameMAST) + ".xyz"
        coordsFile = pathDir +  "/crest_conformers.xyz"
        if os.path.exists(pathXYZ):
            with open(pathXYZ , 'r') as file:
                numAtoms =  int(file.readline().strip())
        else:
            sys.exit()
        
        xyzHash = xyzExtractor(coordsFile , pathNameMAST , cutoffKey, numAtoms )
        for key , vals in xyzHash.items():
            coords = vals["coordinates"]
            geomComWriter(key, coords , outputDir , nprocs = 16 , mem = 25 , 
                          chk = str(key) + ".chk" ,  geomLine = geomOpt , netCharge = netCharge , spin = spin)

if __name__ == "__main__":
    masterDir = str(sys.argv[1])    
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir): 
        os.makedirs(outputDir)
    main(masterDir , outputDir)
    while True:
        link = binaryInput(f"Select if you want to add a --link--\n[0]   No Link\n[1]    Add Link\n")
        if link:
            linkStr = input(f"Write out the input line for your --link--\n")
            linkName = input(f"Enter the title for the --link--: ")
            netCharge = int(input("Please enter the net charge for these jobs: "))
            spin = int(input("Please enter the spin for these jobs: (2s+1)"))
            comFiles = glob.glob(outputDir + "/*.com")
            for com in comFiles:
                addLink(com , linkStr , linkName , netCharge , spin)
        else:
            break
