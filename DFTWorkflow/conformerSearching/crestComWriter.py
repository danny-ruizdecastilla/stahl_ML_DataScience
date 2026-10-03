import os
import re
import sys
import glob
#Danny Ruiz de Castilla
#writes geom optimization com files from crest 
def boxGen(list_):
    lines = []
    for i , line_ in enumerate(list_):
        str = f"{line_} == [{i}]"
        lines.append(str)
    maxLine = max(len(line) for line in lines)
    Box = "+" + "-"*(maxLine) + "+\n"
    for line in lines:
        Box += "| " + line.ljust(maxLine) + " |\n"
    Box +=  "+" + "-"*(maxLine) + "+\n" 
    return Box
def listInputs(prompt:str):
    while True:
        partitionInput = input(prompt)
        partitionList = [part.strip() for part in partitionInput.split(",")]
        if len(partitionList) == 0 or any(part == '' for part in partitionList):
            partitionList = []
            break
        else:
            print("Input List:", partitionList)
            conf = int(input("Confirm your strings entry by Typing 1: "))
            if conf == 1:
                break
    return partitionList
def energyCutoff(energiesFile):
    energiesDict = {}
    with open(energiesFile , 'r') as file:
        for idx, line in enumerate(file):
            inputs = line.split("        ")
            energiesDict[idx] = float(inputs[-1].strip())
    if len(list(energiesDict.keys())) >= 80:
        energyCutoff = 3.6
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
def getSolventInputs(linkStr):
    solventHash = {}
    solventNames = ["stoichiometry" , "solventname" , "eps" , "epsinf" , "hbondacidity" , "hbondbasicity" , "SurfaceTensionAtInterface" , "CarbonAromaticity", "ElectronegativeHalogenicity"]
    for name in solventNames:
        val = input(f"Please enter the {name} for {linkStr}: ")
        solventHash[name] = val
    return solventHash
def xyzExtractor(coordsFile , pathNameMAST ,numComs, numAtoms , basisGen = False ):
    if numComs == 0:
        numComs +=1
    comCount = 0
    confHash = {}
    termMiddle = False
    atomSet = set()
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
                atomSet.add(atom)
                coords = [num.strip() for num in lineOptions if is_float(num.strip())]
                coordHash[str(idx) + "," + atom] = coords
                #print("atom List:" , len(coordHash.keys()))
            elif is_float(line.strip()):
                #Energy level
                energyLevel = float(line.strip())
        if not termMiddle:
            confHash[pathNameMAST + "_conf_" + str(comCount)] = {"EnergyLevel" : energyLevel , "coordinates" : coordHash}
    confHash = conformerDownsize(confHash , 100)
    if basisGen:
        groupHash = {}
        atomList = list(atomSet)
        atomBox = boxGen(atomList)
        while True:
            #Goal is to get user to group the atoms into corresponding basis sets, iteratively generate a new box showing decreasing number of atoms as they are grouped into basis sets
            prompt = f"Please enter the indices of the atoms you want to group into a basis set\n{atomBox}\n"
            groupedIdx = listInputs(prompt)
            #drop the grouped atoms from the atomList and atomBox
            groupedAtoms = [atomList[int(idx)] for idx in groupedIdx]
            basisName = input(f"Please enter the basis set name for the following atoms: {groupedAtoms}\n").strip()
            #Check for ECP
            while True:
                ecp = input(f"Does this basis set require an ECP? [0] for No and [1] for Yes\n").strip()
                if ecp not in ("0", "1"):
                    print("Input must be '0' or '1'")
                else:
                    break
            if ecp == "1":
                ecpName = input(f"Please enter the ECP name for this basis set: {groupedAtoms}\n").strip()
                groupHash[basisName + " " + ecpName] = groupedAtoms
            else:
                groupHash[basisName] = groupedAtoms
            for atom in groupedAtoms:
                atomList.remove(atom)
            if len(atomList) == 0:
                break
            atomBox = boxGen(atomList)
        return confHash , groupHash
    return confHash , None
def conformerDownsize(xyzHash , popThresh):
    population = list(xyzHash.keys())
    if len(population) >= popThresh:
        #We must downsize by reducing the number of keys
        energyList = []
        for key , val in xyzHash.items():
            try:
                energy = val["EnergyLevel"]
                energyList.append(energy)
            except:
                print(key , "no energy")
        energyGroups = groupThresh(energyList ,0.00001)
        allowedEnergies = []
        for group in energyGroups:
            energy = min(group)
            allowedEnergies.append(float(energy))
        newHash = {key: conf for key, conf in xyzHash.items() if float(conf["EnergyLevel"]) in allowedEnergies}
        return newHash
    else:
        return xyzHash
def geomComWriter(saveStr , coords , output ,  **kwargs):
    comFile = str(output) + "/" + str(saveStr) + ".com"
    solventHash = kwargs["solventHash"]
    metalHash = kwargs.get("metalBasis", None)
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
        if metalHash is not None and "Gen" in InputGeomLine:
            metalKeys = [metal for metal in metalHash.keys()]
            for metalKey in metalKeys: #metalKey looks like this "LanL2DZ LanL2DZ"
                valence = metalKey.split(" ")[-1]
                metalList = metalHash[metalKey]
                metalStr = ""
                for metal in metalList:
                    metalStr += f"{metal} "
                metalStr += " 0"
                f.write(f"{metalStr}\n{valence}\n")
                f.write(f"****\n")
            f.write(f"\n")
            if "ECP" in InputGeomLine:
                for metalKey in metalKeys:
                    if metalKey.split(" ") > 1:
                        core = metalKey.split(" ")[0]
                        metalList = metalHash[metalKey]
                        metalStr = ""
                        for metal in metalList:
                            metalStr += f"{metal} "
                        metalStr += " 0"
                        f.write(f"{metalStr}\n{core}\n")
                    else:
                        continue
                f.write(f"\n")
        if solventHash is not None:
            for key , val in solventHash.items():
                f.write(f"{key}={val}\n")
            f.write("\n")
def binaryInput(inputStr):
    binaryVal = input(f"{inputStr}").strip()
    if binaryVal not in ("0", "1"):
        raise ValueError("Input must be '0' or '1'")
    return bool(int(binaryVal))
def addLink(comFile , linkStr , linkName , charge , spin , **kwargs):
    solventHash = kwargs["solventHash"]
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
        if solventHash is not None:
            for key , val in solventHash.items():
                f.write(f"{key}={val}\n")
            f.write("\n")
def main(masterDir , outputDir):
    pathAvailables = glob.glob(masterDir + "/*/crest.energies")
    print(pathAvailables)
    geomOpt = str(input("Please enter the geometry optimization line: "))
    if re.search(r"Solvent=Generic,read",geomOpt, re.IGNORECASE):
        solventHash_ = getSolventInputs(" Energy Minimization")
    else:
        solventHash_ = None
    if "Gen" in geomOpt:
        basisGen = True
    netCharge = int(input("Please enter the net charge for these jobs: "))
    spin = int(input("Please enter the spin for these jobs: (2s+1): "))
    for path in pathAvailables:
        print("257")
        pathNameMAST = str(path.split("/")[-2].strip())
        cutoffKey = energyCutoff(path)
        pathDir = path.split("crest.energies")[0]
        pathXYZ = pathDir + "/" + str(pathNameMAST) + ".xyz"
        print(pathXYZ)
        coordsFile = pathDir +  "/crest_conformers.xyz"
        if os.path.exists(pathXYZ):
            with open(pathXYZ , 'r') as file:
                numAtoms =  int(file.readline().strip())
        else:
            print("exiting")
            sys.exit()
        print("266")
        xyzHash  , areMetals = xyzExtractor(coordsFile , pathNameMAST , cutoffKey, numAtoms , basisGen )
        if areMetals is not None:
            for key , vals in xyzHash.items():
                print("270")
                coords = vals["coordinates"]
                geomComWriter(key, coords , outputDir , solventHash = solventHash_, nprocs = 16 , mem = 25 , 
                            chk = str(key) + ".chk" ,  geomLine = geomOpt , netCharge = netCharge , spin = spin , metalBasis = areMetals)
        else:
            for key , vals in xyzHash.items():
                coords = vals["coordinates"]
                geomComWriter(key, coords , outputDir , solventHash = solventHash_, nprocs = 16 , mem = 25 , 
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
            if re.search(r"Solvent=Generic,read",linkStr, re.IGNORECASE):
                solventHash_ = getSolventInputs(" Energy Minimization")
            else:
                solventHash_ = None
            linkName = input(f"Enter the title for the --link--: ")
            netCharge = int(input("Please enter the net charge for these jobs: "))
            spin = int(input("Please enter the spin for these jobs: (2s+1)"))
            comFiles = glob.glob(outputDir + "/*.com")
            for com in comFiles:
                addLink(com , linkStr , linkName , netCharge , spin , solventHash = solventHash_)
        else:
            break