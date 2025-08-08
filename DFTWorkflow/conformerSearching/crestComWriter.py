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
    if len(list(energiesDict.keys())) >= 50:
        energyCutoff = float(input(f"Highest energy deviation from L.E.C. is {list(energiesDict.values())[-1]}, Enter an energy cutoff for conformers: "))
        finalDict = {key: val for key, val in energiesDict.items() if val <= energyCutoff}
        cutoffKey = list(finalDict.keys())[-1]
    else:
        cutoffKey = len(list(energiesDict.keys())) -1
    return cutoffKey
def xyzExtractor(coordsFile , pathNameMAST , numComs , outputDir , numAtoms ):
    comCount = 0
    confHash = {}
    with open(coordsFile , 'r') as file:
        coordHash = None
        for idx , line in enumerate(file):
            if line.strip() == str(numAtoms):#new 3d conformer
                if coordHash is not None:
                    confHash[pathNameMAST + "_conf_" + str(comCount)] = coordHash
                    coordHash = {}
                else:
                    coordHash = {}
                comCount +=1
            elif line.split("         ")[0].strip()[:1].isalpha():
                atom = str(line.split("         ")[0].strip())
                coords = str(line.split(("         ")[-1]))
                coordHash[atom] = [coords]
            if comCount == numComs:
                break
    return confHash
def geomComWriter(saveStr , coords , theory , output , **kwargs):
    comFile = str(output) + "/" + str(saveStr) + ".com"
    try:
        nprocs = int(kwargs["nprocs"])
        mem = int(kwargs["mem"])
        chk = str(kwargs["chk"])
        symmetry = str(kwargs['symmetry'])
        correction = str(kwargs['corrections'])
        netCharge = int(kwargs["netCharge"])
        spin = int(kwargs["spin"])
    except KeyError as e:
        raise ValueError(f"Missing required parameter: {e}")
    if correction == "True":
        empiricalDispersion = 'empiricaldispersion=GD3BJ'
    with open(comFile, "w") as f:
        print("writing")
        f.write(f"%nprocs={nprocs}\n")
        f.write(f"%mem={mem}GB\n")
        f.write(f"%chk={chk}\n")
        if empiricalDispersion is not None:
            f.write(f"#{theory} symmetry={symmetry} {empiricalDispersion} scf=tight int=(grid=superfine) opt=calcfc freq\n\n")
        else:
            f.write(f"#{theory} symmetry={symmetry} scf=tight int=(grid=superfine) opt=calcfc freq\n\n")
        f.write(" Energy Minimization\n\n")
        f.write(f"{netCharge} {spin}\n\n")
        for atom , coordinates in coords.items():
            f.write(f"{atom} {coordinates[0]}")
        f.write("\n\n")
    

def main(masterDir , outputDir):
    pathAvailables = glob.glob(masterDir + "/*/crest.energies")
    theory = input("Please enter the level of theory for geom. optimization: ")
    netCharge = int(input("Please enter the net charge for these jobs: "))
    spin = int(input("Please enter the spin for these jobs: "))
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
            sys.close()
        
        xyzHash = xyzExtractor(coordsFile , pathNameMAST , cutoffKey  , outputDir , numAtoms )
        for key , coords in xyzHash.items():
            geomComWriter(key, coords , theory , outputDir , nprocs = 16 , mem = 25 , 
                          chk = str(key) + ".chk" , symmetry = 'tight' , corrections = 'True' , netCharge = netCharge , spin = spin)

if __name__ == "__main__":
    masterDir = str(sys.argv[1])    
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir): 
        os.makedirs(outputDir)
    main(masterDir , outputDir)