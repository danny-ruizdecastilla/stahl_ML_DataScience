import os
import sys
import glob
import shutil
import re
from pathlib import Path
#Danny Ruiz de Castilla
#writes geom optimization com files from crest 
def energyCutoff(energiesFile):
    energiesDict = {}
    with open(energiesFile , 'r') as file:
        for idx, line in enumerate(file):
            inputs = line.split("        ")
            energiesDict[idx] = float(inputs[-1].strip())
    if len(list(energiesDict.keys())) >= 50:
        energyCutoff = float(input(f"Highest energy deviation from L.E.C. is {list(energiesDict.values()[-1])}"))
        finalDict = {key: val for key, val in energiesDict.items() if val <= energyCutoff}
        cutoffKey = list(finalDict.keys())[-1]
    else:
        cutoffKey = len(list(energiesDict.keys())) -1
    return cutoffKey
def xyzExtractor(coordsFile , pathNameMAST , numComs , outputDir , numAtoms ):
    comCount = 0
    confHash = {}
    with open(coordsFile , 'r') as file:
        for idx , line in enumerate(file):
            if line.strip() == str(numAtoms):
                if coordHash:
                    confHash[pathNameMAST + "conf_" + comCount] = coordHash
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

def main(masterDir , outputDir):
    pathAvailables = glob.glob(masterDir + "/*/crest.energies")
    for path in pathAvailables:
        pathNameMAST = str(path.split("/")[-2].strip())
        cutoffKey = energyCutoff(path)
        pathDir = path.split("crest.energies")[0]
        pathXYZ = pathDir + "/" + str(pathNameMAST) + ".xyz"
        coordsFile = pathDir = "/crest_conformers.xyz"
        if os.path.exists(pathXYZ):
            with open(pathXYZ , 'r') as file:
                numAtoms =  int(file.readline().strip())
        else:
            sys.close()
        
        xyzHash = xyzExtractor(coordsFile , pathNameMAST , cutoffKey  , outputDir , numAtoms )
        

if __name__ == "__main__":
    masterDir = str(sys.argv[1])    
    outputDir = str(sys.argv[2])
    main(masterDir , outputDir)