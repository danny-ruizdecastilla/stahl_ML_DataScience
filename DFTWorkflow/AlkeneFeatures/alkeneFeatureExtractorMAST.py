#Master Feature Extractor
import os 
import sys
import glob
import tkinter
import pandas as pd 
import numpy as np
import chemdraw
import base64
import re
from pathlib import Path
from morfeus import BuriedVolume
from rdkit import Chem
from networkx import Graph
from itertools import islice
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from DFTWorkflow.cleanLogs import basicTerm
from DFTWorkflow.ionComGenerator import locateinLog
from DFTWorkflow.AlkeneFeatures.alkeneSubstitution import eVszAlkenes
from dimensionalityReduction.reactivityFeatures import boxGen
from breadthFirstSearch.radialBasedCorrelation import getCC
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
def extractNBOOccupancies(logFile , nboStr , charge:int):
    logName = str(logFile.name.split(".")[0])
    nboIdx = locateinLog(logFile , logName + f"_{nboStr}" , 1 )
    nboIdx2 = locateinLog(logFile , logName + f"_{nboStr}" , 2 )
    print(nboIdx , nboIdx2)
    nbo7Hash = {}
    with open(logFile, 'r') as f:
        extractingDensities = False
        for i, line in enumerate(islice(f, nboIdx, nboIdx2), start=nboIdx):
            if "Atom No    Charge" in line and "Density" not in line:
                extractingDensities = True 
            if extractingDensities:
                if "* Total * " in line:
                    print(line)
                    extractingDensities = False
                else:
                    print(line)
def extractShiftsByIdx(logFile: str, extract1:str, extract2:str,location1:str ,location2, idxList  , m , b):
    lowerInd = locateinLog(logFile , extract1, location1)
    upperInd = locateinLog(logFile , extract2, location2 )
    if upperInd == "Poison" or lowerInd == "Poison":
        raise ValueError(f"The Log file {logFile} did not terminate properly")
    atomShifts = {}
    with open(logFile , "r") as f:
        for idx, line in enumerate(f):
            if idx > int(lowerInd)  and idx < int(upperInd):
                if "   Isotropic =  " in line:
                    atomInd = line.split("   Isotropic =  ")[0].strip()
                    atomIdx = int(atomInd.split(" ")[0].strip())
                    if atomIdx in idxList:
                        shiftsTots = line.split("   Isotropic =  ")[-1].strip()
                        shifts = shiftsTots.split("Anisotropy")[0]
                        atomName = str(atomInd.split(" ")[-1].strip())
                        atomShifts[atomIdx] = [ atomName,(b - float(shifts.strip())) / m]
    return atomShifts  
def getAtomCoordsRobust(logFile , xyzStr , commaSplit:int , locationIdx  ):
    #Extracts atom coordinates into a dict from a log file
    atomCoords = {}
    lowerIdx = locateinLog(logFile , xyzStr, locationIdx )
    upperIdx = locateinLog(logFile, "The archive entry for this job was punched." , locationIdx)
    if not "Poison" in [lowerIdx , upperIdx]:
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
    else:
        return "Poison"
def extractEnergies(logFile , linkStr, energyStr):
    patterns = {
        'electronic':r'SCF Done', 
        'gibbs': r'Sum of electronic and thermal Free Energies',
        'enthalpy': r'Sum of electronic and thermal Enthalpies',
        'zpe': r'Sum of electronic and zero-point Energies'
    }
    try:
        energyType = patterns[energyStr]
    except:
        raise ValueError(f"Unknown energy type: {energyStr}. "
                f"Available types: {list(patterns.keys())}")
    logName = str(logFile.name.split(".")[0])
    firstIdx = locateinLog(logFile , logName + f"_{linkStr}" , "earliest" )
    secondIdx = locateinLog(logFile , logName + f"_{linkStr}" , "latest" )
    with open(logFile, 'r') as f:
        for i, line in enumerate(f):
            if firstIdx <= i < secondIdx and energyType in line:
                energyLevel = float(line.strip().split("=")[-1].split("A.U.")[0].strip())
                return energyLevel
        return "Poison"
def getBoltzmannWeightsGauss(logDirs, temperature, energyStr ):
    R = 1.987204e-3
    HARTREE_TO_KCAL = 627.5094740631
    results = []
    
    for log in logDirs:
        if not log.exists():
            print(f"Warning: File {str(log)} not found. Skipping.")
            continue        
        try:
            energy = extractEnergies(log, "energy", energyStr)
            if energy is not "Poison":
                results.append({'logID': Path(log.name).stem,  'E_Ha': energy})
            else:
                print(f"Warning: Could not extract {energyStr} energy from {log}")     
        except Exception as e:
            print(f"Error processing {str(log)}: {str(e)}")
            continue
    if len(results) == 0:
        raise ValueError(f"No energies could be extracted from any log files for the log types {logDirs[-1]}")
    df = pd.DataFrame(results)
    df['E_kCal'] = df['E_Ha'] * HARTREE_TO_KCAL
    minE = df['E_kCal'].min()
    df['rE_kCal'] = df['E_kCal'] - minE
    df['boltzmannFacts'] = np.exp(-df['rE_kCal'] / (R * temperature))
    normFacts = df['boltzmannFacts'].sum()
    df['boltzWeights'] = df['boltzmannFacts'] / normFacts

    df = df.sort_values('E_kCal').reset_index(drop=True)
    
    return df[['logID', 'E_Ha', 'rE_kCal', 'boltzWeights']]
def hCount(molec, wildCards):
    hydrogens = 0
    for wildcard in wildCards:
        atom = molec.GetAtomWithIdx(wildcard)
        atom1 = str(atom.GetSymbol())
        if atom1 == "H":
            hydrogens += 1
    return hydrogens

def getAlkenes(substratesHash , smilesHash , featureHash ):
    featuresMASTDF = pd.DataFrame()
    for id , smilesStr in smilesHash.items():
        hashList = []
        cc , molec = getCC(smilesStr)
        conformerFiles = substratesHash[id]
        boltzmannDF = getBoltzmannWeightsGauss(conformerFiles, 298, "electronic")
        idHash = {"SMILES" : smilesStr , "ID" : id}  
        hashList.append(idHash)   
        featureList = list(featureHash.keys())  
        if "C13_shift" in featureList:
            C13_C1 = []
            C13_C2 = []
            
            for name in list(boltzmannDF["logID"]):
                fileStr = f"{name}.log"
                conformer  = next((f for f in conformerFiles if fileStr in f.name), None)
                alkeneHash = extractShiftsByIdx(conformer, "SCF GIAO Magnetic shielding tensor (ppm):", "Eigenvalues:  ", "earliest", "latest", [c + 1 for c in cc] , 1.0081 , 195.6683)
                #print(alkeneHash)
                if len(C13_C1) == 0:
                    #C1 and C2 not yet defined 
                    atom0 = alkeneHash[int(cc[0]+ 1)][1]
                    atom1 = alkeneHash[int(cc[1]+ 1)][1]
                    if atom0 >= atom1:#C1 is smallest shift, C2 is largest shift
                        C1 = cc[1] + 1
                        C2 = cc[0] + 1
                    else:
                        C1 = cc[0] + 1
                        C2 = cc[1] + 1
                C13_C1.append(alkeneHash[C1][1])
                C13_C2.append(alkeneHash[C2][1])
            weights = boltzmannDF["boltzWeights"]
            C13_Cmx = ( (C13_C2 * weights).sum()    / weights.sum()  )
            C13_Cmn = ( (C13_C1 * weights).sum()    / weights.sum()  )
            C13_delta = C13_Cmx - C13_Cmn
            C13_Mean = np.mean([C13_Cmx,C13_Cmn])
            c13Hash = {"C13_Cmx" : C13_Cmx , "C13_Cmn" : C13_Cmn , "C13_delta" : C13_delta , "C13_Mean" : C13_Mean}
            hashList.append(c13Hash)
        else:
            C1 = cc[1] + 1
            C2 = cc[0] + 1 

        if "NBO7" in featureList:
            nboNeutralStr = featureList["NBO7"][0]

            neutralHash = alkeneNBOExtractor(conformerFiles , min([C1,C2]) , max([C1,C2]) , "electronic", id , smilesStr)
            hashList.append(neutralHash)
        if "fukuiParameter" in featureList:
            nboNeutralStr = featureList["NBO7"][0]
            nboCationStr = featureList["fukuiParameter"][0]
            nboAnionstr = featureList["fukuiParameter"][1]
            neutralHash = alkeneNBOExtractor(conformerFiles , min([C1,C2]) , max([C1,C2]) , "electronic", id , smilesStr)
            cationHash = alkeneNBOExtractor(conformerFiles , min([C1,C2]) , max([C1,C2]) , "electronic", id , smilesStr)
            anionHash = alkeneNBOExtractor(conformerFiles , min([C1,C2]) , max([C1,C2]) , "electronic", id , smilesStr)

        if "%Vbur" in featureList:
            Vbur_Cmn_2 = []
            Vbur_Cmn_3 = []
            Vbur_Cmx_2 = []
            Vbur_Cmx_3 = []
            
            for name in list(boltzmannDF["logID"]):
                fileStr = f"{name}.log"
                conformer  = next((f for f in conformerFiles if fileStr in f.name), None)

                coordHash = getAtomCoordsRobust(str(conformer) , "GINC-COMPUTE" , 5  , 2)
                elements = []
                coordinates = []
                for _, coords in coordHash.items():
                    elements.append(str(coords[0]))
                    coordinates.append(np.array(coords[2:5]))
                coordinates = np.array(coordinates, dtype=float)
                vburC1_2 = BuriedVolume(elements, coordinates, int(C1-1), include_hs=True, radius=2.0).fraction_buried_volume
                vburC1_3 = BuriedVolume(elements, coordinates, int(C1-1), include_hs=True, radius=3.5).fraction_buried_volume
                vburC2_2 = BuriedVolume(elements, coordinates, int(C2-1), include_hs=True, radius=2.0).fraction_buried_volume
                vburC2_3 = BuriedVolume(elements, coordinates, int(C2-1), include_hs=True, radius=3.5).fraction_buried_volume

                Vbur_Cmn_2.append(vburC1_2)
                Vbur_Cmx_2.append(vburC2_2)
                Vbur_Cmn_3.append(vburC1_3)
                Vbur_Cmx_3.append(vburC2_3)

            weights = boltzmannDF["boltzWeights"]
            Vbur_Cmx_2Ang = ( (Vbur_Cmx_2 * weights).sum()    / weights.sum()  )
            Vbur_Cmn_2Ang= ( (Vbur_Cmn_2 * weights).sum()    / weights.sum()  )    
            Vbur_Cmx_3Ang = ( (Vbur_Cmx_3 * weights).sum()    / weights.sum()  )
            Vbur_Cmn_3Ang= ( (Vbur_Cmn_3 * weights).sum()    / weights.sum()  )               

            Vbur_2Ang_delta = Vbur_Cmx_2Ang - Vbur_Cmn_2Ang
            Vbur_3Ang_delta = Vbur_Cmx_3Ang - Vbur_Cmn_3Ang
            Vbur_2Ang_mean =  np.mean([Vbur_Cmx_2Ang, Vbur_Cmn_2Ang])
            Vbur_3Ang_mean = np.mean([Vbur_Cmx_3Ang, Vbur_Cmn_3Ang])
            BurVolHash = {"2Ang_Cmx" :Vbur_Cmx_2Ang , "2Ang_Cmn" :Vbur_Cmn_2Ang ,
                            "3Ang_Cmx" :Vbur_Cmx_3Ang , "3Ang_Cmn" :Vbur_Cmn_3Ang ,
                            "2Ang_delta" :Vbur_2Ang_delta , "3Ang_delta" :Vbur_3Ang_delta ,
                            "3Ang_mean" :Vbur_3Ang_mean , "2Ang_mean" :Vbur_2Ang_mean }
            hashList.append(BurVolHash)
        if "EvsZ" in featureList:
            conformer = conformerFiles[0]
            coordHash = getAtomCoordsRobust(str(conformer) , "GINC-COMPUTE" , 5 , 2 )
            g = Graph()
            #print("CC" , CC)
            for bond in molec.GetBonds():
                start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                g.add_edge(start, end)
            CminNeighbors = list(g.neighbors(int(C1-1)))
            CminNeighbors.remove(int(C2-1))
            CminHcount = hCount(molec, CminNeighbors)
            CmaxNeighbors = list(g.neighbors(int(C2-1)))
            CmaxNeighbors.remove(int(C1-1))
            CmaxHcount = hCount(molec, CmaxNeighbors)

            hydrogens = CmaxHcount + CminHcount
            if hydrogens <=2: #Easy check for E vs Z 
                alkeneType = eVszAlkenes(molec , g , {int(C1-1): CminNeighbors} , {int(C2-1) : CmaxNeighbors} ,  coordHash)# 0 if neither E or Z, -1 if Z, 1 if E
            else:
                alkeneType = 0
            alkeneHash = {"HCount" : hydrogens , "EvsZ" : alkeneType}
            hashList.append(alkeneHash)
        if "Dist." in featureList:
            distList = []
            for name in list(boltzmannDF["logID"]):
                fileStr = f"{name}.log"
                conformer  = next((f for f in conformerFiles if fileStr in f.name), None)

                coordHash = getAtomCoordsRobust(str(conformer) , "GINC-COMPUTE" , 5 , 2 )
                idx = 0
                for _ , coords in coordHash.items():
                    idx += 1
                    if idx == C1:
                        c1_coords = coords[2:5]
                    elif idx == C2:
                        c2_coords = coords[2:5]
                dist1 = np.linalg.norm(np.array(c1_coords , dtype=float) - np.array(c2_coords , dtype=float))
                distList.append(dist1)
            distMAST = ( (distList * weights).sum()    / weights.sum()  )  
            distHash = {"C1C2Dist" : distMAST} 
            hashList.append(distHash)
        masterDict = {}
        for d in hashList:
            masterDict.update(d)
        df_row = pd.DataFrame([masterDict]) 
        featuresMASTDF = pd.concat([featuresMASTDF, df_row], ignore_index=True)

    return featuresMASTDF
def compartmentalization(logDir , outputDir , substrateFile):
    substrateDF = pd.read_csv(substrateFile)
    cols = list(substrateDF.columns)
    boxCols = boxGen(cols)
    while True:
        colId = input(f"Here are the columns for the file {substrateFile}\n {boxCols}\n Enter the ID number for the column corresponding to the SMILES strings: ").strip()
        try:
            smilesId = int(colId)
            smilesCol = cols[smilesId]
            smilesMAST = substrateDF[smilesCol]
            break
        except:
            print("Try again, enter an integer")
    while True:
        colId = input(f"Here are the columns for the file {substrateFile}\n {boxCols}\n Enter the ID number for the column corresponding to the substrate ID strings: ").strip()
        try:
            substrateId = int(colId)
            subId = cols[substrateId]
            idMAST = substrateDF[subId]
            break
        except:
            print("Try again, enter an integer")
    logPaths = Path(logDir)
    logFiles = list(logPaths.glob('*.log'))  
    firstLog = logFiles[0] if logFiles else None 
    fileSplit = input(f"{str(firstLog)} Enter the string iteral that seperates the common name with the conf. type : ")

    substrateHash = {}
    smilesHash = {}
    for log in logFiles:
        termError = basicTerm(log , "Error termination" , "Normal termination")
        if not termError:
            substrates = list(substrateHash.keys())
            fileID = log.name.split(fileSplit)[0]
            if fileID in substrates:
                substrateHash[fileID].append(log)
            else:
                substrateHash[fileID] = [log]
                pos = idMAST[idMAST == fileID].index[0]
                smilesHash[fileID] = smilesMAST[pos]
    extractNum = input(f"Enter the number corresponding to which substructre you want to extract information from:\n [0] Alkenes\n")
    if int(extractNum) == 0:
        localStrs = ["C13_shift" , "NBO7" , "fukuiParameters" , "%Vbur" , "EvsZ" , "Dist." ]
        localDescriptorsInput = boxGen(localStrs)
        featureList = listInputs(f"Enter the indexes corresponding to the features you would like to extract\n{localDescriptorsInput}")
        featuresMAST = {}
        for idx in featureList:
            feature = localStrs[int(idx)]
            if feature == "fukuiParameters":
                cationStr = input(f"Please enter the --link-- string for the cationic molecule: ")
                anionStr = input(f"Please enter the --link-- string for teh anionic molecule: ")
                feature["fukuiParameters"] = [cationStr , anionStr]
            elif feature == "NBO7":
                neutralStr = input(f"Please enter the --link-- string for the neutral molecule: ")
                feature["NBO7"] = [neutralStr]
            else:
                featuresMAST[feature] = [feature]

        substratesMAST = getAlkenes(substrateHash , smilesHash , featuresMAST)
        outputFile = Path(outputDir) / "alkeneFeaturesMAST.csv"
        substratesMAST.to_csv(outputFile , index=False )

if __name__ == "__main__":
    logDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir ): 
        os.makedirs(outputDir) 
    substrateCSV = str(sys.argv[3])
    compartmentalization(logDir , outputDir, substrateCSV)
