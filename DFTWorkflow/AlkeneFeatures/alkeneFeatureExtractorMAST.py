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
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
from DFTWorkflow.cleanLogs import basicTerm
from DFTWorkflow.ionComGenerator import locateinLog
from DFTWorkflow.AlkeneFeatures.alkeneStericDougnut import alkeneSemiCylinders
from DFTWorkflow.AlkeneFeatures.alkeneSubstitution import eVszAlkenes
from dimensionalityReduction.reactivityFeatures import boxGen
from breadthFirstSearch.radialBasedCorrelation import getCC
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
def extractNBOOccupancies(logFile , nboStr , charge:int):
    logName = str(logFile.name.split(".")[0])
    #print(logName)
    nboIdx = locateinLog(logFile , logName + f"{nboStr}" , 1 )
    nboIdx2 = locateinLog(logFile , logName + f"{nboStr}" , 2 )
    nbo7Hash = {}
    bondsHash = {}
    if charge == 0:
        with open(logFile, 'r') as f:
            extractingDensities = False
            getBondOccupancies = False
            for i, line in enumerate(islice(f, nboIdx, nboIdx2), start=nboIdx):
                if "Atom No    Charge" in line and "Density" not in line:
                    extractingDensities = True 
                if extractingDensities:
                    if "* Total * " in line:
                        extractingDensities = False
                    else:
                        try:
                            occupancies = line.strip()
                            occupancyLines = occupancies.split("    ")
                            atomNums = occupancyLines[0].split("   ")
                            atomInd = atomNums[0].split()
                            nums = []
                            for num in occupancyLines[1:]:
                                num = float(''.join(num.split()))
                                nums.append(num)
                            nbo7Hash[int(atomInd[-1])] = [str(atomInd[0]) , np.max(nums)]
                        except:
                            continue
                elif "Molecular unit" in line:
                    getBondOccupancies = True

                if getBondOccupancies:
                    if ". RY " in line:
                        getBondOccupancies = False
                    if " BD " in line or " BD*"  in line: #collect bonding occupancies
                        lineList = line.strip().split("    ")
                        occupancy = float(lineList[1].strip())
                        energy = float(lineList[-1].split("  ")[0])
                        bondID = int(lineList[0].split("(")[-1].split(")")[0].strip())
                        atoms = lineList[0].split(")")[-1].split("-")
                        atom1 = "".join(atoms[0].split())
                        atom2 = "".join(atoms[1].split())
                        bondHash = {"occupancy" : occupancy , "bondType" : bondID , "atoms" : [atom1 , atom2 ] , "energy" : energy}
                        if " BD " in line:
                            bondType = "bonding_"
                        elif " BD*" in line:
                            bondType = "antibonding_"
                        bondStr = bondType + atom1 + "_" + atom2 + "_" + str(bondID)
                        bondsHash[bondStr] = bondHash
                            
    elif charge != 0:
        with open(logFile, 'r') as f:
            extractingDensities = False
            getBondOccupancies = False
            for i, line in enumerate(islice(f, nboIdx, nboIdx2), start=nboIdx):
                if "Atom No    Charge" in line and "Density" in line:
                    extractingDensities = True 
                if extractingDensities:
                    if "* Total * " in line:
                        extractingDensities = False
                    else:
                        try:
                            occupancies = line.strip()
                            occupancyLines = occupancies.split("    ")
                            atomNums = occupancyLines[0].split("   ")
                            atomInd = atomNums[0].split()
                            nums = []
                            for num in occupancyLines[1:]:
                                num = float(''.join(num.split()))
                                nums.append(num)
                            nbo7Hash[int(atomInd[-1])] = [str(atomInd[0]) , np.max(nums)] 
                        except:
                            continue
    return nbo7Hash, bondsHash
def getAlkeneNBOInfo(logList , C1 , C2 , energyStr , logNameMAST , smiles , nboStr , charge, logEnergyStr):
    weightsDF = getBoltzmannWeightsGauss(logList , 298 , energyStr , logEnergyStr)
    Cmin_NBO = []
    Cmax_NBO = []
    NBO_delta = []
    NBO_mean = []
    piBond = []
    piAntiBond = []
    piEnergy = []
    antiPiEnergy = []
    piBondEnergyWeights = []
    localWeights = weightsDF["boltzWeights"]
    for idx , row in weightsDF.iterrows():
        logFile = row["logID"]
        logPath = f"{logFile}.log"
        conformer  = next((f for f in logList if logPath in f.name), None)
        chargeHash  , bondHash = extractNBOOccupancies(conformer , nboStr , charge)
        alkeneHash = {key: chargeHash[key] for key in chargeHash if key == C1 or key == C2}
        c1NBO = float(alkeneHash[C1][1])
        c2NBO = float(alkeneHash[C2][1])
        if c2NBO >= c1NBO:
            Cmax_NBO.append(c2NBO)
            Cmin_NBO.append(c1NBO)
        else:
            Cmax_NBO.append(c1NBO)
            Cmin_NBO.append(c2NBO)

        NBO_delta.append(abs(c1NBO-c2NBO))
        NBO_mean.append(float(np.mean([c1NBO,c2NBO])))
        if charge == 0:
            try:
                alkeneBondStr = f"bonding_C{C1}_C{C2}_2"
                alkeneAntiBondStr = f"antibonding_C{C1}_C{C2}_2"
                piOccupancy = float(bondHash[alkeneBondStr]["occupancy"])
                piEnergy_ = float(bondHash[alkeneBondStr]["energy"])
                piBond.append(piOccupancy)
                piEnergy.append(piEnergy_)
                piAntiOccupancy = float(bondHash[alkeneAntiBondStr]["occupancy"])
                piAntiEnergy = float(bondHash[alkeneAntiBondStr]["energy"])
                piAntiBond.append(piAntiOccupancy)
                antiPiEnergy.append(piAntiEnergy)
                localWeight = localWeights[idx]
                piBondEnergyWeights.append(localWeight)
            except:
                continue
    def weightedAvg(values, weights):
        weightTot = sum(weights)
        if weightTot == 0:
            return 0
        return sum(v * w for v, w in zip(values, weights)) / weightTot
    finalHash = {}
    finalHash["ID"] = logNameMAST
    finalHash["SMILES"] = smiles
    finalHash["NBO_mxAlk"] =weightedAvg(Cmin_NBO , list(weightsDF["boltzWeights"]))
    finalHash["NBO_mnAlk"] =weightedAvg(Cmax_NBO , list(weightsDF["boltzWeights"]))
    finalHash["NBO_Mean"] =weightedAvg(NBO_mean , list(weightsDF["boltzWeights"]))
    finalHash["NBO_Delta"] =weightedAvg(NBO_delta , list(weightsDF["boltzWeights"]))
    if charge == 0 and len(piBond) != 0:
        finalHash["piBond"] = weightedAvg(piBond, piBondEnergyWeights)
        finalHash["piEnergy"] = weightedAvg(piEnergy, piBondEnergyWeights)
        finalHash["antiPiBond"] = weightedAvg(piAntiBond, piBondEnergyWeights)
        finalHash["antiPiEnergy"] = weightedAvg(antiPiEnergy, piBondEnergyWeights)
    elif charge == 0 and len(piBond) == 0:
        finalHash["piBond"] = "Nan"
        finalHash["piEnergy"] ="Nan"
        finalHash["antiPiBond"] ="Nan"
        finalHash["antiPiEnergy"] = "Nan"  
    return finalHash
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
            atomStr[2:5] = map(float, atomStr[2:5])
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
    firstIdx = locateinLog(logFile , logName + f"{linkStr}" , "earliest" )
    secondIdx = locateinLog(logFile , logName + f"{linkStr}" , "latest" )
    with open(logFile, 'r') as f:
        for i, line in enumerate(f):
            if int(firstIdx) <= i < secondIdx and energyType in line:
                energyLevel = float(line.strip().split("=")[-1].split("A.U.")[0].strip())
                return energyLevel
        return "Poison"
def getBoltzmannWeightsGauss(logDirs, temperature, energyStr , logEnergyStr):
    R = 1.987204e-3
    HARTREE_TO_KCAL = 627.5094740631
    results = []
    
    for log in logDirs:
        if not log.exists():
            print(f"Warning: File {str(log)} not found. Skipping.")
            continue        
        try:
            energy = extractEnergies(log, logEnergyStr, energyStr)
            if energy != "Poison":
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

def getAlkenes(substratesHash , smilesHash , featureHash, logEnergyStr ):
    featuresMASTDF = pd.DataFrame()
    for id , smilesStr in smilesHash.items():
        print(id,smilesStr)
        hashList = []
        cc , molec = getCC(smilesStr)
        conformerFiles = substratesHash[id]
        boltzmannDF = getBoltzmannWeightsGauss(conformerFiles, 298, "electronic" , logEnergyStr)
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
            C13_delta = abs(C13_Cmx - C13_Cmn)
            C13_Mean = np.mean([C13_Cmx,C13_Cmn])
            c13Hash = {"C13_Cmx" : C13_Cmx , "C13_Cmn" : C13_Cmn , "C13_delta" : C13_delta , "C13_Mean" : C13_Mean}
            hashList.append(c13Hash)
        else:
            C1 = cc[1] + 1
            C2 = cc[0] + 1 

        if "NBO7" in featureList:
            nboNeutralStr = featureHash["NBO7"][0]

            neutralHash = getAlkeneNBOInfo(conformerFiles , min([C1,C2]) , max([C1,C2]) , "electronic", id , smilesStr , nboNeutralStr , 0 , logEnergyStr)
            hashList.append(neutralHash)
        if "fukuiParameters" in featureList:
            nboNeutralStr = featureHash["NBO7"][0]
            nboCationStr = featureHash["fukuiParameters"][0]
            nboAnionstr = featureHash["fukuiParameters"][1]
            neutralHash = getAlkeneNBOInfo(conformerFiles , min([C1,C2]) , max([C1,C2]) , "electronic", id , smilesStr , nboNeutralStr , 0, logEnergyStr)
            cationHash = getAlkeneNBOInfo(conformerFiles , min([C1,C2]) , max([C1,C2]) , "electronic", id , smilesStr , nboCationStr , 1 , logEnergyStr)
            anionHash = getAlkeneNBOInfo(conformerFiles , min([C1,C2]) , max([C1,C2]) , "electronic", id , smilesStr, nboAnionstr , -1, logEnergyStr)
            q0_mx = neutralHash["NBO_mxAlk"]
            q1_mx = cationHash["NBO_mxAlk"]
            q1__mx = anionHash["NBO_mxAlk"]
            q0_mn = neutralHash["NBO_mnAlk"]
            q1_mn = cationHash["NBO_mnAlk"]
            q1__mn = anionHash["NBO_mnAlk"]
            f_minus_Mx = float(q0_mx) - float(q1__mx)
            f_plus_Mx = float(q1_mx) - float(q0_mx)
            f_neut_Mx = 0.5 * ((q1_mx) - float(q1__mx))
            f_minus_Mn = float(q0_mn) - float(q1__mn)
            f_plus_Mn = float(q1_mn) - float(q0_mn)
            f_neut_Mn = 0.5 * ((q1_mn) - float(q1__mn))
            hashList.append({"f_neg_mxAlk" : f_minus_Mx , "f_pos_mxAlk" : f_plus_Mx , "f_neut_mxAlk" : f_neut_Mx ,
                            "f_neg_mnAlk" : f_minus_Mn , "f_pos_mnAlk" : f_plus_Mn , "f_neut_mnAlk" : f_neut_Mn , 
                            "f_neg_Delta" : abs(f_minus_Mx - f_minus_Mn) , "f_pos_Delta" : abs(f_plus_Mx-f_plus_Mn) , "f_neut_Delta" : abs(f_neut_Mx-f_neut_Mn) ,
                            "f_neg_meanAlk" : float(f_minus_Mx+f_minus_Mn)/2, "f_pos_meanAlk" : float(f_plus_Mx+f_plus_Mn)/2 , "f_neut_meanAlk" : float(f_neut_Mx+f_neut_Mn)/2 })
        if "%Vbur" in featureList:

            radList = [1.5,2.0,2.5,3.0,3.5,4.0]
            weights = boltzmannDF["boltzWeights"].to_numpy()
            weight_sum = weights.sum()

            # Storage
            Vbur_Cmin = {r: [] for r in radList}
            Vbur_Cmax = {r: [] for r in radList}

            for name in boltzmannDF["logID"]:
                fileStr = f"{name}.log"
                conformer = next(
                    (f for f in conformerFiles if fileStr in f.name),
                    None
                )
                if conformer is None:
                    continue

                coordHash = getAtomCoordsRobust(
                    str(conformer),
                    "GINC-COMPUTE",
                    5,
                    1
                )
                elements = []
                coordinates = []
                for _, coords in coordHash.items():
                    elements.append(str(coords[0]))
                    coordinates.append(np.array(coords[2:5]))

                for rad in radList:
                    vburC1 = BuriedVolume(
                        elements,
                        coordinates,
                        int(C1 - 1),
                        include_hs=True,
                        radius=rad
                    ).fraction_buried_volume

                    vburC2 = BuriedVolume(
                        elements,
                        coordinates,
                        int(C2 - 1),
                        include_hs=True,
                        radius=rad
                    ).fraction_buried_volume

                    Vbur_Cmin[rad].append(vburC1)
                    Vbur_Cmax[rad].append(vburC2)

            # Boltzmann-weighted averages
            BurVolHash = {}
            for rad in radList:
                Vbur_Cmin_arr = np.array(Vbur_Cmin[rad])
                Vbur_Cmax_arr = np.array(Vbur_Cmax[rad])
                Vburr_Cmin = (Vbur_Cmin_arr * weights).sum() / weight_sum
                Vburr_Cmax = (Vbur_Cmax_arr * weights).sum() / weight_sum
                Vburr_mean = (Vburr_Cmin + Vburr_Cmax)/2
                Vburr_delta = abs(Vburr_Cmin - Vburr_Cmax)
                BurVolHash[f"{rad}_Ang_Vburr_Cmn"] = Vburr_Cmin
                BurVolHash[f"{rad}_Ang_Vburr_Cmx"] = Vburr_Cmax
                BurVolHash[f"{rad}_Ang_Vburr_mean"] = Vburr_mean
                BurVolHash[f"{rad}_Ang_Vburr_delta"] = Vburr_delta

            hashList.append(BurVolHash)

        if "%VburSemiCylinders" in featureList:
            radList = [1.5,2.0,2.5,3.0,3.5,4.0]

            maxSemiHash = {r: [] for r in radList}
            minSemiHash = {r: [] for r in radList}
            maxSemiHashOrth = {r: [] for r in radList}
            minSemiHashOrth = {r: [] for r in radList}         
            CburrHash = {r: [] for r in radList}
            g = Graph()
            for bond in molec.GetBonds():
                start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                g.add_edge(start, end)
            CminNeighbors = list(g.neighbors(int(C1-1)))
            CminNeighbors.remove(int(C2-1))
            CmaxNeighbors = list(g.neighbors(int(C2-1)))
            CmaxNeighbors.remove(int(C1-1))
            #print(CmaxNeighbors)
            #print(CminNeighbors)
            #C1Hash = {"0" : [x,y,z] , "1" : [[x,y,z] , [x,y,z]]  , "idx" : idx}
            for name in list(boltzmannDF["logID"]):
                CmaxContacts = []
                CminContacts = []
                fileStr = f"{name}.log"
                conformer  = next((f for f in conformerFiles if fileStr in f.name), None)

                coordHash = getAtomCoordsRobust(str(conformer) , "GINC-COMPUTE" , 5  , 1)
                idx = 0
                for _ , coords in coordHash.items():
                    idx +=1
                    if idx == C1:
                        c1_coords = np.array(coords[2:5])
                        c1Idx = idx
                    elif idx == C2:
                        c2_coords = coords[2:5]
                        c2Idx = idx 
                    if (idx-1) in CminNeighbors:
                        crds = np.array(coords[2:5])
                        CminContacts.append(crds) 
                    if (idx-1) in CmaxNeighbors:
                        crds = np.array(coords[2:5])
                        CmaxContacts.append(crds) 
                C1Hash = {"0" : c1_coords , "1" : CminContacts , "idx" : c1Idx}
                C2Hash = {"0" : c2_coords , "1" : CmaxContacts , "idx" : c2Idx}
                
                for rad in radList:
                    mainCylinder = alkeneSemiCylinders(C1Hash , C2Hash , rad)
                    mainCylinder.getAtoms(coordHash ,{"Nan" : "Nan"} , False)
                    maxSemi_Pi , minSemi_Pi , maxSemi_Orth, minSemi_Orth = mainCylinder.getBurriedVolume(True , False)
                    #Cburr = mainCylinder.getBurriedVolume(False , False)
                    #CburrHash[rad].append(Cburr[0])
                    maxSemiHash[rad].append(maxSemi_Pi)
                    minSemiHash[rad].append(minSemi_Pi)
                    maxSemiHashOrth[rad].append(maxSemi_Orth)
                    minSemiHashOrth[rad].append(minSemi_Orth)
                    

            weights = boltzmannDF["boltzWeights"]
            weight_sum = weights.sum()

            Vbur_MaxSemi = { r: (maxSemiHash[r] * weights).sum() / weight_sum for r in radList}
            Vbur_MinSemi = { r: (minSemiHash[r] * weights).sum() / weight_sum for r in radList}
            Vbur_MaxSemi_Orth = { r: (maxSemiHashOrth[r] * weights).sum() / weight_sum for r in radList}
            Vbur_MinSemi_Orth = { r: (minSemiHashOrth[r] * weights).sum() / weight_sum for r in radList}
            #Vburr_Cylinder = { r: (CburrHash[r] * weights).sum() / weight_sum  for r in radList } 
            Vburr_SemiCylinders ={}
            for rad in radList:
                Vburr_SemiCylinders[f"Vbur_MaxSemi_Pi_{rad}"] = Vbur_MaxSemi[rad]
                Vburr_SemiCylinders[f"Vbur_MinSemi_Pi_{rad}"] = Vbur_MinSemi[rad]
                Vburr_SemiCylinders[f"Vbur_MaxSemi_Orth_{rad}"] = Vbur_MaxSemi_Orth[rad]
                Vburr_SemiCylinders[f"Vbur_MinSemi_Orth_{rad}"] = Vbur_MinSemi_Orth[rad]

            hashList.append(Vburr_SemiCylinders)

        if "EvsZ" in featureList:
            conformer = conformerFiles[0]
            coordHash = getAtomCoordsRobust(str(conformer) , "GINC-COMPUTE" , 5 , 1 )
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

                coordHash = getAtomCoordsRobust(str(conformer) , "GINC-COMPUTE" , 5 , 1 )
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
    logFiles = [
        p for p in logPaths.glob('*.log')
        if '.' not in p.stem]
    firstLog = logFiles[0] if logFiles else None 
    fileSplit = input(f"{str(firstLog)} Enter the string iteral that seperates the common name with the conf. type : ")

    substrateHash = {}
    smilesHash = {}
    for log in logFiles:
        termError = basicTerm(log , "Error termination" , "Normal termination")
        if not termError:
            substrates = list(substrateHash.keys())
            fileID = log.name.split(fileSplit)[0]
            try:
                pos = idMAST[idMAST == fileID].index[0]
                if fileID in substrates:
                    substrateHash[fileID].append(log)
                else:
                    substrateHash[fileID] = [log]
                    pos = idMAST[idMAST == fileID].index[0]
                    smilesHash[fileID] = smilesMAST[pos]
            except:
                continue
    extractNum = input(f"Enter the number corresponding to which substructre you want to extract information from:\n [0] Alkenes\n")
    if int(extractNum) == 0:
        localStrs = ["C13_shift" , "NBO7" , "fukuiParameters" , "%Vbur" , "EvsZ" , "Dist.", "%VburSemiCylinders" ]
        localDescriptorsInput = boxGen(localStrs)
        featureList = listInputs(f"Enter the indexes corresponding to the features you would like to extract\n{localDescriptorsInput}")
        featuresMAST = {}
        for idx in featureList:
            feature = localStrs[int(idx)]
            if feature == "fukuiParameters":
                cationStr = input(f"Please enter the --link-- string for the cationic molecule: ")
                anionStr = input(f"Please enter the --link-- string for the anionic molecule: ")
                featuresMAST["fukuiParameters"] = [cationStr , anionStr]
            elif feature == "NBO7":
                neutralStr = input(f"Please enter the --link-- string for the neutral molecule: ")
                featuresMAST["NBO7"] = [neutralStr]
            else:
                featuresMAST[feature] = [feature]
        logEnergyStr = input(f"Please enter the .log Energy string for these jobs: ")
        substratesMAST = getAlkenes(substrateHash , smilesHash , featuresMAST , logEnergyStr)
        outputFile = Path(outputDir) / "alkeneFeaturesMAST.csv"
        substratesMAST.to_csv(outputFile , index=False )

if __name__ == "__main__":
    logDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir ): 
        os.makedirs(outputDir) 
    substrateCSV = str(sys.argv[3])
    compartmentalization(logDir , outputDir, substrateCSV)
