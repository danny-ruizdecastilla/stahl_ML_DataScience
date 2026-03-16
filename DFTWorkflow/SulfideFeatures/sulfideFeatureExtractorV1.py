#Danny's version of the sulfide feature extractor 
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
from morfeus import SASA , Sterimol, BuriedVolume
from rdkit import Chem
from networkx import Graph
from itertools import islice
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
from DFTWorkflow.ionComGenerator import locateinLog
from DFTWorkflow.dftFeatureExtractorMAST import compartmentalization ,  time_to_seconds , extractNBOOccupancies , getBoltzmannWeightsGauss, getAtomCoordsRobust , extractShiftsByIdx , getGlobalGreeks
def getS(smiles):
    molec = Chem.MolFromSmiles(smiles)
    sulfideIdx = []
    aromaticSulfides = []
    for atom in molec.GetAtoms():
        if atom.GetAtomicNum() == 16:
            aromatic = atom.GetIsAromatic()
            sulfideID = atom.GetIdx()
            if aromatic:
                aromaticSulfides.append(sulfideID)
                continue
            else:
                sulfideIdx.append(sulfideID)
    return sulfideIdx , aromaticSulfides , molec 
def weightedAvg(values, weights):
    weightTot = sum(weights)
    if weightTot == 0:
        return 0
    return sum(v * w for v, w in zip(values, weights)) / weightTot
def getSulfidesNBOInfo(logList , sulfurIdx ,  energyStr , logNameMAST , smiles , nboStr , charge, logEnergyStr):
    weightsDF = getBoltzmannWeightsGauss(logList , 298 , energyStr , logEnergyStr)
    lonePairMax = []
    lonePairMin = []    
    sulfurNBO = []
    lonePairEnergyMax = []
    lonePairEnergyMin = []
    for idx , row in weightsDF.iterrows():
        logFile = row["logID"]
        logPath = f"{logFile}.log"
        conformer  = next((f for f in logList if logPath in f.name), None)
        chargeHash  , bondHash , lonePairHash = extractNBOOccupancies(conformer , nboStr , charge)
        #print(lonePairHash)
        sulfurCharge = float(chargeHash[sulfurIdx][1])
        sulfurNBO.append(sulfurCharge)
        sulfurLP1 = lonePairHash[f"S_{sulfurIdx}_lp_1"]
        sulfurLP2 = lonePairHash[f"S_{sulfurIdx}_lp_2"]
        lonePairMax.append(sulfurLP1["occupancy"])
        lonePairMin.append(sulfurLP2["occupancy"])
        lonePairEnergyMax.append(sulfurLP1["energy"])
        lonePairEnergyMin.append(sulfurLP2["energy"])
    finalHash = {}
    finalHash["ID"] = logNameMAST
    finalHash["SMILES"] = smiles
    finalHash["SulfurNBO"] =weightedAvg(sulfurNBO , list(weightsDF["boltzWeights"]))
    finalHash["lonePair1"] =weightedAvg(lonePairMax , list(weightsDF["boltzWeights"]))
    finalHash["lonePair2"] =weightedAvg(lonePairMin, list(weightsDF["boltzWeights"]))
    finalHash["lonePairEnergy1"] =weightedAvg(lonePairEnergyMax , list(weightsDF["boltzWeights"]))
    finalHash["lonePairEnergy2"] =weightedAvg(lonePairEnergyMin  , list(weightsDF["boltzWeights"]))

    return finalHash

def getSulfides(substratesHash , smilesHash , featureHash, logEnergyStr ):
    featuresMASTDF = pd.DataFrame()
    for id , smilesStr in smilesHash.items():
        print(id,smilesStr)
        hashList = []
        sulfides , aromaticSulfides , molec  = getS(smilesStr)
        conformerFiles = substratesHash[id]
        boltzmannDF = getBoltzmannWeightsGauss(conformerFiles, 298, "electronic" , logEnergyStr)
        idHash = {"SMILES" : smilesStr , "ID" : id}  
        hashList.append(idHash)   
        featureList = list(featureHash.keys())  
        #["SulfurShift" , "NBO7" , "fukuiParameters" , "%Vbur" ,  "globalFeatures" , "Sterimol" ]
        if "SulfurShift" in featureList:
            sulfurShifts = []
            for name in list(boltzmannDF["logID"]):
                fileStr = f"{name}.log"
                conformer  = next((f for f in conformerFiles if fileStr in f.name), None)
                sulfurShiftsHash = extractShiftsByIdx(conformer, "SCF GIAO Magnetic shielding tensor (ppm):", "Eigenvalues:  ", "earliest", "latest", 
                                                      [s +1 for s in sulfides], 1, 0)
                sulfurShifts.append(sulfurShiftsHash[sulfides[0] + 1][1])
            weights = boltzmannDF["boltzWeights"]
            sulfurIsotropic = ( (sulfurShifts * weights).sum()    / weights.sum()  )
            sulfurIsotropicHash = {"sulfurIsotropicShift" : sulfurIsotropic}
            hashList.append(sulfurIsotropicHash)

        if "NBO7" in featureList:
            nboNeutralStr = featureHash["NBO7"][0]

            neutralHash = getSulfidesNBOInfo(conformerFiles ,  [s +1 for s in sulfides][0], "electronic", id , smilesStr , nboNeutralStr , 0 , logEnergyStr)
            hashList.append(neutralHash)

        if "fukuiParameters" in featureList:
            nboNeutralStr = featureHash["NBO7"][0]
            nboCationStr = featureHash["fukuiParameters"][0]
            nboAnionStr = featureHash["fukuiParameters"][1]
            fukPos = []
            fukNeg = []
            fukNeut = []
            for name in list(boltzmannDF["logID"]):
                fileStr = f"{name}.log"
                conformer  = next((f for f in conformerFiles if fileStr in f.name), None)
                chargeHash  , bondHash , lonePairHash = extractNBOOccupancies(conformer , nboNeutralStr , 0)
                neutralCharge = float(chargeHash[(sulfides[0]+1)][1])
                chargeHash  , bondHash , lonePairHash = extractNBOOccupancies(conformer , nboCationStr, 1)
                cationCharge = float(chargeHash[(sulfides[0]+1)][1])
                chargeHash  , bondHash , lonePairHash = extractNBOOccupancies(conformer , nboAnionStr , -1)
                anionCharge = float(chargeHash[(sulfides[0]+1)][1])

                f_minus = float(neutralCharge) - float(anionCharge)
                f_plus= float(cationCharge) - float(neutralCharge)
                f_neut= 0.5 * ((cationCharge) - float(anionCharge))
                fukPos.append(f_plus)
                fukNeg.append(f_minus)
                fukNeut.append(f_neut)

            weights = boltzmannDF["boltzWeights"]
            avgFukPos = weightedAvg(fukPos ,weights )
            avgFukNeg = weightedAvg(fukNeg ,weights )
            avgFukNeut = weightedAvg(fukNeut ,weights )

            hashList.append({"f_neg_Sulfide" : avgFukNeg , "f_pos_Sulfide" : avgFukPos, "f_neut_mxAlk" : avgFukNeut })

        if "%Vbur" in featureList:

            radList = [1.5,2.0,2.5,3.0,3.5,4.0]
            weights = boltzmannDF["boltzWeights"].to_numpy()
            weight_sum = weights.sum()

            # Storage
            Vbur_Sulfur = {r: [] for r in radList}

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
                    vburS = BuriedVolume(
                        elements,
                        coordinates,
                        int(sulfides[0]),
                        include_hs=True,
                        radius=rad
                    ).fraction_buried_volume

                    Vbur_Sulfur[rad].append(vburS)

            # Boltzmann-weighted averages
            BurVolHash = {}
            for rad in radList:
                Vbur_Sulfur_arr = np.array(Vbur_Sulfur[rad])
                Vburr_Sulf = (Vbur_Sulfur_arr * weights).sum() / weight_sum
                BurVolHash[f"{rad}_Ang_Vburr_Sulfure"] = Vburr_Sulf

            hashList.append(BurVolHash)

        if "globalFeatures" in featureList:
            nboNeutralStr = featureHash["globalFeatures"][0]
            features = {
                "sasaVol": [],
                "sasaArea": [],
                "sasaSurface": [],
                "homo": [],
                "lumo": [],
                "eta": [],
                "mu": [],
                "omega": [],
                "cpuTime": [],
                "wallTime": [],
                "dipole": [],
                "polarizability": []
            }
            cpuStr = "Job cpu time:"
            walltimeStr = "Elapsed time:"
            weights = boltzmannDF["boltzWeights"]
            for name in list(boltzmannDF["logID"]):
                fileStr = f"{name}.log"
                conformer  = next((f for f in conformerFiles if fileStr in f.name), None)
                coordHash = getAtomCoordsRobust(str(conformer) , "GINC-COMPUTE" , 5 , 1 )
                elements = []
                coordinates = []
                for _, coords in coordHash.items():
                    elements.append(str(coords[0]))
                    coordinates.append(np.array(coords[2:5]))
                sasa_ = SASA(elements,coordinates) 
                sphericity = np.cbrt((36*np.pi*sasa_.volume**2))/sasa_.area
                area = sasa_.area
                vol = sasa_.volume
                homoE , lumoE,  muMolec , etaMolec , omegaMolec = getGlobalGreeks(conformer ,nboNeutralStr  , "Alpha  occ. eigenvalues --" , "Alpha virt. eigenvalues --") 
                features["sasaVol"].append(vol)
                features["sasaArea"].append(area)
                features["sasaSurface"].append(sphericity)
                features["homo"].append(homoE)
                features["lumo"].append(lumoE)
                features["mu"].append(muMolec)
                features["eta"].append(etaMolec)
                features["omega"].append(omegaMolec)
                cpuIdx = locateinLog(str(conformer) , cpuStr , 0 )
                wallTimeIdx = locateinLog(str(conformer) , walltimeStr , 0 )
                dipoleIdx = locateinLog(str(conformer) , "Dipole moment (field-independent basis, Debye):" , 0)
                polarIdx = locateinLog(str(conformer) , "Approx polarizability:   " , 0)
                with open(str(conformer), "r") as f:
                    lines = f.readlines()

                    cpu_line = lines[cpuIdx].split(cpuStr)[-1]
                    wall_line = lines[wallTimeIdx].split(walltimeStr)[-1]

                    features["cpuTime"].append(time_to_seconds(cpu_line))
                    features["wallTime"].append(time_to_seconds(wall_line))

                    dipoleStr = lines[dipoleIdx + 1].split("Tot=")[-1].strip()
                    features["dipole"].append(float(dipoleStr))

                    polarStr = lines[polarIdx].split("Approx polarizability:   ")[-1]
                    features["polarizability"].append(float(polarStr.split()[0]))
            globalRow = {
                key: np.average(values, weights=weights)
                for key, values in features.items()
            }   
            hashList.append(globalRow)
        masterDict = {}
        for d in hashList:
            masterDict.update(d)
        df_row = pd.DataFrame([masterDict]) 
        featuresMASTDF = pd.concat([featuresMASTDF, df_row], ignore_index=True)

    return featuresMASTDF   


if __name__ == "__main__":
    logDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir ): 
        os.makedirs(outputDir) 
    substrateCSV = str(sys.argv[3])
    substrateHash , smilesHash , featuresMAST , logEnergyStr  = compartmentalization(logDir , outputDir, substrateCSV)
    substratesMAST = getSulfides(substrateHash , smilesHash , featuresMAST , logEnergyStr)
    outputFile = Path(outputDir) / "sulfideFeaturesMAST.csv"
    substratesMAST.to_csv(outputFile , index=False )
