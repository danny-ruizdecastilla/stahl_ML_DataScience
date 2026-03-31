#Master Feature Extractor
#Danny Ruiz de Castilla | 11.24.2025
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
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
from DFTWorkflow.ionComGenerator import locateinLog
from DFTWorkflow.AlkeneFeatures.alkeneStericDougnut import alkeneSemiCylinders
from DFTWorkflow.AlkeneFeatures.alkeneSubstitution import eVszAlkenes
from DFTWorkflow.dftFeatureExtractorMAST import compartmentalization , time_to_seconds , extractNBOOccupancies , getBoltzmannWeightsGauss , getAtomCoordsRobust , extractShiftsByIdx , getGlobalGreeks
from breadthFirstSearch.radialBasedCorrelation import getCC

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
        chargeHash  , bondHash , lpHash = extractNBOOccupancies(conformer , nboStr , charge)
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
                alkeneBondStr = f"bonding_C{min([C1,C2])}_C{max([C1,C2])}_2"
                alkeneAntiBondStr = f"antibonding_C{min([C1,C2])}_C{max([C1,C2])}_2"
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
            print(C1,C2)
            neutralHash = getAlkeneNBOInfo(conformerFiles , C1 , C2, "electronic", id , smilesStr , nboNeutralStr , 0 , logEnergyStr)
            hashList.append(neutralHash)
        if "fukuiParameters" in featureList:
            nboNeutralStr = featureHash["NBO7"][0]
            nboCationStr = featureHash["fukuiParameters"][0]
            nboAnionstr = featureHash["fukuiParameters"][1]
            neutralHash = getAlkeneNBOInfo(conformerFiles ,  C1 , C2, "electronic", id , smilesStr , nboNeutralStr , 0, logEnergyStr)
            cationHash = getAlkeneNBOInfo(conformerFiles ,  C1 , C2 , "electronic", id , smilesStr , nboCationStr , 1 , logEnergyStr)
            anionHash = getAlkeneNBOInfo(conformerFiles ,  C1 , C2,"electronic", id , smilesStr, nboAnionstr , -1, logEnergyStr)
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
                alkeneType = 0 #gem substituted
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

        if "Sterimol" in featureList:
            radList = [3.0,4.0,4.5]
            weights = boltzmannDF["boltzWeights"].to_numpy()
            weight_sum = weights.sum()

            # Storage
            ster_B1 = {r: [] for r in radList}
            ster_L = {r: [] for r in radList}
            ster_B5 = {r: [] for r in radList}

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
                    #print(C1 , C2)
                    sterimol_values = Sterimol(elements, coordinates, int(C1), int(C2)) 
                    sterimol_values.bury(method="delete", sphere_radius=float(rad))
                    L = sterimol_values.L_value
                    B1 = sterimol_values.B_1_value
                    B5 = sterimol_values.B_5_value
                    ster_B1[rad].append(B1) 
                    ster_B5[rad].append(B5)
                    ster_L[rad].append(L)
            
            sterHash = {}
            for rad in radList:
                sterB1 = np.array(ster_B1[rad])
                sterB5 = np.array(ster_B5[rad])
                sterL = np.array(ster_L[rad])
                avgSterB1 = (sterB1 * weights).sum() / weight_sum
                avgSterB5 = (sterB5 * weights).sum() / weight_sum
                avgSterL = (sterL * weights).sum() / weight_sum
                sterHash[f"sterB1_{rad}"] = avgSterB1
                sterHash[f"sterB5_{rad}"] = avgSterB5
                sterHash[f"sterL_{rad}"] = avgSterL

            hashList.append(sterHash)
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
    substratesMAST = getAlkenes(substrateHash , smilesHash , featuresMAST , logEnergyStr)
    outputFile = Path(outputDir) / "alkeneFeaturesMAST.csv"
    substratesMAST.to_csv(outputFile , index=False )
