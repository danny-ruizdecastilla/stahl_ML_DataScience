#Danny Ruiz de Castilla
#Generate features from xyz files 
#MM Paper Exclusive
import os
import glob as glob
import numpy as np
import pandas as pd 
from pathlib import Path
import sys
from rdkit import Chem
from itertools import combinations
from morfeus import Sterimol, BuriedVolume
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
from DFTWorkflow.AlkeneFeatures.alkeneStericDougnut import alkeneSemiCylinders
from DFTWorkflow.AlkeneFeatures.alkeneSlicedOranges import alkeneSlicedOranges
sys.path.append(str(parentDir))
def vdw_radius(symbol):
    pt = Chem.GetPeriodicTable()
    return pt.GetRvdw(pt.GetAtomicNumber(symbol))
def alkene2xyz(xyzFile):
    xyz = Path(xyzFile)
    atomHash = {}
    atomCoords = []
    carbonList = []
    with open(xyz, "r") as f:
        for lineNum, line in enumerate(f, start=1):
            if lineNum == 1:
                numAtoms = int(line.strip())
                #atomHash["totAtoms"] = numAtoms
            elif "." in line: #atom coords
                atomIdx = int(lineNum - 3) 
                parts = line.strip().split()
                atomHash[atomIdx] = [str(parts[0])] + [0] +  [float(part) for part in parts[1:4]]
                atomCoords.append(np.array([float(part) for part in parts[1:4]]))
                if str(parts[0]) == "C":
                    carbonList.append(atomIdx)
    coords = np.array(atomCoords)
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]

    distanceMatrix = np.linalg.norm(diff, axis=2)
    existingHash = {"wiscAlks_193" :[(0, 1, np.float64(1.32350208774977))] , "Alkenes_02863" : [(0, 1, np.float64(1.3278428512885099))] , 
                    "Alkenes_02868" : [(0, 1, np.float64(1.3279667991011674))], "Alkenes_02870": [(0, 1, np.float64(1.3271562010984237))],
                      "Alkenes_02878": [(0, 1, np.float64(1.3271024689860988))], "Alkenes_02893" : [(0, 1, np.float64(1.3272129846433844))] }
    idName = xyz.name.split(".")[0]
    if idName in existingHash.keys():
        alkenes = existingHash[idName]
        return alkenes , atomHash , distanceMatrix

    else:
        alkeneMax = 1.430417269
        alkeneMin = 1.308758104
        alkenes = []
        ccDoubleTraversed = {}
        for i, j in combinations(carbonList, 2):
            d = distanceMatrix[i, j]
            if alkeneMin <= d <= alkeneMax:
                alkenes.append((i, j, d))
                ccDoubleTraversed[i] = ccDoubleTraversed.get(i, 0) + 1
                ccDoubleTraversed[j] = ccDoubleTraversed.get(j, 0) + 1
        #print(ccDoubleTraversed)
        #Remove cyclical fragments from heterocycles 
        
        if len(alkenes) != 1:
            #we need to do some traversals to remove aromatics 
            removeAroms = []
            for bond in alkenes:
                if ccDoubleTraversed[bond[0]] >= 2 or ccDoubleTraversed[bond[1]] >=2:
                    continue
                else:
                    removeAroms.append(bond)
            return removeAroms ,  atomHash , distanceMatrix
        else:    
            return alkenes , atomHash , distanceMatrix
def alkene2Feats(xyzList):
    featuresMASTDF = pd.DataFrame()
    for xyz in xyzList:
        idName = Path(xyz).name.split(".")[0]
        print(idName)
        hashList = [{"ID" : idName}]
        alkene , atomHash , distanceMatrix = alkene2xyz(xyz)
        #Get alkene Contacts 
        #CHash = {"0" : c1_coords , "1" : CminContacts , "idx" : c1Idx}
        CAlphaHash = {}
        CBetaHash = {}

        for bond in alkene: 
            source = bond[0]
            target = bond[1]
            allLengths = distanceMatrix[source].copy()
            allLengths[[source, target]] = np.inf
            nearestSource = [(int(idx), float(allLengths[idx])) for idx in np.argsort(allLengths)[:2]]
            sourceEndsVdw = 0 
            for atom in nearestSource:
                atomID = atomHash[atom[0]][0]

                sourceEndsVdw += vdw_radius(atomID)

            allLengths2 = distanceMatrix[target].copy()
            allLengths2[[target, source]] = np.inf
            nearestTarget = [(int(idx), float(allLengths2[idx])) for idx in np.argsort(allLengths2)[:2]]
            targetEndsVdw = 0 
            for atom in nearestTarget:
                atomID = atomHash[atom[0]][0]
                targetEndsVdw += vdw_radius(atomID)
            if targetEndsVdw > sourceEndsVdw:
                #Target is alpha Source is beta 
                CAlphaHash["0"] = np.array(atomHash[target][2:5])
                CAlphaHash["idx"] = target 
                CAlphaHash["1"] = [np.array(atomHash[nearestTarget[0][0]][2:5]) , np.array(atomHash[nearestTarget[1][0]][2:5])]

                CBetaHash["0"] = np.array(atomHash[source][2:5])
                CBetaHash["idx"] = source

                length1 = np.linalg.norm(np.array(atomHash[nearestTarget[0][0]][2:5])  - np.array(atomHash[nearestSource[0][0]][2:5]))
                length2 = np.linalg.norm(np.array(atomHash[nearestTarget[0][0]][2:5])  - np.array(atomHash[nearestSource[1][0]][2:5]))
                if length1 > length2:
                    CBetaHash["1"] = [np.array(atomHash[nearestSource[1][0]][2:5]) , np.array(atomHash[nearestSource[0][0]][2:5])]
                else:
                    CBetaHash["1"] = [np.array(atomHash[nearestSource[0][0]][2:5]) , np.array(atomHash[nearestSource[1][0]][2:5])]
            else:
                #Target is beta Source is alpha 
                CAlphaHash["0"] = np.array(atomHash[source][2:5])
                CAlphaHash["idx"] = source
                CAlphaHash["1"] = [np.array(atomHash[nearestSource[0][0]][2:5]) , np.array(atomHash[nearestSource[1][0]][2:5])]

                CBetaHash["0"] = np.array(atomHash[target][2:5])
                CBetaHash["idx"] = target 
                length1 = np.linalg.norm(np.array(atomHash[nearestSource[0][0]][2:5])  - np.array(atomHash[nearestTarget[0][0]][2:5]))
                length2 = np.linalg.norm(np.array(atomHash[nearestSource[0][0]][2:5])  - np.array(atomHash[nearestTarget[1][0]][2:5]))
                if length1 > length2:
                    CBetaHash["1"] = [np.array(atomHash[nearestTarget[1][0]][2:5]) , np.array(atomHash[nearestTarget[0][0]][2:5])]
                else:
                    CBetaHash["1"] = [np.array(atomHash[nearestTarget[0][0]][2:5]) , np.array(atomHash[nearestTarget[1][0]][2:5])]
        #Done orienting lets calculate features 
        radListSemis = [2.5,3.0,3.5]
        segmentList = ["maxSemiPi" , "minSemiPi" , "maxSemiOrth" , "minSemiOrth" , "maxCap" , "minCap" , "CburrTot"]
        cappedSemisHash = {}
        for rad in radListSemis:
            mainCylinder = alkeneSemiCylinders(CAlphaHash , CBetaHash , rad , 0.15)
            #print(atomHash)
            mainCylinder.getAtoms(atomHash ,{"Nan" : "Nan"} , False)
            segmentHash = mainCylinder.getBurriedVolume(True , False)
            for segment in segmentList:
                cappedSemisHash[f"Vbur_{segment}_{rad}"] = segmentHash[segment]
            cappedSemisHash[f"Vbur_deltaPi_{rad}"] = np.abs(cappedSemisHash[f"Vbur_maxSemiPi_{rad}"] - cappedSemisHash[f"Vbur_minSemiPi_{rad}"])
            cappedSemisHash[f"Vbur_deltaOrth_{rad}"] = np.abs(cappedSemisHash[f"Vbur_maxSemiOrth_{rad}"] - cappedSemisHash[f"Vbur_minSemiOrth_{rad}"])
            cappedSemisHash[f"Vbur_deltaCap_{rad}"] = np.abs(cappedSemisHash[f"Vbur_maxCap_{rad}"] - cappedSemisHash[f"Vbur_minCap_{rad}"])
        hashList.append(cappedSemisHash)

        slicesIdx = {"C1" : "CAlpha" , "C2" : "CBeta"}
        atoms = []
        idxList = []
        symbols = []
        for atomIdx , coords in atomHash.items():
            idxList.append(atomIdx)
            atoms.append(coords[2:5])
            symbols.append(coords[0])
        radListSlices = [2.0,2.5,3.0]
        orangeData = {}
        for rad in radListSlices:
            alkeneOranges = alkeneSlicedOranges(CAlphaHash , CBetaHash , rad)
            alkeneOranges.getAtoms(atoms , idxList , symbols ,{"Nan" : "Nan"} , False)
            orangeSlices = alkeneOranges.slicedOranges(True , False , False ,"alkene")
            for center in ("C1", "C2"): #C1 is Alpha C2 is Beta 
                name = slicesIdx[center]
                for quad in range(4):
                    orangeData[f"vbur{name}_slice{quad}_{rad}"] = orangeSlices[f"{center}_quad_{quad}"]
        hashList.append(orangeData)

        firstContactRad = [2.0,2.5,3.0,3.5,4.0]
        firstContactHash = {}
        for rad in firstContactRad:
            cAlphaContacts = CAlphaHash["1"]
            CBetaContacts = CBetaHash["1"]
            for i in range(len(cAlphaContacts)):
                AlphaContact = atoms.index(list(cAlphaContacts[i]))
                BetaContact = atoms.index(list(CBetaContacts[i]))

                alphaContactVbur = BuriedVolume(symbols,atoms,AlphaContact,include_hs=True,radius=rad).fraction_buried_volume
                betaContactVbur = BuriedVolume(symbols,atoms,BetaContact,include_hs=True,radius=rad).fraction_buried_volume

                firstContactHash[f"C_Alpha_{i+1}_{rad}"] = alphaContactVbur
                firstContactHash[f"C_Beta_{i+1}_{rad}"] = betaContactVbur
            avgFace1 = (firstContactHash[f"C_Alpha_1_{rad}"] + firstContactHash[f"C_Beta_1_{rad}"]) /2 
            avgFace2 = (firstContactHash[f"C_Alpha_2_{rad}"] + firstContactHash[f"C_Beta_2_{rad}"]) /2 
            avg1 = (firstContactHash[f"C_Alpha_1_{rad}"] + firstContactHash[f"C_Alpha_2_{rad}"]) /2 
            avg2 = (firstContactHash[f"C_Beta_1_{rad}"] + firstContactHash[f"C_Beta_2_{rad}"]) /2 
            firstContactHash[f"delta_Faces_{rad}"] = np.abs(avgFace1-avgFace2)
            firstContactHash[f"delta_Ends_{rad}"] = np.abs(avg1-avg2)
        hashList.append(firstContactHash)

        sterimolHash = {}
        sterimolRad = [2.0 , 2.5 , 3.0 , 3.5 , 4.0 , 4.5]
        for rad in sterimolRad:
            sterimolAlpha = Sterimol(symbols, atoms, CAlphaHash["idx"] + 1 , CBetaHash["idx"] + 1) 
            sterimolAlpha.bury(method="delete", sphere_radius=float(rad))
            sterimolBeta = Sterimol(symbols, atoms, CBetaHash["idx"] + 1 , CAlphaHash["idx"] + 1) 
            sterimolBeta.bury(method="delete", sphere_radius=float(rad))

            sterimolHash[f"CAlpha_sterimol_{rad}_L"] = sterimolAlpha.L_value
            sterimolHash[f"CAlpha_sterimol_{rad}_B1"] = sterimolAlpha.B_1_value
            sterimolHash[f"CAlpha_sterimol_{rad}_B5"] = sterimolAlpha.B_5_value

            sterimolHash[f"CBeta_sterimol_{rad}_L"] = sterimolBeta.L_value
            sterimolHash[f"CBeta_sterimol_{rad}_B1"] = sterimolBeta.B_1_value
            sterimolHash[f"CBeta_sterimol_{rad}_B5"] = sterimolBeta.B_5_value
        hashList.append(sterimolHash)
        masterDict = {}
        for d in hashList:
            masterDict.update(d)
        df_row = pd.DataFrame([masterDict]) 
        featuresMASTDF = pd.concat([featuresMASTDF, df_row], ignore_index=True)
    return featuresMASTDF
if __name__ == "__main__":
    xyzdir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir ): 
        os.makedirs(outputDir) 
    xyzFiles = glob.glob(xyzdir + "/*.xyz")
    featureDF = alkene2Feats(xyzFiles)
    featureDF.to_csv(f"{outputDir}/MMFeatureMAST.csv")