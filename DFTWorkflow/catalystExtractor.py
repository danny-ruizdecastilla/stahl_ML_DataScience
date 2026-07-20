#Danny Ruiz de Castilla 
#07.04.2026
import os 
import sys
import pandas as pd 
import numpy as np
import chemdraw
import glob as glob
from pathlib import Path
from rdkit import Chem
from networkx import Graph
from morfeus import BuriedVolume
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from DFTWorkflow.cleanLogs import basicTerm
from dimensionalityReduction.reactivityFeatures import boxGen
from DFTWorkflow.dftFeatureExtractorMAST import time_to_seconds , extractNBOOccupancies , getBoltzmannWeightsGauss ,  getAtomCoordsRobust

def locateCatRows(logDir , splitStr , catalystFile):
    catalystDF = pd.read_csv(catalystFile)
    cols = list(catalystDF.columns)
    boxCols = boxGen(cols)
    while True:
        colId = input(f"Here are the columns for the file {catalystFile}\n {boxCols}\n Enter the ID number for the column corresponding to the SMILES strings: ").strip()
        try:
            smilesId = int(colId)
            smilesCol = cols[smilesId]
            smilesMAST = catalystDF[smilesCol]
            break
        except:
            print("Try again, enter an integer")
    while True:
        colId = input(f"Here are the columns for the file {catalystFile}\n {boxCols}\n Enter the ID number for the column corresponding to the catalyst ID strings: ").strip()
        try:
            catalystId = int(colId)
            subId = cols[catalystId]
            idMAST = catalystDF[subId]
            break
        except:
            print("Try again, enter an integer")
    #returns hash table of log paths that lead to unique catalysts with solvents 
    path = Path(logDir)
    logFiles = list(path.glob("*/*.log"))
    catalystHash = {}
    smilesHash = {}
    for log in logFiles:
        termError = basicTerm(log , "Error termination" , "Normal termination")
        if not termError:
            catalysts = list(catalystHash.keys())
            fileID = log.name.split(splitStr)[0]
            solventID = str(log).split("/")[-2]
            try:
                pos = idMAST[idMAST == fileID].index[0]
                if f"{fileID}_{solventID}" in catalysts:
                    catalystHash[f"{fileID}_{solventID}"].append(log)
                else:
                    catalystHash[f"{fileID}_{solventID}"] = [log]
                    pos = idMAST[idMAST == fileID].index[0]
                    smilesHash[fileID] = smilesMAST[pos]
            except:
                continue
    return catalystHash , smilesHash
def getChiralKetone(smiles):
    molec = Chem.MolFromSmiles(smiles)

    for bond in molec.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if (atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE):
            #Carbonyl with C at atom1
            neighbors = atom1.GetNeighbors()
            carbonCount = 0
            for nbr in neighbors:
                atomNum = nbr.GetAtomicNum()
                if int(atomNum) == 6:
                    carbonCount +=1
            if carbonCount == 2 :# 2 carbon neighbors this is the carbonyl we want
                c = bond.GetBeginAtomIdx() #carbon first then oxygen
                o = bond.GetEndAtomIdx()
                return [c , o] , molec
        elif (atom1.GetAtomicNum() == 8 and atom2.GetAtomicNum() == 6 and bond.GetBondType() == Chem.BondType.DOUBLE):
            #Carbonyl with C at atom1
            neighbors = atom2.GetNeighbors()
            carbonCount = 0
            for nbr in neighbors:
                atomNum = nbr.GetAtomicNum()
                if int(atomNum) == 6:
                    carbonCount +=1
            if carbonCount == 2 :# 2 carbon neighbors this is the carbonyl we want
                o = bond.GetBeginAtomIdx() #carbon first then oxygen
                c = bond.GetEndAtomIdx()
                return [c , o] , molec      
    return ["Poison" , "Poison"] , molec
def getKetoneNBOInfo(weightsDF,logList  ,  carbonyl , nboStr , charge , catalystID , smiles ):
    carbonyl = [c + 1 for c in carbonyl]
    cNBO = []
    oNBO = []

    piBond = []
    piAntiBond = []
    piEnergy = []
    antiPiEnergy = []

    oLP1 = []
    oLP2 = []
    oLP1E = []
    oLP2E = []

    piBondEnergyWeights = []
    localWeights = weightsDF["boltzWeights"]

    for idx , row in weightsDF.iterrows():
        logFile = row["logID"]
        logPath = f"{logFile}.log"
        conformer  = next((f for f in logList if logPath in f.name), None)
        chargeHash  , bondHash , lonePairHash = extractNBOOccupancies(conformer , nboStr , charge)
        carbonylHash = {key: chargeHash[key] for key in chargeHash if key == carbonyl[0] or key == carbonyl[1]}

        cNBO.append(float(carbonylHash[carbonyl[0]][1]))
        oNBO.append(float(carbonylHash[carbonyl[1]][1]))
        
        oxygenLP1 = lonePairHash[f"O_{carbonyl[1]}_lp_1"]
        oxygenLP2 = lonePairHash[f"O_{carbonyl[1]}_lp_2"]
        oLP1.append(oxygenLP1["occupancy"])
        oLP2.append(oxygenLP2["occupancy"])
        oLP1E.append(oxygenLP1["energy"])
        oLP2E.append(oxygenLP2["energy"])

        if charge == 0:
            try:
                carbonylBondStr = f"bonding_O{carbonyl[1]}_C{carbonyl[0]}_2"
                carbonylAntiBondStr = f"antibonding_O{carbonyl[1]}_C{carbonyl[0]}_2"
                piOccupancy = float(bondHash[carbonylBondStr]["occupancy"])
                piEnergy_ = float(bondHash[carbonylAntiBondStr]["energy"])
                piBond.append(piOccupancy)
                piEnergy.append(piEnergy_)
                piAntiOccupancy = float(bondHash[carbonylBondStr]["occupancy"])
                piAntiEnergy = float(bondHash[carbonylAntiBondStr]["energy"])
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
    finalHash["ID"] = catalystID
    finalHash["SMILES"] = smiles
    finalHash["carbonNBO"] =weightedAvg(cNBO , list(weightsDF["boltzWeights"]))
    finalHash["oxygenNBO"] =weightedAvg(oNBO , list(weightsDF["boltzWeights"]))
    finalHash["oxygenLonePair1"] =weightedAvg(oLP1 , list(weightsDF["boltzWeights"]))
    finalHash["oxygenLonePair2"] =weightedAvg(oLP2, list(weightsDF["boltzWeights"]))
    finalHash["oxygenLonePairEnergy1"] =weightedAvg(oLP1E, list(weightsDF["boltzWeights"]))
    finalHash["oxygenLonePairEnergy2"] =weightedAvg(oLP2E , list(weightsDF["boltzWeights"]))
    return finalHash

def getChiralKetoneFeatures(catalystHash , smilesHash  , logEnergyStr , nboStr ):
    featuresMASTDF = pd.DataFrame()
    for catalystID , catPaths in catalystHash.items():
        print(catalystID)
        catID = catalystID.split("_")[1]
        idName = f"cat_{catID}"
        solvent = catalystID.split("_")[2]
        smilesStr = smilesHash[idName]
        carbonyl , molec = getChiralKetone(smilesStr)
        boltzmannDF = getBoltzmannWeightsGauss(catPaths, 298, "electronic" , f"{logEnergyStr}_{solvent}")
        #NBO Features 
        neutralHash = getKetoneNBOInfo(boltzmannDF,catPaths ,  carbonyl,  nboStr , 0 , catalystID , smilesStr)
        oVbur_3 = []
        weights = boltzmannDF["boltzWeights"].to_numpy()
        minWeightIdx = max(enumerate(weights), key=lambda x: x[1])[0]
        weight_sum = weights.sum()
        for name in boltzmannDF["logID"]:
            print(name)
            fileStr = f"{name}.log"
            conformer = next((f for f in catPaths if fileStr in f.name),None)
            if conformer is None:
                continue
            coordHash = getAtomCoordsRobust(str(conformer),"GINC-COMPUTE",5,1)
            elements = []
            coordinates = []
            for _, coords in coordHash.items():
                elements.append(str(coords[0]))
                coordinates.append(np.array(coords[2:5]))
            vburO = BuriedVolume(elements,coordinates,int(carbonyl[1]),include_hs=True,radius=3.0).fraction_buried_volume
            oVbur_3.append(vburO)
        avgVbur = (np.asarray(oVbur_3)*weights).sum()/weight_sum 
        Vbur_lowE = oVbur_3[minWeightIdx]

        final = {"Solvent" : solvent , "Ox_vbur_3" : avgVbur , "Ox_vbur_3_lowE" : Vbur_lowE}
        masterDict = {}
        for d in [neutralHash , final]:
            masterDict.update(d)
        df_row = pd.DataFrame([masterDict]) 
        featuresMASTDF = pd.concat([featuresMASTDF, df_row], ignore_index=True)
    return featuresMASTDF


if __name__ == "__main__":
    logDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir ): 
        os.makedirs(outputDir) 
    catalystCSV = str(sys.argv[3])
    splitStr = str(sys.argv[4])
    logEnergyStr = str(sys.argv[5])
    catalystHash , smilesHash  = locateCatRows(logDir , splitStr , catalystCSV)
    substratesMAST = getChiralKetoneFeatures(catalystHash , smilesHash  , logEnergyStr , "nbo7_neutral" )
    outputFile = Path(outputDir) / "catalystFeaturesMAST.csv"
    substratesMAST.to_csv(outputFile , index=False )