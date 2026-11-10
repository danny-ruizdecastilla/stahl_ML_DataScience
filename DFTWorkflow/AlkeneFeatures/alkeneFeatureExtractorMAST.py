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
from DFTWorkflow.AlkeneFeatures.alkeneSlicedOranges import alkeneSlicedOranges
from DFTWorkflow.AlkeneFeatures.alkeneSubstitution import eVszAlkenes
from DFTWorkflow.dftFeatureExtractorMAST import compartmentalization , time_to_seconds , extractNBOOccupancies , getBoltzmannWeightsGauss , getAtomCoordsRobust , extractShiftsByIdx , getGlobalGreeks , extractCHELPGCharges , extractHirshfeldCharges , extractFullCharges
from breadthFirstSearch.radialBasedCorrelation import getCC
from topologicalStericEffects.alkeneTopologyCalc import calcSingle
def weightedAvg(values, weights):
    weightTot = sum(weights)
    if weightTot == 0:
        return 0
    return sum(v * w for v, w in zip(values, weights)) / weightTot
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
    finalHash = {}
    finalHash["ID"] = logNameMAST
    finalHash["SMILES"] = smiles
    finalHash["NBO_mxAlk"] =weightedAvg(Cmin_NBO , list(weightsDF["boltzWeights"]))
    finalHash["NBO_mnAlk"] =weightedAvg(Cmax_NBO , list(weightsDF["boltzWeights"]))
    finalHash["NBO_Mean"] =weightedAvg(NBO_mean , list(weightsDF["boltzWeights"]))
    finalHash["NBO_Delta"] =weightedAvg(NBO_delta , list(weightsDF["boltzWeights"]))

    lowEidx = max(enumerate(list(weightsDF["boltzWeights"])), key=lambda x: x[1])[0]
    Cmin_NBO_lowE = Cmin_NBO[lowEidx]
    Cmax_NBO_lowE = Cmax_NBO[lowEidx]
    finalHash["NBO_mxAlk_lowE"] = Cmin_NBO_lowE
    finalHash["NBO_mnAlk_lowE"] =Cmax_NBO_lowE
    finalHash["NBO_Mean_lowE"] = np.mean([Cmin_NBO_lowE, Cmax_NBO_lowE])
    finalHash["NBO_Delta_lowE"] =abs(Cmin_NBO_lowE - Cmax_NBO_lowE)
    if charge == 0 and len(piBond) != 0:
        piBondlowEIdx = max(enumerate(piBondEnergyWeights), key=lambda x: x[1])[0]
        finalHash["piBond"] = weightedAvg(piBond, piBondEnergyWeights)
        finalHash["piEnergy"] = weightedAvg(piEnergy, piBondEnergyWeights)
        finalHash["antiPiBond"] = weightedAvg(piAntiBond, piBondEnergyWeights)
        finalHash["antiPiEnergy"] = weightedAvg(antiPiEnergy, piBondEnergyWeights)

        finalHash["piBond_lowE"] = piBond[piBondlowEIdx]
        finalHash["piEnergy_lowE"] = piEnergy[piBondlowEIdx]
        finalHash["antiPiBond_lowE"] = piAntiBond[piBondlowEIdx]
        finalHash["antiPiEnergy_lowE"] = antiPiEnergy[piBondlowEIdx]
    elif charge == 0 and len(piBond) == 0:
        finalHash["piBond"] = "Nan"
        finalHash["piEnergy"] ="Nan"
        finalHash["antiPiBond"] ="Nan"
        finalHash["antiPiEnergy"] = "Nan"  
    return finalHash
def getAlkeneCHELPGInfo(logList , CC , energyStr , logNameMAST , smiles , logEnergyStr , chelpgStr , endIdx):
    weightsDF = getBoltzmannWeightsGauss(logList , 298 , energyStr , logEnergyStr)
    Cmin_CHELPG = []
    Cmax_CHELPG = []
    CHELPG_delta = []
    CHELPG_mean = []
    Cmin_chelpgElectro = []
    Cmax_chelpgElectro = []
    electroDelta = []
    electroMean = []
    for idx , row in weightsDF.iterrows():
        logFile = row["logID"]
        logPath = f"{logFile}.log"
        conformer  = next((f for f in logList if logPath in f.name), None)
        chelpgHash = extractCHELPGCharges(conformer , chelpgStr , endIdx)
        alkeneHash = {key: chelpgHash[key] for key in chelpgHash if key == CC[0] or key == CC[1]}
        c1CHELPG = float(alkeneHash[CC[0]][1])
        c2CHELPG = float(alkeneHash[CC[1]][1])
        if c2CHELPG >= c1CHELPG:
            Cmax_CHELPG.append(c2CHELPG)
            Cmin_CHELPG.append(c1CHELPG)
            Cmin_chelpgElectro.append(float(alkeneHash[CC[0]][2]))
            Cmax_chelpgElectro.append(float(alkeneHash[CC[1]][2]))
        else:
            Cmax_CHELPG.append(c1CHELPG)
            Cmin_CHELPG.append(c2CHELPG)
            Cmin_chelpgElectro.append(float(alkeneHash[CC[1]][2]))
            Cmax_chelpgElectro.append(float(alkeneHash[CC[0]][2]))

        CHELPG_delta.append(abs(c1CHELPG-c2CHELPG))
        CHELPG_mean.append(float(np.mean([c1CHELPG,c2CHELPG])))
        electroDelta.append(abs(float(alkeneHash[CC[1]][2])-float(alkeneHash[CC[0]][2])))
        electroMean.append(float(np.mean([float(alkeneHash[CC[0]][2]), float(alkeneHash[CC[1]][2])])))
    finalHash = {}
    finalHash["ID"] = logNameMAST
    finalHash["SMILES"] = smiles
    finalHash["CHELPG_mxAlk"] =weightedAvg(Cmin_CHELPG , list(weightsDF["boltzWeights"]))
    finalHash["CHELPG_mnAlk"] =weightedAvg(Cmax_CHELPG , list(weightsDF["boltzWeights"]))
    finalHash["CHELPG_Mean"] =weightedAvg(CHELPG_mean , list(weightsDF["boltzWeights"]))
    finalHash["CHELPG_Delta"] =weightedAvg(CHELPG_delta , list(weightsDF["boltzWeights"]))
    finalHash["CHELPG_mxElectro"] =weightedAvg(Cmin_chelpgElectro , list(weightsDF["boltzWeights"]))
    finalHash["CHELPG_mnElectro"] =weightedAvg(Cmax_chelpgElectro , list(weightsDF["boltzWeights"]))
    finalHash["CHELPG_ElectroMean"] =weightedAvg(electroMean , list(weightsDF["boltzWeights"]))
    finalHash["CHELPG_ElectroDelta"] =weightedAvg(electroDelta , list(weightsDF["boltzWeights"]))
    lowEidx = max(enumerate(list(weightsDF["boltzWeights"])), key=lambda x: x[1])[0]
    Cmin_CHELPG_lowE = Cmin_CHELPG[lowEidx]
    Cmax_CHELPG_lowE = Cmax_CHELPG[lowEidx]
    Cmin_chelpgElectro_lowE = Cmin_chelpgElectro[lowEidx]
    Cmax_chelpgElectro_lowE = Cmax_chelpgElectro[lowEidx]
    finalHash["CHELPG_mxAlk_lowE"] = Cmax_CHELPG_lowE
    finalHash["CHELPG_mnAlk_lowE"] =Cmin_CHELPG_lowE
    finalHash["CHELPG_Mean_lowE"] = np.mean([Cmin_CHELPG_lowE, Cmax_CHELPG_lowE])
    finalHash["CHELPG_Delta_lowE"] =abs(Cmin_CHELPG_lowE - Cmax_CHELPG_lowE)
    finalHash["CHELPG_mxElectro_lowE"] = Cmax_chelpgElectro_lowE
    finalHash["CHELPG_mnElectro_lowE"] = Cmin_chelpgElectro_lowE
    finalHash["CHELPG_ElectroMean_lowE"] = np.mean([Cmin_chelpgElectro_lowE, Cmax_chelpgElectro_lowE])
    finalHash["CHELPG_ElectroDelta_lowE"] =abs(Cmin_chelpgElectro_lowE - Cmax_chelpgElectro_lowE)
    return finalHash
def getAlkeneHirshInfo(logList , CC , energyStr , logNameMAST , smiles , logEnergyStr , hirshStr , endIdx):
    weightsDF = getBoltzmannWeightsGauss(logList , 298 , energyStr , logEnergyStr)
    Cmin_Hirsh = []
    Cmax_Hirsh = []
    Hirsh_delta = []
    Hirsh_mean = []
    Cmin_SpinUp = []    
    Cmin_SpinDown = []
    Cmax_SpinUp = []
    Cmax_SpinDown = []
    for idx , row in weightsDF.iterrows():
        logFile = row["logID"]
        logPath = f"{logFile}.log"
        conformer  = next((f for f in logList if logPath in f.name), None)
        hirshHash = extractHirshfeldCharges(conformer , hirshStr ,endIdx)
        alkeneHash = {key: hirshHash[key] for key in hirshHash if key == CC[0] or key == CC[1]}
        c1Hirsh = float(alkeneHash[CC[0]][3])
        c2Hirsh = float(alkeneHash[CC[1]][3])
        if c2Hirsh >= c1Hirsh:
            Cmax_Hirsh.append(c2Hirsh)
            Cmin_Hirsh.append(c1Hirsh)
            Cmin_SpinUp.append(float(alkeneHash[CC[0]][1]))
            Cmin_SpinDown.append(float(alkeneHash[CC[0]][2]))
            Cmax_SpinUp.append(float(alkeneHash[CC[1]][1]))
            Cmax_SpinDown.append(float(alkeneHash[CC[1]][2]))
        else:
            Cmax_Hirsh.append(c1Hirsh)
            Cmin_Hirsh.append(c2Hirsh)
            Cmin_SpinUp.append(float(alkeneHash[CC[0]][1]))
            Cmin_SpinDown.append(float(alkeneHash[CC[0]][2]))
            Cmax_SpinUp.append(float(alkeneHash[CC[1]][1]))
            Cmax_SpinDown.append(float(alkeneHash[CC[1]][2]))

        Hirsh_delta.append(abs(c1Hirsh-c2Hirsh))
        Hirsh_mean.append(float(np.mean([c1Hirsh,c2Hirsh])))
    finalHash = {}
    finalHash["ID"] = logNameMAST
    finalHash["SMILES"] = smiles
    finalHash["Hirsch_mxAlk"] =weightedAvg(Cmax_Hirsh , list(weightsDF["boltzWeights"]))
    finalHash["Hirsch_mnAlk"] =weightedAvg(Cmin_Hirsh , list(weightsDF["boltzWeights"]))
    finalHash["Hirsch_Mean"] =weightedAvg(Hirsh_mean , list(weightsDF["boltzWeights"]))
    finalHash["Hirsch_Delta"] =weightedAvg(Hirsh_delta , list(weightsDF["boltzWeights"]))
    finalHash["Hirsch_mxSpinUp"] =weightedAvg(Cmax_SpinUp , list(weightsDF["boltzWeights"]))
    finalHash["Hirsch_mnSpinUp"] =weightedAvg(Cmin_SpinUp , list(weightsDF["boltzWeights"]))
    finalHash["Hirsch_mxSpinDown"] =weightedAvg(Cmax_SpinDown , list(weightsDF["boltzWeights"]))
    finalHash["Hirsch_mnSpinDown"] =weightedAvg(Cmin_SpinDown , list(weightsDF["boltzWeights"]))
    lowEidx = max(enumerate(list(weightsDF["boltzWeights"])), key=lambda x: x[1])[0]
    Cmin_Hirsh_lowE = Cmin_Hirsh[lowEidx]
    Cmax_Hirsh_lowE = Cmax_Hirsh[lowEidx]
    Cmin_SpinUp_lowE = Cmin_SpinUp[lowEidx]
    Cmax_SpinUp_lowE = Cmax_SpinUp[lowEidx]
    Cmin_SpinDown_lowE = Cmin_SpinDown[lowEidx]
    Cmax_SpinDown_lowE = Cmax_SpinDown[lowEidx]
    finalHash["Hirsch_mxAlk_lowE"] = Cmax_Hirsh_lowE
    finalHash["Hirsch_mnAlk_lowE"] =Cmin_Hirsh_lowE
    finalHash["Hirsch_Mean_lowE"] = np.mean([Cmin_Hirsh_lowE, Cmax_Hirsh_lowE])
    finalHash["Hirsch_Delta_lowE"] =abs(Cmin_Hirsh_lowE - Cmax_Hirsh_lowE)
    finalHash["Hirsch_mxSpinUp_lowE"] = Cmax_SpinUp_lowE
    finalHash["Hirsch_mnSpinUp_lowE"] = Cmin_SpinUp_lowE
    finalHash["Hirsch_mxSpinDown_lowE"] = Cmax_SpinDown_lowE
    finalHash["Hirsch_mnSpinDown_lowE"] = Cmin_SpinDown_lowE
    return finalHash
def getAlkeneFullInfo(logList , CC , energyStr , logNameMAST , smiles , logEnergyStr , endIdx , fullStr , charge ):
    weightsDF = getBoltzmannWeightsGauss(logList , 298 , energyStr , logEnergyStr)
    mullikenCmx = []
    mullikenCmn = []
    mullikenDelta = []
    mullikenMean = []
    cMxPx = []
    cMnPx = []
    cMxPy = []
    cMnPy = []
    cMxPz = []
    cMnPz = []
    cMxPxE = []
    cMnPxE = []
    cMxPyE = []
    cMnPyE = []
    cMxPzE = []
    cMnPzE = []
    for idx , row in weightsDF.iterrows():
        logFile = row["logID"]
        logPath = f"{logFile}.log"
        conformer  = next((f for f in logList if logPath in f.name), None)
        print(conformer)
        fullHash = extractFullCharges(conformer , fullStr , endIdx , charge  )
        alkeneHash = {key: fullHash[key] for key in fullHash if key == CC[0] or key == CC[1]}
        c1Mulliken = float(alkeneHash[CC[0]]["mulliken"])
        c2Mulliken = float(alkeneHash[CC[1]]["mulliken"])
        if c2Mulliken >= c1Mulliken:
            mullikenCmx.append(c2Mulliken)
            mullikenCmn.append(c1Mulliken)
            if charge == 0:
                cMnP = alkeneHash[CC[0]]["twoOrbitals"]
                cMnPIdx = list(cMnP.keys())
                cMnPx.append(cMnP[cMnPIdx[0]][0])
                cMnPy.append(cMnP[cMnPIdx[1]][0])
                cMnPz.append(cMnP[cMnPIdx[2]][0])
                cMnPxE.append(cMnP[cMnPIdx[0]][1])
                cMnPyE.append(cMnP[cMnPIdx[1]][1])
                cMnPzE.append(cMnP[cMnPIdx[2]][1])
                CmxP = alkeneHash[CC[1]]["twoOrbitals"]
                cMPIdx = list(CmxP.keys())
                cMxPx.append(CmxP[cMPIdx[0]][0])
                cMxPy.append(CmxP[cMPIdx[1]][0])
                cMxPz.append(CmxP[cMPIdx[2]][0])
                cMxPxE.append(CmxP[cMPIdx[0]][1])
                cMxPyE.append(CmxP[cMPIdx[1]][1])
                cMxPzE.append(CmxP[cMPIdx[2]][1])
        else:
            mullikenCmx.append(c1Mulliken)
            mullikenCmn.append(c2Mulliken)
            if charge == 0:
                cMnP = alkeneHash[CC[1]]["twoOrbitals"]
                cMnPIdx = list(cMnP.keys())
                cMnPx.append(cMnP[cMnPIdx[0]][0])
                cMnPy.append(cMnP[cMnPIdx[1]][0])
                cMnPz.append(cMnP[cMnPIdx[2]][0])
                cMnPxE.append(cMnP[cMnPIdx[0]][1])
                cMnPyE.append(cMnP[cMnPIdx[1]][1])
                cMnPzE.append(cMnP[cMnPIdx[2]][1])
                CmxP = alkeneHash[CC[0]]["twoOrbitals"]
                cMPIdx = list(CmxP.keys())
                cMxPx.append(CmxP[cMPIdx[0]][0])
                cMxPy.append(CmxP[cMPIdx[1]][0])
                cMxPz.append(CmxP[cMPIdx[2]][0])
                cMxPxE.append(CmxP[cMPIdx[0]][1])
                cMxPyE.append(CmxP[cMPIdx[1]][1])
                cMxPzE.append(CmxP[cMPIdx[2]][1])
        mullikenDelta.append(abs(c1Mulliken-c2Mulliken))
        mullikenMean.append(float(np.mean([c1Mulliken,c2Mulliken])))
    finalHash = {}
    finalHash["ID"] = logNameMAST
    finalHash["SMILES"] = smiles
    finalHash["Mulliken_mxAlk"] =weightedAvg(mullikenCmx , list(weightsDF["boltzWeights"]))
    finalHash["Mulliken_mnAlk"] =weightedAvg(mullikenCmn , list(weightsDF["boltzWeights"]))
    finalHash["Mulliken_Mean"] =weightedAvg(mullikenMean , list(weightsDF["boltzWeights"]))
    finalHash["Mulliken_Delta"] =weightedAvg(mullikenDelta , list(weightsDF["boltzWeights"]))
    if charge == 0:
        finalHash["Mulliken_mxPx"] =weightedAvg(cMxPx , list(weightsDF["boltzWeights"]))
        finalHash["Mulliken_mnPx"] =weightedAvg(cMnPx , list(weightsDF["boltzWeights"]))
        finalHash["Mulliken_mxPy"] =weightedAvg(cMxPy , list(weightsDF["boltzWeights"]))
        finalHash["Mulliken_mnPy"] =weightedAvg(cMnPy , list(weightsDF["boltzWeights"]))
        finalHash["Mulliken_mxPz"] =weightedAvg(cMxPz , list(weightsDF["boltzWeights"]))
        finalHash["Mulliken_mnPz"] =weightedAvg(cMnPz , list(weightsDF["boltzWeights"]))
        finalHash["Mulliken_mxPxE"] =weightedAvg(cMxPxE , list(weightsDF["boltzWeights"]))
        finalHash["Mulliken_mnPxE"] =weightedAvg(cMnPxE , list(weightsDF["boltzWeights"]))
        finalHash["Mulliken_mxPyE"] =weightedAvg(cMxPyE , list(weightsDF["boltzWeights"]))
        finalHash["Mulliken_mnPyE"] =weightedAvg(cMnPyE , list(weightsDF["boltzWeights"]))
        finalHash["Mulliken_mxPzE"] =weightedAvg(cMxPzE , list(weightsDF["boltzWeights"]))
        finalHash["Mulliken_mnPzE"] =weightedAvg(cMnPzE , list(weightsDF["boltzWeights"]))
    lowEidx = max(enumerate(list(weightsDF["boltzWeights"])), key=lambda x: x[1])[0]
    Cmin_Mulliken_lowE = mullikenCmn[lowEidx]
    Cmax_Mulliken_lowE = mullikenCmx[lowEidx]
    finalHash["Mulliken_mxAlk_lowE"] = Cmax_Mulliken_lowE
    finalHash["Mulliken_mnAlk_lowE"] =Cmin_Mulliken_lowE
    finalHash["Mulliken_Mean_lowE"] = np.mean([Cmin_Mulliken_lowE, Cmax_Mulliken_lowE])
    finalHash["Mulliken_Delta_lowE"] =abs(Cmin_Mulliken_lowE - Cmax_Mulliken_lowE)
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
            minWeightID = max(enumerate(weights), key=lambda x: x[1])[0]
            C13_Cmx_min = C13_C2[minWeightID]
            C13_Cmn_min = C13_C1[minWeightID]
            C13_delta_min = abs(C13_Cmx_min - C13_Cmn_min)
            C13_Mean_min = np.mean([C13_Cmx_min,C13_Cmn_min])

            C13_Cmx = ( (C13_C2 * weights).sum()    / weights.sum()  )
            C13_Cmn = ( (C13_C1 * weights).sum()    / weights.sum()  )
            C13_delta = abs(C13_Cmx - C13_Cmn)
            C13_Mean = np.mean([C13_Cmx,C13_Cmn])
            c13Hash = {"C13_Cmx" : C13_Cmx , "C13_Cmn" : C13_Cmn , "C13_delta" : C13_delta , "C13_Mean" : C13_Mean , 
                       "C13_Cmx_lowE" : C13_Cmx_min , "C13_Cmn_lowE" : C13_Cmn_min , "C13_delta_lowE" : C13_delta_min , "C13_Mean_lowE" : C13_Mean_min}
            hashList.append(c13Hash)
        else:
            C1 = cc[1] + 1
            C2 = cc[0] + 1 

        if "NBO7" in featureList:
            nboNeutralStr = featureHash["NBO7"][0]
            #print(C1,C2)
            neutralHash = getAlkeneNBOInfo(conformerFiles , C1 , C2, "electronic", id , smilesStr , nboNeutralStr , 0 , logEnergyStr)
            hashList.append(neutralHash)
        if "CHELPG" in featureList:
            chelpgStr = featureHash["CHELPG"][0]
            chelpgHash = getAlkeneCHELPGInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , chelpgStr , 1)
            hashList.append(chelpgHash)
        if "Hirshfeld" in featureList:
            hirshStr = featureHash["Hirshfeld"][0]
            hirshHash = getAlkeneHirshInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , hirshStr , 1) 
            hashList.append(hirshHash)
        if "FullCharges" in featureList:
            fullStr = featureHash["FullCharges"][0]
            fullHash = getAlkeneFullInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , 1 , fullStr , 0) #(logList , CC , energyStr , logNameMAST , smiles , logEnergyStr , endIdx , fullStr)
            hashList.append(fullHash)
        if "fukuiParameters" in featureList:
            electroCationStr = featureHash["fukuiParameters"][0]
            electroAnionstr = featureHash["fukuiParameters"][1]
            electroNeutralStr = featureHash["fukuiParameters"][2]
            if electroNeutralStr == "NBO7":
                nboNeutralStr = featureHash[f"{electroNeutralStr}"][0]
                neutralHash = getAlkeneNBOInfo(conformerFiles ,  C1 , C2, "electronic", id , smilesStr , nboNeutralStr , 0, logEnergyStr)
                cationHash = getAlkeneNBOInfo(conformerFiles ,  C1 , C2 , "electronic", id , smilesStr , electroCationStr , 1 , logEnergyStr)
                anionHash = getAlkeneNBOInfo(conformerFiles ,  C1 , C2,"electronic", id , smilesStr, electroAnionstr , -1, logEnergyStr)
                q0_mx = neutralHash["NBO_mxAlk"]
                q1_mx = cationHash["NBO_mxAlk"]
                q1__mx = anionHash["NBO_mxAlk"]
                q0_mn = neutralHash["NBO_mnAlk"]
                q1_mn = cationHash["NBO_mnAlk"]
                q1__mn = anionHash["NBO_mnAlk"]
            elif electroNeutralStr == "CHELPG":
                chelpgStr = featureHash[f"{electroNeutralStr}"][0]
                neutralHash = getAlkeneCHELPGInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , chelpgStr , 1)
                cationHash = getAlkeneCHELPGInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , electroCationStr , 2)
                anionHash = getAlkeneCHELPGInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , electroAnionstr , 3)
                q0_mx = neutralHash["CHELPG_mxAlk"]
                q1_mx = cationHash["CHELPG_mxAlk"]
                q1__mx = anionHash["CHELPG_mxAlk"]
                q0_mn = neutralHash["CHELPG_mnAlk"]
                q1_mn = cationHash["CHELPG_mnAlk"]
                q1__mn = anionHash["CHELPG_mnAlk"]
            elif electroNeutralStr == "Hirshfeld":
                hirshStr = featureHash[f"{electroNeutralStr}"][0]
                neutralHash = getAlkeneHirshInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , hirshStr , 1)
                cationHash = getAlkeneHirshInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , electroCationStr , 2)
                anionHash = getAlkeneHirshInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , electroAnionstr , 3)
                q0_mx = neutralHash["Hirsch_mxAlk"]
                q1_mx = cationHash["Hirsch_mxAlk"]
                q1__mx = anionHash["Hirsch_mxAlk"]
                q0_mn = neutralHash["Hirsch_mnAlk"]
                q1_mn = cationHash["Hirsch_mnAlk"]
                q1__mn = anionHash["Hirsch_mnAlk"]
            elif electroNeutralStr == "FullCharges":
                fullStr = featureHash[f"{electroNeutralStr}"][0]
                neutralHash = getAlkeneFullInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , 1 , fullStr , 0)
                cationHash = getAlkeneFullInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , 2 , electroCationStr , 1)
                anionHash = getAlkeneFullInfo(conformerFiles , [C1,C2] , "electronic", id , smilesStr , logEnergyStr , 3 , electroAnionstr , -1)
                q0_mx = neutralHash["Mulliken_mxAlk"]
                q1_mx = cationHash["Mulliken_mxAlk"]
                q1__mx = anionHash["Mulliken_mxAlk"]
                q0_mn = neutralHash["Mulliken_mnAlk"]
                q1_mn = cationHash["Mulliken_mnAlk"]
                q1__mn = anionHash["Mulliken_mnAlk"]
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

            radList = [2.0,2.5 , 3.0 , 3.5]
            weights = boltzmannDF["boltzWeights"].to_numpy()
            minWeightIdx = max(enumerate(weights), key=lambda x: x[1])[0]
            weight_sum = weights.sum()

            # Storage
            vBurStorage = {C : {r : [] for r in radList} for C in ("Cmn" , "Cmx")}

            for name in boltzmannDF["logID"]:
                fileStr = f"{name}.log"
                conformer = next((f for f in conformerFiles if fileStr in f.name),None)
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
                    vburC1 = BuriedVolume(elements,coordinates,int(C1 - 1),include_hs=True,radius=rad).fraction_buried_volume
                    vburC2 = BuriedVolume(elements,coordinates,int(C2 - 1),include_hs=True,radius=rad).fraction_buried_volume
                    
                    vBurStorage["Cmn"][rad].append(vburC1)
                    vBurStorage["Cmx"][rad].append(vburC2)

            vBurAvg = {atom : {r : (np.asarray(vBurStorage[atom][r])*weights).sum()/weight_sum for r in radList} for atom in ("Cmn","Cmx")}
            vBurAvg_lowE = {atom : {r : vBurStorage[atom][r][minWeightIdx] for r in radList} for atom in ("Cmn","Cmx")}

            BurVolHash = {}
            for rad in radList:
                vBurList = []
                vBurlowE = []
                for atom in ("Cmn","Cmx"):
                    atom1 = vBurAvg[atom][rad]
                    atom1_lowE = vBurAvg_lowE[atom][rad]
                    vBurlowE.append(atom1_lowE)
                    vBurList.append(atom1)
                    BurVolHash[f"{rad}_Ang_Vburr_{atom}"] = atom1
                    BurVolHash[f"{rad}_Ang_Vburr_{atom}_lowE"] = atom1_lowE
                    if len(vBurList) == 2:
                        delta = np.abs(vBurList[0] - vBurList[1])
                        mean = np.mean(vBurList)
                        BurVolHash[f"{rad}_Ang_Vburr_mean"] = delta
                        BurVolHash[f"{rad}_Ang_Vburr_delta"] = mean
                        deltaLowE = np.abs(vBurlowE[0] - vBurlowE[1])
                        meanLowE = np.mean(vBurlowE)
                        BurVolHash[f"{rad}_Ang_Vburr_mean_lowE"] = deltaLowE
                        BurVolHash[f"{rad}_Ang_Vburr_delta_lowE"] = meanLowE
            hashList.append(BurVolHash)
        if "firstContactVbur" in featureList:
            radList = [2.0,2.5,3.0,3.5,4.0]
            vBurHash = {}
            weights = boltzmannDF["boltzWeights"]
            firstContact_lowEIdx = max(enumerate(weights), key=lambda x: x[1])[0]
            for i in range(len(cc)):
                alkeneCarbon = cc[i]
                otherCarbon = cc[i-1]
                #print(otherCarbon)
                atom = molec.GetAtomWithIdx(alkeneCarbon)
                neighbors = atom.GetNeighbors()
                neighborLst = [nbr.GetIdx() for nbr in neighbors if nbr.GetIdx() != otherCarbon]
                #print(neighborLst)
                vBurHash[f"C_{alkeneCarbon}"] = {neighbor : {rad :[] for rad in radList}for neighbor in neighborLst}
            idxCount = 0
            for name in list(boltzmannDF["logID"]):
                fileStr = f"{name}.log"
                conformer  = next((f for f in conformerFiles if fileStr in f.name), None)
                coordHash = getAtomCoordsRobust(str(conformer) , "GINC-COMPUTE" , 5  , 1)
                elements = []
                coordinates = []
                for _, coords in coordHash.items():
                    elements.append(str(coords[0]))
                    coordinates.append(np.array(coords[2:5]))
                if idxCount == firstContact_lowEIdx:
                    c1Vbur = BuriedVolume(elements,coordinates,int(cc[0]),include_hs=True,radius=2.0).fraction_buried_volume
                    c2Vbur = BuriedVolume(elements,coordinates,int(cc[1]),include_hs=True,radius=2.0).fraction_buried_volume
                    if c1Vbur >= c2Vbur:
                        orientationHash = {f"C_{cc[0]}" : "C_1" , f"C_{cc[1]}" : "C_2"}
                    else:
                        orientationHash = {f"C_{cc[0]}" : "C_2" , f"C_{cc[1]}" : "C_1"}
                    c1 = next((k for k, v in orientationHash.items() if v == "C_1"), None)
                    c2 = next((k for k, v in orientationHash.items() if v == "C_2"), None)
                    c1Contacts = list(vBurHash[c1].keys())
                    c2Contacts = list(vBurHash[c2].keys())
                    length1 = np.linalg.norm(coordinates[c1Contacts[0]] - coordinates[c2Contacts[0]])
                    length2 = np.linalg.norm(coordinates[c1Contacts[0]] - coordinates[c2Contacts[1]])
                    if len(set(c1Contacts + c2Contacts)) == 3:
                        #some tricyclic things going on
                        intersection = list(set(c1Contacts) & set(c2Contacts))[0]
                        bondHash = {"face1" : [intersection] , "face2" : [contact for contact in (c1Contacts + c2Contacts) if contact != intersection]}
                    elif len(set(c1Contacts + c2Contacts)) == 2:
                        intersection = list(set(c1Contacts) & set(c2Contacts))[0]
                        intersection2 = list(set(c1Contacts) & set(c2Contacts))[1]
                        bondHash = {"face1" : [intersection] , "face2" : [intersection2]}
                    else:

                        if length1 < length2:
                            #1 and 1 orientations are good 
                            bondHash = {"face1" : [c1Contacts[0] , c2Contacts[0]] , "face2" : [c1Contacts[1] , c2Contacts[1]]}
                        else:
                            #1 and 2 orientations need to be switched
                            bondHash = {"face1" : [c1Contacts[0] , c2Contacts[1]] , "face2" : [c1Contacts[1] , c2Contacts[0]]}
                idxCount +=1
                for atom in vBurHash:
                    for neighbor in vBurHash[atom]:
                        for rad in vBurHash[atom][neighbor]:
                            vBurHash[atom][neighbor][rad].append(BuriedVolume(elements,coordinates,int(neighbor),include_hs=True,radius=rad).fraction_buried_volume)
            vBurContactsAvg = {atom : {contact : {rad : ((np.asarray(vBurHash[atom][contact][rad])*weights).sum() / weights.sum()) for rad in vBurHash[atom][contact]} for contact in vBurHash[atom]} for atom in vBurHash}
            vBurContactslowE = {atom : {contact : {rad : vBurHash[atom][contact][rad][firstContact_lowEIdx] for rad in vBurHash[atom][contact]} for contact in vBurHash[atom]} for atom in vBurHash}
            vBurContacts = {}
            #print(bondHash)
            for atom in vBurHash:
                orientation = orientationHash[atom]
                #print(atom)
                for contact in vBurContactsAvg[atom]:
                    #print(contact)
                    if contact in bondHash["face1"]:
                        face = "1"
                    else:
                        face = "2"
                    for rad in radList:
                        vBurContacts[f"{orientation}_{face}_{rad}"] = vBurContactsAvg[atom][contact][rad]
                        vBurContacts[f"{orientation}_{face}_{rad}_lowE"] = vBurContactslowE[atom][contact][rad]
            #print(vBurContacts)
            for rad in radList:
                avgFace1 = (vBurContacts[f"C_1_1_{rad}"] + vBurContacts[f"C_2_1_{rad}"]) /2 
                avgFace2 = (vBurContacts[f"C_1_2_{rad}"] + vBurContacts[f"C_2_2_{rad}"]) /2 
                avg1 = (vBurContacts[f"C_1_1_{rad}"] + vBurContacts[f"C_1_2_{rad}"]) /2 
                avg2 = (vBurContacts[f"C_2_1_{rad}"] + vBurContacts[f"C_2_2_{rad}"]) /2 
                vBurContacts[f"delta_Faces_{rad}"] = np.abs(avgFace1-avgFace2)
                vBurContacts[f"delta_Ends_{rad}"] = np.abs(avg1-avg2)


            hashList.append(vBurContacts)      
        if "orangeSlices" in featureList:
            radList = [2.0,2.5,3.0]

            orangeData = {atom: {quad: {r: [] for r in radList} for quad in range(4)} for atom in ("C1", "C2")}
            g = Graph()
            for bond in molec.GetBonds():
                start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                g.add_edge(start, end)
            CminNeighbors = list(g.neighbors(int(C1-1)))
            CminNeighbors.remove(int(C2-1))
            CmaxNeighbors = list(g.neighbors(int(C2-1)))
            CmaxNeighbors.remove(int(C1-1))

            for name in list(boltzmannDF["logID"]):
                CmaxContacts = []
                CminContacts = []
                fileStr = f"{name}.log"
                conformer  = next((f for f in conformerFiles if fileStr in f.name), None)
                coordHash = getAtomCoordsRobust(str(conformer) , "GINC-COMPUTE" , 5  , 1)

                atoms = []
                idxList = []
                symbols = []
                for atomIdx , coords in coordHash.items():
                    idxList.append(atomIdx)
                    atoms.append(coords[2:5])
                    symbols.append(coords[0])
                    atomIdx +=1
                    if atomIdx == C1:
                        c1_coords = np.array(coords[2:5])
                        c1Idx = atomIdx
                    elif atomIdx == C2:
                        c2_coords = coords[2:5]
                        c2Idx = atomIdx
                    if (atomIdx-1) in CminNeighbors:
                        crds = np.array(coords[2:5])
                        CminContacts.append(crds) 
                    if (atomIdx-1) in CmaxNeighbors:
                        crds = np.array(coords[2:5])
                        CmaxContacts.append(crds) 
                C1Hash = {"0" : c1_coords , "1" : CminContacts , "idx" : c1Idx}
                C2Hash = {"0" : c2_coords , "1" : CmaxContacts , "idx" : c2Idx}
                for rad in radList:
                    alkeneOranges = alkeneSlicedOranges(C1Hash , C2Hash , rad)
                    alkeneOranges.getAtoms(atoms , idxList , symbols ,{"Nan" : "Nan"} , False)
                    orangeSlices = alkeneOranges.slicedOranges(True , False , False ,"alkene")

                    for center in ("C1", "C2"):
                        for quad in range(4):
                            orangeData[center][quad][rad].append(
                                orangeSlices[f"{center}_quad_{quad}"]
                            )
            weights = boltzmannDF["boltzWeights"]
            orangeSlices_lowEIdx = max(enumerate(weights), key=lambda x: x[1])[0]
            weight_sum = weights.sum()

            vburAvg = {
                center: {
                    quad: {
                        r: (
                            np.asarray(orangeData[center][quad][r]) * weights
                        ).sum() / weight_sum
                        for r in radList
                    }
                    for quad in range(4)
                }
                for center in ("C1", "C2")
            }
            vburLowE = {
                center: {
                    quad: {
                        r: orangeData[center][quad][r][orangeSlices_lowEIdx]
                        for r in radList
                    }
                    for quad in range(4)
                }
                for center in ("C1", "C2")
            }
            vburOrangeSlices = {}

            for center in ("C1", "C2"):
                for quad in range(4):
                    for rad in radList:

                        vburOrangeSlices[
                            f"vbur{center}_slice{quad}_{rad}"
                        ] = vburAvg[center][quad][rad]

                        vburOrangeSlices[
                            f"vbur{center}_slice{quad}_{rad}_lowE"
                        ] = vburLowE[center][quad][rad]
            hashList.append(vburOrangeSlices)
        if "%VburSemiCylinders" in featureList:
            radList = [2.5 , 3.0 , 3.5]
            segmentList = ["maxSemiPi" , "minSemiPi" , "maxSemiOrth" , "minSemiOrth" , "maxCap" , "minCap" , "CburrTot"]
            cappedSemis = {segment : {r : [] for r in radList} for segment in segmentList}
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
                    mainCylinder = alkeneSemiCylinders(C1Hash , C2Hash , rad , 0.15)
                    mainCylinder.getAtoms(coordHash ,{"Nan" : "Nan"} , False)
                    segmentHash = mainCylinder.getBurriedVolume(True , False)
                    #Cburr = mainCylinder.getBurriedVolume(False , False)
                    for segment in segmentList:
                        cappedSemis[segment][rad].append(segmentHash[segment])

            weights = boltzmannDF["boltzWeights"]
            SemiCylinder_lowEIdx = max(enumerate(weights), key=lambda x: x[1])[0]
            weight_sum = weights.sum()
            #print(cappedSemis)
            semiCylinders = {segment : {rad : (np.asarray(cappedSemis[segment][rad]) * weights).sum() for rad in radList } for segment in segmentList}
            semiCylinders_lowE = {segment : {rad : cappedSemis[segment][rad][SemiCylinder_lowEIdx] for rad in radList } for segment in segmentList}

            Vburr_SemiCylinders ={}
            for rad in radList:
                piList = []
                piList_lowE = []
                orthList = []
                orthList_lowE = []
                capList = []
                capList_lowE = []
                for segment in segmentList:
                    Vburr_SemiCylinders[f"Vbur_{segment}_{rad}"] = semiCylinders[segment][rad]
                    Vburr_SemiCylinders[f"Vbur_{segment}_{rad}_lowE"] = semiCylinders_lowE[segment][rad]
                    if "Pi" in segment:
                        piList.append(semiCylinders[segment][rad])
                        piList_lowE.append(semiCylinders_lowE[segment][rad])
                        if len(piList) == 2 and len(piList_lowE) == 2:
                            Vburr_SemiCylinders[f"Vbur_deltaPi_{rad}"] = np.abs(piList[0] - piList[1])
                            Vburr_SemiCylinders[f"Vbur_deltaPi_{rad}_lowE"] = np.abs(piList_lowE[0] - piList_lowE[1])
                    elif "Orth" in segment:
                        orthList.append(semiCylinders[segment][rad])
                        orthList_lowE.append(semiCylinders_lowE[segment][rad])
                        if len(orthList) == 2 and len(orthList_lowE) == 2:
                            Vburr_SemiCylinders[f"Vbur_deltaOrth_{rad}"] = np.abs(orthList[0] - orthList[1])
                            Vburr_SemiCylinders[f"Vbur_deltaOrth_{rad}_lowE"] = np.abs(orthList_lowE[0] - orthList_lowE[1])
                    elif "Cap" in segment:
                        capList.append(semiCylinders[segment][rad])
                        capList_lowE.append(semiCylinders_lowE[segment][rad])
                        if len(capList) == 2 and len(capList_lowE) == 2:
                            Vburr_SemiCylinders[f"Vbur_deltaCap_{rad}"] = np.abs(capList[0] - capList[1])
                            Vburr_SemiCylinders[f"Vbur_deltaCap_{rad}_lowE"] = np.abs(capList_lowE[0] - capList_lowE[1])
                    else:
                        continue
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
            alkeneDist_lowE_idx = max(enumerate(weights), key = lambda x:x[1])[0]
            alkeneDist_lowE = distList[alkeneDist_lowE_idx]
            distMAST = ( (distList * weights).sum()    / weights.sum()  )  
            distHash = {"C1C2Dist" : distMAST , "C1C2Dist_lowE" : alkeneDist_lowE} 
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
        if "TSEI" in featureList:
            maxPeriod = featureHash["TSEI"][0]
            if maxPeriod != 1:
                periodRange = np.arange(1, maxPeriod)
            else:
                periodRange = [1]
            for period in periodRange:
                tselList = calcSingle(smilesStr , period , True)
                maxTopology = max(tselList)
                minTopology = min(tselList)
                deltaTopology = np.abs(maxTopology - minTopology)
                meanTopology = (maxTopology + minTopology) / 2 
                hashList.append({f"maxTopology_{period}" : maxTopology , 
                                 f"minTopology_{period}" : minTopology , 
                                 f"deltaTopology_{period}" : deltaTopology , 
                                 f"meanTopology_{period}" : meanTopology , })

        if "Sterimol" in featureList:
            radList = [2.0 , 2.5, 3.0,4.0,4.5]
            weights = boltzmannDF["boltzWeights"].to_numpy()
            weight_sum = weights.sum()
            sterimol_MinIdx = max(enumerate(weights), key=lambda x: x[1])[0]
            # Storage
            sterimolHash = {c : {rad : { ang : [] for ang in ["B1" , "B5" , "L"]} for rad in radList} for c in ["C1" , "C2"]}

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
                    cStr = 0
                    for i in range(len([C1 , C2])):
                        cStr +=1
                        cH = f"C" + str(cStr)
                        alkeneCarbon = [C1 , C2][i]
                        otherCarbon = [C1 , C2][i-1]
                        sterimol_values = Sterimol(elements, coordinates, int(alkeneCarbon), int(otherCarbon)) 
                        sterimol_values.bury(method="delete", sphere_radius=float(rad))
                        L = sterimol_values.L_value
                        B1 = sterimol_values.B_1_value
                        B5 = sterimol_values.B_5_value
                        sterimolHash[cH][rad]["B1"].append(B1)
                        sterimolHash[cH][rad]["L"].append(L)
                        sterimolHash[cH][rad]["B5"].append(B5)
            sterimolLowE = {c : {rad : { ang : sterimolHash[c][rad][ang][sterimol_MinIdx] for ang in ["B1" , "B5" , "L"]} for rad in radList} for c in ["C1" , "C2"]}
            sterimolBoltzAvg = {c : {rad : { ang : (np.asarray(sterimolHash[c][rad][ang])*weights).sum() for ang in ["B1" , "B5" , "L"]} for rad in radList} for c in ["C1" , "C2"]}
            sterHash = {}
            for c in ["C1" , "C2"]:
                for rad in radList:
                    for ang in ["B1" , "B5" , "L"]:
                        sterHash[f"{c}_sterimol_{rad}_{ang}"] = sterimolBoltzAvg[c][rad][ang]
                        sterHash[f"{c}_sterimol_{rad}_{ang}_lowE"] = sterimolLowE[c][rad][ang]

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
