import os 
import sys
import glob
import pandas as pd 
import numpy as np
import chemdraw
import base64
from pathlib import Path
from morfeus import BuriedVolume
from rdkit import Chem
from networkx import Graph
from itertools import islice
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
from DFTWorkflow.ionComGenerator import locateinLog
from DFTWorkflow.AlkeneFeatures.alkeneStericDougnut import alkeneSemiCylinders
from dimensionalityReduction.reactivityFeatures import boxGen
from breadthFirstSearch.radialBasedCorrelation import getCC
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
#Extract Steric Features from Crest raw XTB coordinates 

class CrestChemical:
    def __init__(self, name, smiles , strucFile, energyFile):
        self.name = name
        self.smiles = smiles
        self.strucFile = strucFile
        with open(strucFile, 'r') as f:
            counts = 0
            first = f.readline().strip()
            for line in f:
                line = line.strip()
                if line == first:
                    counts +=1
        self.NumAtoms = first
        self.NumConfs = counts 
        self.energyFile = energyFile
        self.conformersHash = {}
    def getConformers(self):
        def is_float(s):
            try:
                float(s)
                return True
            except ValueError:
                return False
        with open(self.strucFile , 'r' ) as f:
            collectingCoordinates = False
            for idx , line in enumerate(f):
                if not collectingCoordinates and line.strip() == self.NumAtoms:   #new conformer!
                    coordHash = {"atoms" : [] , "coordinates" : []}
                    collectingCoordinates = True
                    atomCount = 0
                elif collectingCoordinates and is_float(line.strip()):
                    energyLevel = float(line.strip())
                elif collectingCoordinates and line.split("         ")[0].strip()[:1].isalpha():
                    lineOptions = line.strip().split("    ")
                    atom = str(lineOptions[0])
                    coords = [float(num.strip()) for num in lineOptions if is_float(num.strip())]
                    coordHash["atoms"].append(atom)
                    coordHash["coordinates"].append(np.array(coords))
                    atomCount +=1
                    if atomCount == int(self.NumAtoms):
                        self.conformersHash[energyLevel] = coordHash 
                        collectingCoordinates = False
        return self.conformersHash
    def getBoltzmannWeights(self , hartree , temperature ):
        if hartree:
            R = 1.987204e-3
            HARTREE_TO_KCAL = 627.5094740631
            energies = list(self.conformersHash.keys())
            kCalE = [e * HARTREE_TO_KCAL for e in energies ]
            minE = min(kCalE)
            rE_kCal = [e - minE for e in kCalE ]
            boltzmannFacts =[np.exp((-1 * e) / (R * temperature)) for e in rE_kCal ]
            self.boltzmannFactors = boltzmannFacts 
            self.relEnergies = rE_kCal
    def getBurriedVol(self , atoms , radList ): #atoms = { idx : ["C"  , "alkene"] , ix : ["C" , "alkene2"]}
        burriedVolumeHash = {}
        boltzSum = sum(self.boltzmannFactors)
        for rad in radList:
            vBurrList = []
            for atom , name in atoms.items():
                symbol = name[0]
                nameStr = name[1]
                saveStr = f"{rad}_Ang_Vburr_{symbol}{nameStr}"
                enum = 0
                totVburr = 0
                for energy , coordsHash in self.conformersHash.items():
                    boltzFact = self.boltzmannFactors[enum]
                    vBur = BuriedVolume(coordsHash["atoms"],coordsHash["coordinates"],atom,include_hs=True,radius=rad).fraction_buried_volume
                    enum +=1
                    weightedVbur = boltzFact * vBur
                    tutVburr += weightedVbur
                avgVbur = totVburr / boltzSum
                burriedVolumeHash[saveStr] = avgVbur
                vBurrList.append(avgVbur)
            meanVburr = np.mean(vBurrList)
            deltaVburr = abs(vBurrList[0] - vBurrList[1])
            burriedVolumeHash[f"{rad}_Ang_Vburr_mean"] = meanVburr
            burriedVolumeHash[f"{rad}_Ang_Vburr_delta"] = deltaVburr
        return burriedVolumeHash
    def getBurriedCylinder(self , alkeneIndex , radList , g ):
        burriedVolumeHash = {}
        boltzSum = sum(self.boltzmannFactors)
        for rad in radList:
            enum = 0
            for energy , coordsHash in self.conformersHash.items():
                boltzFact = self.boltzmannFactors[enum]
                coordHash = {}
                for i in range(self.NumAtoms):
                    atom = coordsHash["atoms"][i]
                    coordinates = coordsHash["coordinates"][i]
                    coordHash[i] = [atom , "0" ]
                    for coord in coordinates:
                        coordHash[i].append(coord)
                    if i in alkeneIndex:
                        contacts  =  list(g.neighbors(i))
                        
                mainCylinder = alkeneSemiCylinders(C1Hash , C2Hash , rad)
                mainCylinder.getAtoms(coordHash)
                maxSemi , minSemi = mainCylinder.getBurriedVolume(True)
                Cburr = mainCylinder.getBurriedVolume(False)
                enum +=1
                weightedVbur = boltzFact * vBur
                tutVburr += weightedVbur
            avgVbur = totVburr / boltzSum
            burriedVolumeHash[saveStr] = avgVbur

        return burriedVolumeHash