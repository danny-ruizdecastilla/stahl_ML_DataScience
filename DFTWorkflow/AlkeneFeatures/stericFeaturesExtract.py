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
from dimensionalityReduction.reactivityFeatures import boxGen
from breadthFirstSearch.radialBasedCorrelation import getCC
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
#Extract Steric Features from Crest raw XTB coordinates 

class CrestChemical:
    def __init__(self, name, strucFile, energyFile):
        self.name = name
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
                    coordHash = {}
                    collectingCoordinates = True
                    atomCount = 0
                elif collectingCoordinates and is_float(line.strip()):
                    energyLevel = float(line.strip())
                elif collectingCoordinates and line.split("         ")[0].strip()[:1].isalpha():
                    lineOptions = line.strip().split("    ")
                    atom = str(lineOptions[0])
                    coords = [num.strip() for num in lineOptions if is_float(num.strip())]
                    coordHash[str(idx) + "," + atom] = coords
                    atomCount +=1
                    if atomCount == int(self.NumAtoms):
                        self.conformersHash[energyLevel] = coordHash 
                        collectingCoordinates = False
        return self.conformersHash


class Conformer:
    def __init__(self, parent, xyz, energy):
        if parent is None:
            raise ValueError("Conformer must have a parent CrestChemical")
        self.parent = parent
        self.xyz = xyz
        self.energy = energy