import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
from rdkit.Chem.PandasTools import LoadSDF
import ast
import networkx as nx
from networkx import Graph
import matplotlib.pyplot as plt
from itertools import combinations
from collections import deque 
import random
def minFirstSearch(atomList , g , rejections):
    graphSubstitutions = {}
    for j in range (len(atomList)): #Need to extend network to each carbon in the Alkene 
        #print(j)
        C1 = atomList[j] #leading Carbon 
        print("carbon:" , C1)
        atom = C1
        branching = {}
        atomConnections = [C1]
        while True:
            #print("atom" , atom)
            #print(list(branching.values()))
            contacts = list(g.neighbors(atom))
            eligibles = [contact for contact in contacts if contact not in rejections]
            #print("eligibles" , eligibles)
            branching[atom] = len(eligibles)
            if len(eligibles) == 0:
                remainings = list(branching.values())
                if all(value == 0 for value in remainings):
                    break
                else:
                    remainingAtoms = [key for key, value in branching.items() if value != 0]
                    nextAtom = min(remainingAtoms)
                    #atomConnections.append(nextAtom)
                    rejections.append(atom)
                    atom = nextAtom
            else:
                
                nextAtom = min(eligibles)
                atomConnections.append(nextAtom)
                print("nextAtom" , nextAtom)
                rejections.append(atom)
                atom = nextAtom
            graphSubstitutions[C1] = atomConnections
    return graphSubstitutions
def molecWeight(molec, index):
    ptable = Chem.GetPeriodicTable()
    atom = molec.GetAtomWithIdx(index)
    atomicMass = ptable.GetAtomicWeight(atom.GetAtomicNum())
    return atomicMass
def maxKeyOfMaxValue(dict_):
    maxVal = max(dict_.values())
    keys = [k for k, v in dict_.items() if v == maxVal]
    return max(keys)
def getMaxWeight(molec, idxList):
    weights = {}

    for idx in idxList:
        weight = molecWeight(molec, idx)
        weights[idx] = weight
    maxKey = maxKeyOfMaxValue(weights)
    return maxKey
def maxPathCompare(molec, g , atomList  , rejections):
    pathDict = {}
    for idx in atomList:
            atomicMass = molecWeight(molec , idx)
            pathDict[idx] = deque([ atomicMass , idx])
    termPaths = []
    while True:
        paths = list(pathDict.keys())

        paths = [path for path in paths if path not in termPaths]
        if len(paths) == 0:
            break
        else:
            currentWeights = []
            for atom in paths:
                currentPath = pathDict[atom]
                currentAtom = currentPath[-1]
                contacts = list(g.neighbors(currentAtom))
                eligibles = [contact for contact in contacts if contact not in rejections]
                if len(eligibles) == 0:
                    termPaths.append(atom)
                else:
                    nextAtom = getMaxWeight( molec, eligibles)
                    rejections.append(nextAtom)
                    currentPath.append(nextAtom)
                    atomicMass = molecWeight(molec, nextAtom)
                    newWeight = atomicMass + currentPath[0]
                    currentPath.appendleft(newWeight)
                    currentWeights.append(newWeight)
            # Break if all weights are unique
            if len(currentWeights) == len(set(currentWeights)):
                break
    
    return pathDict