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
    firstPrinciples = [dq[0] for dq in pathDict.values()]
    pathSame = all(val == firstPrinciples[0] for val in firstPrinciples)

    if not pathSame:

        return pathDict
    else:

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
def checkWeights(hash1 , hash2 , tolerance):
    print(hash1)
    print(hash2)
    idx1 = list(hash1.keys())[-1]
    idx2 = list(hash2.keys())[-1]
    weight1 = float(hash1[idx1]["PrevWeight"])
    weight2 = float(hash2[idx2]["PrevWeight"])
    atom1 = float(hash1[idx1]["atomWeight"])
    atom2 = float(hash2[idx2]["atomWeight"])
    diff = abs(weight1 - weight2)
    diff1 = abs(atom1 - atom2)
    if diff >= tolerance or diff1 >= tolerance:
        return True
    else:
        return False
def getEligibleAtoms(molec ,graph,  contacts , badAtoms):
    currentEligibles = [contact for contact in contacts if contact not in badAtoms]
    #print(currentEligibles)
    if currentEligibles != 0:
        newEligibles = {}
        contactWeight = 0
        protonCount = 0
        for atom in currentEligibles:
            symbol = molec.GetAtomWithIdx(atom)
            if symbol != "H": #avoid protons
                contactWeight += molecWeight(molec, atom)
                newEligibles[atom] = molecWeight(molec, atom)
            else:
                protonCount += 1
        if protonCount == 3: #terminal carbon
            terminal = True
        else:
            terminal = False

        return contactWeight, newEligibles , terminal
    else:
        return 12.01 , contacts , True
def reversePaths(atomHash): 
    nextoptions = []
    for atom , info in atomHash.items():
        availableContacts = list(info["openContacts"])
        if len(availableContacts) != 0:
            nextoptions.append(atom)
    if len(nextoptions) == 0:
        return None , None , None
    else:
        nextKey = nextoptions[-1]
        nextContacts = atomHash[nextKey]["openContacts"]
        nextAtom = max(nextContacts.items(), key=lambda x: x[1])
        nextWeight = float(nextContacts[nextAtom])
        del atomHash[nextKey]["openContacts"][nextAtom]
        return nextAtom , atomHash , nextWeight
    
def EvsZCompare(molec, g , path1 , path2 , badAtoms1,badAtoms2, diffTolerance):
    print(badAtoms1)
    isDifferent = checkWeights(path1, path2, diffTolerance)
    
    if isDifferent:
        return path1 , path2  
        print("isDifferentExit")
    else: 
        currentAtom1 = list(path1.keys())[-1]
        print(currentAtom1)
        contacts1 = list(g.neighbors(currentAtom1))
        contactWeights1 , eligibleAtoms1 , isTerminal1 = getEligibleAtoms(molec, g, contacts1 , badAtoms1)
        #print(eligibleAtoms1)
        if isTerminal1:#get nextAtom from another path
            nextAtom1 , path1 , nextWeight1 = reversePaths(path1)
            if nextAtom1 is None:
                print("nextAtom1Exit")
                return path1 , path2
        else:
            nextAtom1 , nextWeight1 = max(eligibleAtoms1.items(), key=lambda x: x[1])
            del eligibleAtoms1[nextAtom1]
        badAtoms1.append(nextAtom1)
        
        currentAtom2 = list(path2.keys())[-1]
        contacts2 = list(g.neighbors(currentAtom2))
        contactWeights2 , eligibleAtoms2 , isTerminal2  = getEligibleAtoms(molec,g ,  contacts2 , badAtoms2)
        if isTerminal2:#get nextAtom from another path
            nextAtom2 , path2 , nextWeight2 = reversePaths(path2)    
            if nextAtom2 is None:
                print("nextAtom2Exit")
                return path1 , path2
        else:
            nextAtom2  , nextWeight2 = max(eligibleAtoms2.items(), key=lambda x: x[1])
            del eligibleAtoms2[nextAtom2]
        badAtoms2.append(nextAtom2)

        path1[nextAtom1] = {"openContacts": eligibleAtoms1 , 
                            "PrevWeight" : contactWeights1 , "atomWeight" : nextWeight1 }
        path2[nextAtom2] = {"openContacts": eligibleAtoms2 , 
                            "PrevWeight" : contactWeights2 , "atomWeight" : nextWeight2 }
        print("nextLoop")
        return EvsZCompare(molec, g , path1 , path2 , badAtoms1,badAtoms2, diffTolerance)
