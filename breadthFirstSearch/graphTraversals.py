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
def getListWeight(molec, atomIdxList):
    totalMass = 0.0
    for idx in atomIdxList:
        atom = molec.GetAtomWithIdx(idx)
        totalMass += atom.GetMass()
    return totalMass
def SmartNextStep(molec , g, currentAtom , contacts , badAtoms):
    currentEligibles = [contact for contact in contacts if contact not in badAtoms]
    if len(currentEligibles) == 0:
        #Terminate 
        return currentAtom , badAtoms
    else:
        atomHash = {}
        for contact in currentEligibles:
            weight = molecWeight(molec, contact)
            atomHash[contact] = weight
        atomHash = dict(sorted(atomHash.items(), key=lambda item: item[1], reverse=True))
        maxValue = max(atomHash.values())
        eligibleKeys = [key for key, val in atomHash.items() if val == maxValue]
        eligibles = []
        for key in eligibleKeys:
            contacts1 = list(g.neighbors(key))
            mass = getListWeight(molec , contacts1)
            eligibles.append(mass)
        maxIdx, maxVal = max(enumerate(eligibles), key=lambda x: x[1])
        badAtoms.append(currentAtom)
        return maxIdx , badAtoms
    
def reversePaths(atomHash): 
    #print("reversePath" , atomHash)
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
        #print("nextContacts" , nextContacts)
        nextAtom , nextWeight = max(nextContacts.items(), key=lambda x: x[1])
        del atomHash[nextKey]["openContacts"][nextAtom]
        return nextAtom , atomHash , nextWeight
    
def getListWeight(molec, atomIdxList):
    totalMass = 0.0
    for idx in atomIdxList:
        atom = molec.GetAtomWithIdx(idx)
        totalMass += atom.GetMass()
    return totalMass
def SmartNextStep(molec, g, currentAtom, contacts, badAtoms):
    # Filter out already-visited atoms
    currentEligibles = [contact for contact in contacts if contact not in badAtoms]
    
    if not currentEligibles:
        # No valid next atoms — terminate
        print(f"⚠️ No eligible atoms from {currentAtom}, terminating.")
        return currentAtom, badAtoms

    # Build atom weight mapping
    atomHash = {contact: molecWeight(molec, contact) for contact in currentEligibles}
    atomHash = dict(sorted(atomHash.items(), key=lambda item: item[1], reverse=True))

    # Get atoms with maximum single-atom weight
    maxValue = max(atomHash.values())
    eligibleKeys = [key for key, val in atomHash.items() if val == maxValue]

    # Evaluate neighborhood weights of these equally heavy atoms
    eligibles = []
    for key in eligibleKeys:
        contacts1 = list(g.neighbors(key))
        mass = getListWeight(molec, contacts1)
        eligibles.append(mass)

    # Choose the atom with the largest neighborhood weight
    maxIdx, maxVal = max(enumerate(eligibles), key=lambda x: x[1])
    nextAtom = eligibleKeys[maxIdx] 

    badAtoms.append(currentAtom)

    print(f"Next step from {currentAtom} → {nextAtom} (mass={maxVal:.3f})")
    return nextAtom, badAtoms
def heteroAtomContacts(molec, contacts1, contacts2 , badAtoms1 , badAtoms2):
    hetero1 = []
    hetero2 = []
    for idx in contacts1:
        atom = molec.GetAtomWithIdx(idx)
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H"):
            hetero1.append(idx)
    
    for idx in contacts2:
        atom = molec.GetAtomWithIdx(idx)
        symbol = atom.GetSymbol()
        if symbol not in ("C", "H"):
            hetero2.append(idx) 
    if len(hetero1) == len(hetero2) == 0:
        weights1 = getListWeight( molec, [contact for contact in contacts1 if contact not in badAtoms1])
        weights2 = getListWeight(molec, [contact for contact in contacts2 if contact not in badAtoms2])
        print(contacts1)
        print(contacts2)
        print(weights1 , weights2)
        if weights1 > weights2:
            return True , False
        elif weights1 < weights2:
            return False, True
        else:
            return False , False
    else:
        return False, False
def maxPathCompare(molec, graph, path1, path2, badAtoms1, badAtoms2, endSearch, depth=0):
    indent = "  " * depth  # for readable indentation per recursion level
    print(f"{indent}▶️ Depth {depth}: path1={path1}, path2={path2}, endSearch={endSearch}")

    weight1 = getListWeight(molec, path1)
    weight2 = getListWeight(molec, path2)
    print(f"{indent}Weights → path1: {weight1:.3f}, path2: {weight2:.3f}")

    # Compare path weights
    if weight1 > weight2:
        print(f"{indent}path1 heavier — returning [path1, path2]")
        return [path1, path2], False

    elif weight1 < weight2:
        print(f"{indent}path2 heavier — returning [path2, path1]")
        return [path2, path1], False

    # Equal weight — need to keep exploring
    else:

        if endSearch:
            return [path1, path2], True

        currentAtom1 = path1[-1]
        currentAtom2 = path2[-1]

        contacts1 = list(graph.neighbors(currentAtom1))
        print("contacts1" , contacts1)
        
        contacts2 = list(graph.neighbors(currentAtom2))
        print("contacts2" , contacts2)
        
        termPath1 , termPath2 = heteroAtomContacts(molec, contacts1 , contacts2 , badAtoms1 , badAtoms2)
        if termPath1:
            return [path1, path2], False
    
        elif termPath2:
            return [path2, path1], False


        
        else:
            nextAtom1, badAtoms1 = SmartNextStep(molec, graph, currentAtom1, contacts1, badAtoms1)
            nextAtom2, badAtoms2 = SmartNextStep(molec, graph, currentAtom2, contacts2, badAtoms2)
    
            if nextAtom1 == currentAtom1 or nextAtom2 == currentAtom2:
                print(f"{indent}No new atoms, ending search.")
                return [path1, path2], True
        
            else:
                # Extend paths and recurse
                new_path1 = path1 + [nextAtom1]
                new_path2 = path2 + [nextAtom2]
                print(f"{indent}Recursing deeper with new paths...")
                return maxPathCompare(molec, graph, new_path1, new_path2, badAtoms1, badAtoms2, endSearch, depth+1)






