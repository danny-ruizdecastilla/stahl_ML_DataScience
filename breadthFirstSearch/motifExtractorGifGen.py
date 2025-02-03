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
import random
import chemdraw
import plotly.graph_objs as go
from plotly_gif import GIF, two_d_time_series
from PIL import Image

#Given an Alkene and a SMILES returns a gif of the pathways taken to get the radius 
def getCC(smiles):
    #returns the index of C1 and C2 double bond 
    molec = Chem.MolFromSmiles(smiles)
    for bond in molec.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        double = 0
        #print(smiles)
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6 and bond.GetBondType() == Chem.BondType.DOUBLE:
            c1 = bond.GetBeginAtomIdx()
            c2 = bond.GetEndAtomIdx()
            CC = [c1 , c2]
            double +=1
    if double >1:
        return ["Error" , "Error"] , molec
    else:
        return CC , molec 
def getAdjencyMatrix(graph, radius):
    adjMatrix = nx.adjacency_matrix(graph).toarray()
    pathMatrix = np.linalg.matrix_power(adjMatrix, radius)
    return  pathMatrix 
def getUpperLimits(pathMatrix , C , graph , radius):
    limitList = []
    scope = pathMatrix[C] #column of contacts at the C (starting carbon) of interest
    for i in range (len(scope)):
        atom = scope[i]
        if atom == 1:
            #This is the hard limit of the scope
            limitList.append(i)
    if len(limitList) == 0:
        return ["Poison"] #Entire molecule can be kept 
    else:
        #maxLimit = len(limitList)
        #Trim down limitList 
        newlimitList = []
        for i in range (len(limitList)):
            atom2 = limitList[i]
            #print("atom2" , atom2)
            atom1atom2Paths = list(nx.all_simple_paths(graph, source=C, target=atom2))
            #print(atom1atom2Paths)
            withinRad = any(len(path) -1  < radius for path in atom1atom2Paths)
            if not withinRad:
                newlimitList.append(atom2)
            #print("70")
        return newlimitList 
def randomWalk(atom1 , prevAtom , graph):
    contacts = list(graph.neighbors(atom1))
    contacts = [c for c in contacts if c != prevAtom]
    if len(contacts) == 0:
        return "Nan"
    else:
        atom2 = random.choice(contacts)
        return atom2 
def getContacts(atom1, prevAtom, graph, limitList, CC, molec):
    #bond1 = getBondIndex(molec, prevAtom , atom1)
    #bondList = []
    #bondList.append(bond1)
    molecList = [prevAtom,atom1 ]
    while True:
        nextAtom = randomWalk(atom1 , prevAtom, graph)
        if nextAtom == "Nan" or nextAtom in CC:
            break 
        elif nextAtom in limitList:
            molecList.append(nextAtom)
            #nextBond = getBondIndex(molec, nextAtom , atom1)
            #bondList.append(nextBond)
            break
        else:    
            molecList.append(nextAtom)
            prevAtom = atom1
            atom1 = nextAtom 
            #nextBond = getBondIndex(molec, nextAtom , atom1)
            #bondList.append(nextBond)
    return molecList 
def getBondIndex(molec, atom1, atom2):
    bond = molec.GetBondBetweenAtoms(atom1, atom2)
    if bond is not None:
        return bond.GetIdx()
    return None
def motifExtract(smiles, radius):
    CC , molec = getCC(smiles)

    if not CC[0] == 'Error':
        g = Graph()
        #print("CC" , CC)
        for bond in molec.GetBonds():
            start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            g.add_edge(start, end)
        limitMatrix = getAdjencyMatrix( g , radius)
        pathListMAST = []
        
        for j in range (len(CC)):
            C1 = CC[j] #leading Carbon 
            C2 = CC[j-1] #left behind 
            #print("C1" , C1)
            #print("C2" , C2)
            upperLimits  = getUpperLimits(limitMatrix, C1 , g, radius)
            smilesEdges = upperLimits
            upperLimits.append(C2)
            upperLimits.append(C1)
            if upperLimits[0] == "Poison":
                print("Entire Molecule fits in path")
                sys.close()
            else:
                fail = 0
                while True:
                    path  = getContacts(C1 , C2, g, upperLimits , CC , molec)
                    #print(contacts)
                    if not path in pathListMAST:
                        #print(count)
                        pathListMAST.append(path)
                        
                        #bondListMAST.append(bonds)
                        fail = 0
                    else:
                        fail += 1
                    if fail == 20: #Statistically impossible to superceed with mw limit of 400 (need to fix for general form)
                        break
        #Now we have a list of all eligible paths in PathListMAST
        highlightList = [[CC[0] , CC[1]]]
        stepSize = 2
        
        while True:
            tempList = [CC[0] , CC[1]]
            for i in range (len(pathListMAST)):
                path = pathListMAST[i]
                #print(path)
                try:
                    atoms = path[2: stepSize+1]
                    #print(stepSize)
                    #print(path)
                    #print(atoms)
                    for i , atom in enumerate (atoms):
                        if atom not in tempList:
                            tempList.append(atom)

                except IndexError:
                    continue
            highlightList.append(tempList)
            stepSize +=1
            if stepSize == radius +2 :
                break

        return highlightList , molec

def figGenerator(highlight ,bonds,  smile, ind , pngDir):
    molec = chemdraw.Molecule(smile)
    accentColor = "rgb(252,186,63)"

    for id_ in highlight:
        molec.atoms[id_].highlight.show = True
        molec.atoms[id_].highlight.color = accentColor
    for bond in bonds:
        molec.bonds[bond].highlight.show = True
        molec.bonds[bond].highlight.color = accentColor


    drawer = chemdraw.Drawer(molec, title=smile)
    drawer.draw_img(pngDir + "/" + str(ind) + "molec.png")
            
def createGif(pngDir, outputFile, duration=500):

    # Get list of PNG files in directory
    pngFiles = glob.glob(pngDir + "/*.png")
    pngFiles.sort()  # Sort files to ensure consistent order
    
    # Open all PNG files
    images = []
    for png in pngFiles:
        img = Image.open(png)
        images.append(img)
    
    # Save the GIF
    if images:
        images[0].save(
            outputFile,
            save_all=True,
            append_images=images[1:],
            duration=duration,
            loop=0
        )
        print(f"GIF created successfully: {outputFile}")
        print(f"Number of frames: {len(images)}")
    else:
        print("No PNG files found in the specified directory")

def getBondIndexes(atomIndexes , molec):
    pairs = list(combinations(atomIndexes, 2))
    bondList = []
    for i in range (len(pairs)):
        atoms = pairs[i]
        bond = getBondIndex(molec, atoms[0], atoms[1])
        if bond:
            bondList.append(bond)
    return bondList


if __name__ == "__main__":
    alkeneSMILE = "COC1=CC=C(NS(=O)(N(C/C=C\C(C(C)(CCC(C)(CC)C(C)(C)C)C(C(C)(CC)C)C(C)CC(C(C)C)C)CC(C(C(C)C)C(C(C)C)CC)CC)C(C)(C)C)=O)C=C1"

    radius = int(sys.argv[1])
    pngDir = str(sys.argv[2])
    highlightList , molec  = motifExtract( alkeneSMILE , radius)
    print(highlightList)
    for i in range (len(highlightList)):
        bonds = getBondIndexes(highlightList[i] , molec)
        figGenerator(highlightList[i] ,bonds ,  alkeneSMILE , i , pngDir )#list of pngs 
    outputFile =   "alkene.gif"
    createGif(pngDir, outputFile, duration=500)
    

    