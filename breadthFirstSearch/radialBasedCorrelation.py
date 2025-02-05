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
#Danny Ruiz de Castilla 01/12/2024
'''
Goal of this program is to take in a .xlsx and return a series of smaller .xlsx that groups common structures that are common up to some cutoff radial distance 
Eliminate any duplicates that were formed out of the post processing as well
''' 
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
def isCyclic(graph, atom1, atom2):
    atom1atom2Paths = list(nx.all_simple_paths(graph, source=atom1, target=atom2))
    if len(atom1atom2Paths) == 1:
        #Not Cyclic 
        cyclic = False
    else:
        cyclic = True
    return cyclic

def getAdjencyMatrix(graph, radius):
    adjMatrix = nx.adjacency_matrix(graph).toarray()
    pathMatrix = np.linalg.matrix_power(adjMatrix, radius)
    return  pathMatrix 
def getUpperLimits(pathMatrix , C , graph , radius):
    limitList = []
    scope = pathMatrix[C]
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

def getContacts(atom1, prevAtom, graph, limitList, CC):

    molecList = [prevAtom,atom1 ]
    while True:
        nextAtom = randomWalk(atom1 , prevAtom, graph)
        if nextAtom == "Nan" or nextAtom in CC:
            break 
        elif nextAtom in limitList:
            molecList.append(nextAtom)
            break
        else:    
            molecList.append(nextAtom)
            prevAtom = atom1
            atom1 = nextAtom 
    return molecList
def motifExtraction(SMILES, indexSMILES ,radius):
    motifMASTER = []
    smilesIndeces = []
    failCount = 0
    for smileInd, smiles in enumerate(list(SMILES)):
        print("i" , smileInd)
        
        print(smiles)
        try:

            CC , molec = getCC(smiles)
            if CC[0] == "Error":
                continue
            else:
                g = Graph()
                #print("CC" , CC)
                for bond in molec.GetBonds():
                    start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
                    g.add_edge(start, end)
                limitMatrix = getAdjencyMatrix( g , radius)
                contactListMAST = []
                smilesEdgesMAST = []
                for j in range (len(CC)): #Need to extend network to each carbon in the Alkene 
                    #print(j)
                    contactList = []
                    C1 = CC[j] #leading Carbon 
                    C2 = CC[j-1] #left behind 
                    #print("C1" , C1)
                    #print("C2" , C2)
                    upperLimits  = getUpperLimits(limitMatrix, C1 , g, radius)
                    smilesEdges = upperLimits.copy()
                    print("smilesEdges127" , smilesEdges)
                    upperLimits.append(C2)
                    upperLimits.append(C1)
                    #print("upper limits",  upperLimits)
                    #print("Passed Upper Limits")
                    if upperLimits[0] == "Poison":
                        #Entire molecule can be kept 
                        motifMASTER.append(smiles)
                        identSmiles = str(list(indexSMILES)[smileInd])
                        smilesIndeces.append(identSmiles)
                        bypass = True
                        break
                    else: 
                        bypass = False
                        fail = 0
                        while True:
                            print("fail" , fail)
                            contacts= getContacts(C1 , C2, g, upperLimits , CC)
                            #print(contacts)
                            if not contacts in contactList:
                                #print(count)
                                contactList.append(contacts)
                                fail = 0
                            else:
                                fail += 1
                            if fail == 20: #Output should include stars and aromatics, then new outputs can clean up and add either all Hydrogens or preserve double bond network for searches 
                                break
                    network = list(set(num for sublist in contactList for num in sublist))
                    #print("network" , network)
                    contactListMAST.append(network)
                    smilesEdgesMAST.append(smilesEdges)
                smilesEdgesMAST = list(set(num for sublist in smilesEdgesMAST for num in sublist))
                contactListMAST = list(set(num for sublist in contactListMAST for num in sublist))
                smilesEdgesMAST = smallestDistance(g, smilesEdgesMAST , radius , CC)

                if not bypass:
                    #correct the 2 contact Lists and turn to motifs 
                    molec = removeProblemAroms(molec, contactListMAST ,  smilesEdgesMAST)
                    editMolec = Chem.EditableMol(molec)
                    print("smilesEdges163" , smilesEdgesMAST)

                    for i in range (len(smilesEdgesMAST)):
                        starInd = smilesEdges[i]
                        editMolec.ReplaceAtom(starInd, Chem.Atom("*"))
                    molecMAST = editMolec.GetMol()
                    
                    motifSMILES = Chem.MolFragmentToSmiles(molecMAST, contactListMAST, kekuleSmiles=False)

                    #Post Processing of motifs 
                    motifMASTER.append(motifSMILES)
                    identSmiles = str(list(indexSMILES)[smileInd])
                    smilesIndeces.append(identSmiles)
        except AttributeError:
            #Cursed SMILES
            failCount +=1
            continue
    print("failCount" , failCount)
    return motifMASTER , smilesIndeces
#Post Processing of motifs
def smallestDistance(graph, edgeList , cutDist, atoms):
    edgesMAST = []
    for atom in atoms:
        #index of atom, check if distance of all in edgeList is greater than or equal to cutDist

    return edgesMAST
def removeProblemAroms(molec ,  contactList ,smilesEdges ):
    blanketList = contactList.copy()
    
    while True:
        #print("contactListMAST183" , contactList)
        #print(blanketList)
        atomId = random.choice(blanketList)
        atom = molec.GetAtomWithIdx(atomId)
        if atom.GetIsAromatic():
            #check if atom is in contactListMAST
            if atomId in smilesEdges:
                #return network of smilesString
                startAtom = atomId
                atomList = []
                #walk the aromatic path
                prevAtom = -1
                while True:
                    atom = molec.GetAtomWithIdx(atomId)
                    bonds = atom.GetBonds()
                    for bond in bonds:
                        atom1 = bond.GetBeginAtomIdx()  # Index of the first atom in the bond
                        atom2 = bond.GetEndAtomIdx()      # Index of the second atom in the bond
                        bondType = str(bond.GetBondType() )
                        atoms = [atom1 , atom2]
                        if bondType == 'AROMATIC' and prevAtom not in atoms:  # Take the next step 
                            # Next atom
                            if atom1 == atomId:  # Then it is atom_2 that is the next step
                                atomList.append(atom2)
                                prevAtom = atom1
                                atomId = atom2
                            else:
                                atomList.append(atom1)
                                prevAtom = atom2
                                atomId = atom1
                            #print("currentAtom", atomId)
                            break
                    if atomId == startAtom:
                        break 
                atomSet = set(atomList)
                blanketList = list(set(blanketList) - atomSet) 
                missingElements = atomSet - set(contactList)
                #print("contactListMAST220" , contactList)
                if missingElements:
                    #print("missingElements")
                    molec = removeArom(molec, atomList)
            else:
                blanketList.remove(atomId)
        else:
            blanketList.remove(atomId)

        if len(blanketList) == 0:
            break
    return molec
def removeArom(molec, indices):
    editableMol = Chem.RWMol(molec)

    for i in range(len(indices)):
        atom1 = indices[i-1]
        atom2 = indices[i]
        bond = editableMol.GetBondBetweenAtoms(atom1, atom2)
        if bond is not None and bond.GetIsAromatic():
            bond.SetBondType(Chem.BondType.SINGLE)
            bond.SetIsAromatic(False)
    for idx in indices:
        atom = editableMol.GetAtomWithIdx(idx)
        if atom.GetIsAromatic():
            atom.SetIsAromatic(False)

    Chem.SanitizeMol(editableMol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    
    # Return the modified molecule
    return editableMol

def createXLSX(smilesList ,smilesID ,  radius):
    columns = ["SMILES", "OriginalID"]
    # Create the DataFrame
    df = pd.DataFrame({columns[0]: smilesList, columns[1]: smilesID})
    saveDir = mainDir  + "/motifExtractorOutput/"
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    df.to_excel(saveDir + "breadthFirstSearchOuts_" + str(radius) + ".xlsx", index=False) 

if __name__ == "__main__":
    mainDir = str(sys.argv[1])
    radius = int(sys.argv[2])
    inputDF = str(sys.argv[3]) #name of the input DF
    
    dF = pd.read_excel(mainDir + "/" + inputDF + ".xlsx")
    dfMAST = dF.drop_duplicates(subset="SMILES", keep="first") 

    smilesDF = dfMAST['SMILES']
    indexingDF = dfMAST['Compound_Name']

    motifList , smilesID = motifExtraction(smilesDF, indexingDF, radius) #SmilesID keeps track of the original 
 
    createXLSX(motifList , smilesID , radius)
    #visulizations!!!!!

    