import sys
import glob
import numpy as np
import os
from rdkit import Chem
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
import itertools


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
    scope = pathMatrix[C]
    #print("scope" , scope)
    for i in range (len(scope)):
        atom = scope[i]
        if atom == 1:
            #This is the hard limit of the scope
            limitList.append(i)
        if atom == 2:
            #test for symmetry 
            atom1atom2Paths = list(nx.all_simple_paths(graph, source=C, target=i))
            withinRad = any(len(path) -1  < radius for path in atom1atom2Paths)
            if not withinRad:
                limitList.append(i)
    if len(limitList) == 0:
        return ["Poison"] #Entire molecule can be kept 
    else:
        #maxLimit = len(limitList)
        #Trim down limitList 
        newlimitList = []
        #print("oldLimitList" , limitList)
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
def getContacts(atom1, prevAtom, graph, limitList, CC ): #(C1 , C2, g, upperLimits , CC , molec)
    #print("limitList" , limitList)
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
def getBondIndex(molec, atom1, atom2):
    bond = molec.GetBondBetweenAtoms(int(atom1), int(atom2))
    if bond is not None:
        return bond.GetIdx()
    return None
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
def smallestDistance(graph, edgeList , cutDist, atoms , compare):
    edgesMAST = []
    for edgeAtom in edgeList:
        proximityList = []
        for atom in atoms:
            atom1atom2Paths = list(nx.all_simple_paths(graph, source=atom, target=edgeAtom))
            for path in atom1atom2Paths:
                length = len(path)-1
                proximityList.append(length)
        #print("proximityList" , proximityList)
        #print("atoms" , atoms)
        if compare == ">=":
            condition = all(num >= cutDist for num in proximityList)
        elif compare == "<=":
            condition = all(num <= cutDist for num in proximityList)
        else:
            raise ValueError("Invalid compare argument. Use '>=' or '<='.")
        if condition:
            edgesMAST.append(edgeAtom)
    return edgesMAST
def removeProblemAroms(molec ,  contactList ,smilesEdges ):
    blanketList = contactList.copy()
    removeList = []
    removes = False
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
                #print("startAtom" , startAtom)
                atomList = []
                #walk the aromatic path
                prevAtom = -1
                while True:
                    atom = molec.GetAtomWithIdx(atomId)
                    #print("currentAtom" , atomId)
                    bonds = atom.GetBonds()
                    stepChoices = []
                    for bond in bonds:
                        atom1 = bond.GetBeginAtomIdx()  # Index of the first atom in the bond
                        atom2 = bond.GetEndAtomIdx()      # Index of the second atom in the bond
                        bondType = str(bond.GetBondType() )
                        atoms = [atom1 , atom2]
                        if bondType == 'AROMATIC' and prevAtom not in atoms:  # Take the next step 
                            # Next atom
                            if atom1 == atomId:  # Then it is atom_2 that is the next step
                                stepChoices.append([atom1 , atom2])

                            else:
                                stepChoices.append([atom2, atom1])
                    #print("stepChoices" , stepChoices)
                    if len(stepChoices) != 0:
                        
                        atoms = random.choice(stepChoices)
                        prevAtom = atoms[0]
                        atomId = atoms[1]
                        atomList.append(atomId)
                        if atomId == startAtom:
                            break 
                    else:
                        break
                atomSet = set(atomList)
                blanketList = list(set(blanketList) - atomSet) 
                missingElements = atomSet - set(contactList)
                #print("contactListMAST220" , contactList)
                if missingElements:
                    #print("missingElements")
                    #molec = removeArom(molec, atomList)
                    removeList.append(atomList)
                    removes = True
            else:
                blanketList.remove(atomId)
        else:
            blanketList.remove(atomId)

        if len(blanketList) == 0:
            break
    
    if removes:
        removeList = list(set(num for sublist in removeList for num in sublist))
        molec = removeArom(molec, removeList)
    return molec

def getMotif(edgeAtoms , graph , cutDist , subAtoms , molec , contactList ):
    smilesEdges = list(set(num for sublist in edgeAtoms for num in sublist))
    smilesEdges = smallestDistance(graph, smilesEdges , cutDist , subAtoms, ">=")
    molec1 = removeProblemAroms(molec , contactList , smilesEdges)
    editMolec = Chem.EditableMol(molec1)
    for i in range(len(smilesEdges)):
        starInd = smilesEdges[i]
        editMolec.ReplaceAtom(starInd, Chem.Atom("*"))
    molecMAST = editMolec.GetMol()
    motifSMILES = Chem.MolFragmentToSmiles(molecMAST, contactList, kekuleSmiles=False)
    return motifSMILES
def motifExtract(smiles, radius):
    CC , molec = getCC(smiles)

    if not CC[0] == 'Error':
        g = Graph()
        #print("CC" , CC)
        for bond in molec.GetBonds():
            start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            g.add_edge(start, end)
        limitMatrix = getAdjencyMatrix( g , radius)
        contactListMAST = []
        smilesEdgesMAST = []
        allPaths = []
        for j in range (len(CC)):
            contactList = []
            C1 = CC[j] #leading Carbon 
            C2 = CC[j-1] #left behind 
            #print("C1" , C1)
            #print("C2" , C2)
            upperLimits  = getUpperLimits(limitMatrix, C1 , g, radius)
            smilesEdges = upperLimits.copy()
            upperLimits.append(C2)
            upperLimits.append(C1)
            if upperLimits[0] == "Poison":
                print("Entire Molecule fits in path")
                sys.close()
            else:
                fail = 0
                while True:
                    contacts  = getContacts(C1 , C2, g, upperLimits , CC )
                    #print(contacts)
                    if not contacts in contactList:
                        #print(count)
                        contactList.append(contacts)
                        allPaths.append(contacts)
                        fail = 0
                    else:
                        fail += 1
                    if fail == 50: #Statistically impossible to superceed with mw limit of 400 (need to fix for general form)
                        break
            network = list(set(num for sublist in contactList for num in sublist))
            #print("network" , network)
            contactListMAST.append(network)
            smilesEdgesMAST.append(smilesEdges)
        #Now we have a list of all eligible paths in PathListMAST
        contactListMAST = list(set(num for sublist in contactListMAST for num in sublist))
        motifSMILE = getMotif(smilesEdgesMAST , g, radius , CC , molec , contactListMAST )
        highlightDict = getHighLightDict(radius,  CC, allPaths)
        return highlightDict , molec , motifSMILE
def getHighLightDict(cutDist,  initAtoms , contactPaths):
    highLightDictionary = {}
    i = 0
    #print(contactPaths)
    while True:
        if i == 0:
            #highlight only the initAtoms 
            highLightDictionary[255] = list(initAtoms)
        else:
            #print("Highlight i" , i)
            #print("Contact Paths" , contactPaths)
            iOut = []
            #print("iout" , iOut)
            scale = int(len(initAtoms)) -1 

            for contact in contactPaths:
                if len(contact) > (scale + i):
                    atom1 = contact[scale + i ]
                    iOut.append(atom1)
            iOut = list(np.unique(iOut))
            #print("iOut" , iOut)
            iOutMAST = []
            for atom in iOut:
                exists = any(atom in values for values in highLightDictionary.values())
                if not exists:
                    iOutMAST.append(atom)
        
            highLightDictionary[255 - int(i)*alphaScale] = list(iOutMAST)
        i += 1
        if i == cutDist + 1:
            break

    return highLightDictionary

def figGenerator(atomDict ,bondDict , smile):
    molec = chemdraw.Molecule(smile)
    print("NEW FRAME!!!!!!!!!!!!!!!!!!!!")
    ind = len(atomDict)
    for i, key in enumerate(atomDict):
        alpha = float(key)
        atoms = list(atomDict[key])
        #print("alpha" , alpha)
        print("atoms" , atoms)
        colorG = colorGradient[i]
        color = f"rgba({colorG[0]},{colorG[1]},{colorG[2]},{alpha})"
        print("atomColor" , color)
        for id_ in atoms:
            molec.atoms[id_].highlight.show = True
            molec.atoms[id_].highlight.color = color
        drawer = chemdraw.Drawer(molec, title=smile)
    for j , key in enumerate(bondDict):
        alpha = float(key)
        colorG = colorGradient[j]
        color = f"rgba({colorG[0]},{colorG[1]},{colorG[2]},{alpha})"
        print("bondColor" , color)
        bonds = list(bondDict[key])
        print("bonds" , bonds)
        for bond in bonds:
            molec.bonds[bond].highlight.show = True
            molec.bonds[bond].highlight.color = color
        drawer = chemdraw.Drawer(molec, title=smile)
    drawer.draw_img(pngDir + "/" + str(ind) + str(smile) + "molec.png")
def finalImage(smile):
    molec = Chem.MolFromSmiles(smile)
    Chem.rdDepictor.Compute2DCoords(molec)
    Chem.Draw.MolToFile(molec , pngDir + "/subMolec"  + str(smile) + "Finalmolec.png")

def createGif(pngDir, outputFile, duration=500):
    # Get list of PNG files in directory
    pngFiles = glob.glob(pngDir + "/*.png")
    pngFiles.sort()  # Sort files to ensure consistent order
    
    # Open all PNG files
    images = []
    for png in pngFiles:
        img = Image.open(png)
        images.append(img)
    
    if images:
        # Determine target size based on the first image
        target_size = images[0].size
        
        # Resize all images to the target size using LANCZOS filter
        resized_images = [img.resize(target_size, Image.Resampling.LANCZOS) for img in images]
        
        # Save the GIF using the resized images
        resized_images[0].save(
            outputFile,
            save_all=True,
            append_images=resized_images[1:],
            duration=duration,
            loop=0
        )
        print(f"GIF created successfully: {outputFile}")
        print(f"Number of frames: {len(resized_images)}")
    else:
        print("No PNG files found in the specified directory")

def getBondIndexes(dict1 ,alphas, atomIndexes , molec):
    #ouut: { a : [1,2,3] , b : [4,5,6]}
    pairs = list(combinations(atomIndexes, 2))
    bondDict = {}
    dictLen = len(dict1)
    if  255 in dict1:
        mainSub = list(dict1[255])
        pairs = list(combinations(mainSub, 2))
        mainSubBonds = []
        for pair in pairs:
            isPair = getBondIndex(molec, pair[0], pair[1])
            if isPair:
                mainSubBonds.append(isPair)
    
        bondDict[255] = mainSubBonds
    
    if dictLen > 1:
        #Add bonds between radii
        for i in range (len(atomIndexes) - 1):
            mainSubBonds = []
            alpha = alphas[i+1]
            #print("alphas" , alpha)
            atom1 = atomIndexes[i]
            atom2 = atomIndexes[i+1]
            pairs = [[a1,a2] for a1 in atom1 for a2 in atom2]
            for pair in pairs:
                isPair = getBondIndex(molec, pair[0] , pair[1])
                if isPair:
                    mainSubBonds.append(isPair)
            bondDict[alpha] = mainSubBonds
        
    return bondDict

def decomposeDictionaries(dict1):
    dictList = []
    #decomposes dictionaries as so:
    #input:{A : 1 , B : 2, C : 3}
    #output: [{A:1} , {A:1 , B : 2} , {A:1 , B : 2 , C :3}]
    keyRanges = [] #range of keys needed for each subDict
    keyList = list(dict1.keys())
    for i in range (len(dict1)):
        numRange = list(np.arange(0 , i + 1 , 1 ))
        keys = []
        for num in numRange:
            keys.append(keyList[num])
        keyRanges.append(keys)
    for keys in keyRanges:
        dictList.append({k:dict1[k] for k in keys if k in dict1})
    return dictList
def lightGradient(initRGB , numSteps):
    base = np.array(initRGB)
    white = np.array([255 ,255, 255])
    return np.linspace(base , white, numSteps)

if __name__ == "__main__":
    alkeneSMILE = "[H][C@@]12C[C@@](N([C@H](C)C3=CC=CC=C3)[C@@H]2C(OCC4=CC=CC=C4)=O)([H])C=C1"
    #RGB: 153,000,000
    radius = int(sys.argv[1])
    pngDir = str(sys.argv[2])
    alphaScale = float(sys.argv[3])
    rgbInit = [197 , 5 ,12]
    colorGradient = lightGradient(rgbInit , radius + 2)
    highlightDict , molec   , finalSMILE = motifExtract( alkeneSMILE , radius)
    print(highlightDict)
    #print(highlightDict)
    #print("finalSMILE" , finalSMILE)

    dictList = decomposeDictionaries(highlightDict)
    #print(dictList)
    for dict in dictList:
        #print("dict" , dict)
        alphaVals = list(dict.keys())
        atomInds = list(dict.values())
        #print("atomInds" , atomInds)
        bondDict = getBondIndexes(dict, alphaVals, atomInds, molec)
        #print("439" , bondDict)
        figGenerator(dict, bondDict , alkeneSMILE)


    #finalImage(finalSMILE)
    outputFile =  pngDir +  "/" + str(alkeneSMILE) + "alkene.gif"
    createGif(pngDir, outputFile, duration=500)
    

    