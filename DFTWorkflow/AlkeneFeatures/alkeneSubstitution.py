#Danny Ruiz de Castilla 06.22.2025
import sys
import glob
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.PandasTools import LoadSDF
import networkx as nx
from networkx import Graph
import random
import chemdraw
import pandas as pd
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from DFTWorkflow.featureMaping import  createCSV
from breadthFirstSearch.radialBasedCorrelation import getCC,getAdjencyMatrix
from breadthFirstSearch.graphTraversals import minFirstSearch,maxPathCompare

def getSubstitutionList(smilesList):
    hCount = []
    cisTransList = []
    for smiles in smilesList:
        cc , molec = getCC(smiles)
        molec = Chem.AddHs(molec) 
        g = Graph()
        #print("CC" , CC)
        for bond in molec.GetBonds():
            start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            g.add_edge(start, end)
        limitMatrix = getAdjencyMatrix( g , 1)
        ccDict = {}
        for j in range (len(cc)): #Need to extend network to each carbon in the Alkene 
            #print(j)
            contactList = []
            C1 = cc[j] #leading Carbon 
            C2 = cc[j-1] #left behind 
            contacts = limitMatrix[C1]
            for i in range (len(contacts)):
                atom = contacts[i]
                if atom == 1:
                    #This is the hard limit of the scope
                    contactList.append(i)
            if C2 in contactList:
                contactList.remove(C2)
            if len(contactList) == 2:
                ccDict[C1] = contactList
        wildCards = [atom for sublist in ccDict.values() for atom in sublist]
        hydrogens = 0
        for wildcard in wildCards:
            atom = molec.GetAtomWithIdx(wildcard)
            atom1 = str(atom.GetSymbol())
            if atom1 == "H":
                hydrogens += 1
        if hydrogens <=2:
            print(smiles)
            addended = eVszAlkenes(smiles, g, ccDict)
            cisTransList.append(addended)
        else:
            cisTransList.append(0)
        hCount.append(hydrogens)
         
    return hCount , cisTransList
def smiles_to_coords(smiles):
    try:
        m = Chem.MolFromSmiles(smiles)
    except:
        print("could not convert %s to rdkit molecule. Exit!"%(smiles))
        exit()
    try:
        m = Chem.AddHs(m)
    except:
        print("ERROR: could not add hydrogen to rdkit molecule of %s. Exit!"%(smiles))
        exit()
    try:
        AllChem.EmbedMolecule(m)
    except:
        print("ERROR: could not calculate 3D coordinates from rdkit molecule %s. Exit!"%(smiles))
        exit()
    try:
        block = Chem.MolToMolBlock(m)
        blocklines=block.split("\n")
        coords = []
        elements = []
        for line in blocklines[4:]:
            if len(line.split()) == 4:
                break
            elements.append(line.split()[3])
            coords.append([float(line.split()[0]),float(line.split()[1]) , float(line.split()[2])])
        coords = np.array(coords)
        return elements , coords
    except:
        print("Could not convert SMILES to .xyz")
        exit()
     
def eVszAlkenes(molec , graph , ccDict):

    idxMAST = []
    for carbon , startAtoms  in ccDict.items():
        rejections = [carbon , startAtoms[0] , startAtoms[1]]
        pathDict = maxPathCompare(molec, graph, startAtoms, rejections)
        print(pathDict)
        weights = [dq[0] for dq in pathDict.values()]
        if weights[0] == weights[1]:
            #same weights, terminal alkene 
            return 0.25
        else: 
            for atom , dq in pathDict.items():
                lstdq = list(dq)
                midpoint = len(lstdq) // 2
                idxList= lstdq[midpoint:]
                for idx in idxList:
                    idxMAST.append(int(idx))
            idxMAST.append(int(carbon))
    idxSet = list(set(idxMAST))

    smiles1 = Chem.MolFragmentToSmiles(molec, idxSet, kekuleSmiles=False)
    if "." in smiles1:
        raise ValueError(". in the subMolec, incorrect path")
    else:
        print(smiles1)
        slashCount= smiles1.count('/')
        backslashCount = smiles1.count('\\')
        if slashCount == 2:
            #Z alkene with Z brackets 
            return 0.5
        if backslashCount == 2:
            #Z alkene with Z brackets 
            return 0.5
        if slashCount == 1 and backslashCount == 1:
            #E Alkene
            return 0.75

  
def main(featureDir , outputDir):
    featureList = glob.glob(featureDir + "/*.csv")
    for featureFile in featureList:
        df = pd.read_csv(featureFile)
        stringCols = list(df.select_dtypes(include='object').columns)
        if "SMILES" in stringCols:
            smilesList = df["SMILES"]
            hCount , cisVstrans = getSubstitutionList(smilesList)
            df["HCount"] = hCount
            df["cisVstrans"] = cisVstrans
        else:
            raise ValueError("SMILES column not found in DataFrame")
        createCSV(df , outputDir , outputDir.split("/")[-1].split(".")[0])

        
if __name__ == "__main__":
    featureDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    main(featureDir , outputDir)