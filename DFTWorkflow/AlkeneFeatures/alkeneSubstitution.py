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
    EvsZList = []
    for smiles in smilesList:
        cc , molec = getCC(smiles)
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
        hCount.append(hydrogens)
        if hydrogens <=2: #Easy check for E vs Z 
            alkeneType = eVszAlkenes(molec , g , ccDict)# 0 if neither E or Z, -1 if Z, 1 if E
        else:
            alkeneType = 0
        EvsZList.append(alkeneType)          
    return hCount , EvsZList
def smiles_to_coords(smiles):
    try:
        m = Chem.MolFromSmiles(smiles)
    except:
        print("could not convert %s to rdkit molecule. Exit!"%(smiles))
        exit()
    try:
        m = Chem.AddHs(m)
    except:
        print("Warning: could not add hydrogen to rdkit molecule of %s. Exit!"%(smiles))
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
     
def eVszAlkenes(molec , graph , C1Hash , C2Hash , coordinates):
    #C1Hash = {"C1" : [cont1 , cont2]}
    c1Atom = list(C1Hash.keys())[0]
    c1Contacts = C1Hash[c1Atom]

    c2Atom = list(C2Hash.keys())[0]
    c2Contacts = C2Hash[c2Atom]

    pathRank1 , same1 = maxPathCompare(molec, graph , [c1Contacts[0]] , [c1Contacts[1]] , [c1Atom , c1Contacts[0]] , [c1Atom , c1Contacts[1]] , False)
    pathRank2 , same2 = maxPathCompare(molec, graph , [c2Contacts[0]] , [c2Contacts[1]] , [c2Atom , c2Contacts[0]] , [c2Atom , c2Contacts[1]] , False)
    if same1 or same2:
        #Gem substituted, both routs are the same
        return -1
    else:
        maxWildC1 = pathRank1[0][0]

        maxWildC2 = pathRank2[0][0]
        minWildC2 = pathRank2[1][0]
        
        atoms = list(coordinates.keys())
        maxAtomWildC1 = coordinates[atoms[maxWildC1]][2:5]
        maxAtomWildC2 = coordinates[atoms[maxWildC2]][2:5]
        minAtomWildC2 = coordinates[atoms[minWildC2]][2:5]

        distanceC1_maxC2 = np.linalg.norm(np.array(maxAtomWildC1 , dtype=float) - np.array(maxAtomWildC2 , dtype=float))
        distanceC1_minC2 = np.linalg.norm(np.array(maxAtomWildC1 , dtype=float) - np.array(minAtomWildC2 , dtype=float))

        if distanceC1_maxC2 >= distanceC1_minC2:
            #E Alkene
            return 0
        else:
            return 1
            #Z Alkene
def main(featureDir , outputDir):
    featureList = glob.glob(featureDir + "/*.csv")
    for featureFile in featureList:
        df = pd.read_csv(featureFile)
        stringCols = list(df.select_dtypes(include='object').columns)
        if "SMILES" in stringCols:
            smilesList = df["SMILES"]
            hCount , EvsZ = getSubstitutionList(smilesList)
            df["HCount"] = hCount
            df["EvzZ"] = EvsZ
        else:
            raise ValueError("SMILES column not found in DataFrame")
        createCSV(df , outputDir , outputDir.split("/")[-1].split(".")[0])

        
if __name__ == "__main__":
    featureDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    main(featureDir , outputDir)