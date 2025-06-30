#Danny Ruiz de Castilla 06.22.2025
import sys
import glob
import numpy as np
import os
from rdkit import Chem
from rdkit.Chem.PandasTools import LoadSDF
import networkx as nx
from networkx import Graph
import random
import chemdraw
import pandas as pd
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from breadthFirstSearch.radialBasedCorrelation import getCC,getAdjencyMatrix
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

def getSubstitutionList(smilesList):
    substitutionList = []
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
            contactList.pop(C2)
            if len(contactList) == 2:
                ccDict[C1] = contactList
        wildCards = [atom for sublist in ccDict.values() for atom in sublist]
        hydrogens = 0
        for wildcard in wildCards:
            atom = molec.GetAtomWithIdx(wildcard)
            atom1 = str(atom.GetSymbol())
            if atom1 == "H":
                hydrogens += 1
        if hydrogens == 2:
            
            addended = cistransterminal(molec, g, ccDict)
            hydrogens += addended
        substitutionList.append(hydrogens)
         
    return substitutionList

def cistransterminal(molec , graph , ccDict):
    carbonDict = {}
    for carbon , startAtoms  in ccDict.items():
        pathDict = minFirstSearch(startAtoms , graph , list(ccDict.keys()))#gives path 1, path 2 centered on Alkene carbon of interest
        motifList = []
        for atom , path in pathDict.items():
            motifSMILES = Chem.MolFragmentToSmiles(molec, path, kekuleSmiles=False)
            motifList.append(motifSMILES)
        if str(motifList[0][1]) == str(motifList[1][1]) == "H":
            #alkene is terminal 
            return 0
        else:
            carbonDict[carbon] = motifList
    
    

  
def main(featureDir , outputDir):
    featureList = glob.glob("/*.csv")
    for featureDir in featureList:
        df = pd.read_csv(featureDir)
        stringCols = list(df.select_dtypes(include='object').columns)
        if "SMILES" in stringCols:
            smilesList = df["SMILES"]
            substitutionList = getSubstitutionList(smilesList)
            df["alkeneSubstitution"] = substitutionList
        else:
            raise ("SMILES column not found in Dataframe")
        
if __name__ == "__main__":
    featureDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    main(featureDir)