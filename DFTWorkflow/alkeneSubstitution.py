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

def getSubstitutionList(smilesList):
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
            contactList.pop(C2)
            if len(contactList) == 2:
                ccDict[C1] = contactList
        wildCards = [atom for sublist in ccDict.values() for atom in sublist]
        hydrogens = []
        for wildcard in wildCards:

def cisortrans(molec , graph , C1C2):
    positionsDict = {}
    while True:
        
def main(featureDir , outputDir):
    featureList = glob.glob("/*.csv")
    for featureDir in featureList:
        df = pd.read_csv(featureDir)
        stringCols = list(df.select_dtypes(include='object').columns)
        if "SMILES" in stringCols:
            smilesList = df["SMILES"]
            substitutionList = getSubstitutionList(smilesList)
        else:
            
if __name__ == "__main__":
    featureDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    main(featureDir)