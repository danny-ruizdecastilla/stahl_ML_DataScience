import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
from rdkit.Chem.PandasTools import LoadSDF

def initAlkene(smiles , addHs):
    molec = Chem.MolFromSmiles(smiles)
    double = 0
    for bond in molec.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6 and bond.GetBondType() == Chem.BondType.DOUBLE:
            double +=1
            c1 = bond.GetBeginAtomIdx()
            c2 = bond.GetEndAtomIdx()
            CC = [c1 , c2]
    if double != 1:
        raise SystemError("More than 1 Alkene found in system")
    elif addHs:
        mol_with_hs = Chem.AddHs(molec)
        return CC , mol_with_hs
    else:
        return CC , molec 
def breadthFirstSearch(molec, period ,  contactHash , forbiddenAtoms ): #contactHash : { 0 : [startAtom] , 1 : [id1, id2, id3] , 2 : [id4,id5...]}
    currentPeriod = list(contactHash.keys())[-1]
    #print(contactHash)
    if currentPeriod > period:
        return contactHash
    else:
        prevTraversed = set(atom for neighbors in contactHash.values() for atom in neighbors)
        if len(forbiddenAtoms) > 0:
            prevTraversed.update(forbiddenAtoms)
        eligibleIdxs = []
        for idx in contactHash[currentPeriod]:
            atom = molec.GetAtomWithIdx(idx)
            neighbors = atom.GetNeighbors()
            for nbr in neighbors:
                nbrID = nbr.GetIdx()
                if nbrID not in prevTraversed:
                    eligibleIdxs.append(nbrID)
        eligibleIds = set(eligibleIdxs)
        newPeriod = currentPeriod +1 
        contactHash[newPeriod] = list(eligibleIds)
        return breadthFirstSearch(molec, period, contactHash, forbiddenAtoms)