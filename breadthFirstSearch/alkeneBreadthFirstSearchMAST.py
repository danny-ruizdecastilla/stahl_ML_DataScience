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
def breadthFirstSearch(molec, period, contactHash, forbiddenAtoms):
    # contactHash: { 0: [[startAtom]], 1: [[startAtom, id1], [startAtom, id2]], 2: [[startAtom, id1, id4], ...] }
    currentPeriod = list(contactHash.keys())[-1]

    if currentPeriod > period:
        return contactHash
    else:
        prevTraversed = set(atom for paths in contactHash.values() for path in paths for atom in path)
        if len(forbiddenAtoms) > 0:
            prevTraversed.update(forbiddenAtoms)

        eligiblePaths = []
        for path in contactHash[currentPeriod]:
            lastAtom = path[-1]  # the frontier atom of this trajectory
            atom = molec.GetAtomWithIdx(lastAtom)
            neighbors = atom.GetNeighbors()
            for nbr in neighbors:
                nbrID = nbr.GetIdx()
                if nbrID not in prevTraversed:
                    newPath = path + [nbrID]  # extend trajectory
                    eligiblePaths.append(newPath)
                    prevTraversed.add(nbrID)  # mark as visited within this period expansion

        newPeriod = currentPeriod + 1
        contactHash[newPeriod] = eligiblePaths
        return breadthFirstSearch(molec, period, contactHash, forbiddenAtoms)