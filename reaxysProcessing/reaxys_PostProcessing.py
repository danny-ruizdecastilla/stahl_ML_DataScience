import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
from rdkit.Chem.PandasTools import LoadSDF
#Danny Ruiz de Castilla 01/09/2024

'''
Goal of this program: To take any substrates with double nonconjugated double bonds and decompose them into SMILES strings with individual unique double bonds 

Takes in 1 .txt of SMILES strings and returns a .xlsx with cleaned up SMILES and naming structure 

SMILES	                  Compound_Name
C=CC1=CC=CC=C1	          1_1_Alk_000000
CC(C)(C)C1=CC=C(C=C1)C=C  1_1_Alk_000001
C=CC1=CC=C(C)C=C1	      1_1_Alk_000002
C=CC1=CC=C(CCl)C=C1	      1_1_Alk_000003
C=CC1=CC=C(Cl)C=C1	      1_1_Alk_000004
C=CC1=CC=C(Br)C=C1	      1_1_Alk_000005
C=CC1=CC=CC=C1Cl	      1_1_Alk_000006
C=CC1=CC=CC=C1Br	      1_1_Alk_000007
C=CC1=CC=C(F)C=C1	      1_1_Alk_000008
'''

def saveSMILES(smiles , savePath): #savePath is a .dat that temporarily saves stuff
    if not os.path.exists(savePath):
        with open(savePath, "w") as file:
            #print("28")
            file.write(smiles + "\n")
    else:
        with open(savePath, "a") as file:
            #print("33")
            file.write(smiles + "\n")
    with open(savePath , "r+") as file:
        smilesList = [line.strip() for line in file]
    print("Progress:" , smilesList[-1])
def createXLSX(savePath):
    with open(savePath, 'r+') as file:
        smilesList = [line.strip() for line in file]
    #print("Output:" , smilesList[-1])
    smilesArray = np.array(smilesList)
    dim1 = len(smilesList)
    numbers = []
    for i in range(dim1):  # Generate numbers from 1 to 10
        numb = str(i).zfill(6)
        numbs = "Alkene_" + numb
        numbers.append(numbs)
    columns = ["SMILES", "Compound_Name"]
    # Create the DataFrame
    df = pd.DataFrame({columns[0]: smilesArray, columns[1]: numbers})

    df.to_excel(mainDir  + "/reaxys_Alkene_database.xlsx", index=False) 

def checkForFlaws(smiles):
    status = 1
    badStrings = [ '13C' ,'14C' ,  '2H' , '3H' , '15N' , '17O' , '18O'] #if it sees any of these, status turns to 0
    badAtoms = ['Pd' , 'Fe' , 'Ni' , 'Pt' , 'Mn' , 'Cu' , 'Zn']
    if any (bad2 in smiles for bad2 in badAtoms):
        status = 0

    elif any(bad in smiles for bad in badStrings):
        molec2 = Chem.MolFromSmiles(smiles)
        for atom in molec2.GetAtoms():
            if atom.GetAtomicNum() == 6 and atom.GetIsotope() == 13:
                atom.SetIsotope(0)
            if atom.GetAtomicNum() == 6 and atom.GetIsotope() == 14:
                atom.SetIsotope(0)
            if atom.GetAtomicNum() == 1 and atom.GetIsotope() == 2:
                atom.SetIsotope(0)
            if atom.GetAtomicNum() == 1 and atom.GetIsotope() == 3:
                atom.SetIsotope(0)
            if atom.GetAtomicNum() == 7 and atom.GetIsotope() == 15:
                atom.SetIsotope(0)
            if atom.GetAtomicNum() == 8 and atom.GetIsotope() == 17:
                atom.SetIsotope(0)
            if atom.GetAtomicNum() == 8 and atom.GetIsotope() == 18:
                atom.SetIsotope(0)
        smiles = Chem.MolToSmiles(molec2)
        status = 1

    elif '.' in smiles:
        partit = smiles.split('.')
        for i in range (len(partit)):
            newSmiles = partit[i]
            if len(newSmiles) > 4:
                molec = Chem.MolFromSmiles(newSmiles)
                for bond in molec.GetBonds():
                    atom1 = bond.GetBeginAtom()

                    atom2 = bond.GetEndAtom()
                    if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() == 6 and bond.GetBondType() == Chem.BondType.DOUBLE:
                        status = 1
                        smiles = newSmiles

    elif '+' in smiles or '-' in smiles:
        status = 0
    else:
        status = 1         
    #print("76: " + str(smiles))
    return status , smiles

def decompose(smiles , CC): #CC is list of positions for C=C bonds 
    #turn any multipple double bond SMILES to more conformers with individual double bonds 
    newSMILES = []
    print("CC" , CC)
    for i in range (len(CC)):
        molec = Chem.MolFromSmiles(smiles)
        knockoutCC = CC[:i] + CC[i+1:]
        print("knockoutCC" , knockoutCC)
        for j in range (len(knockoutCC)):
            C1 = knockoutCC[j][0]
            C2 = knockoutCC[j][1]
            #connectivity = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondType())for bond in molec.GetBonds()]
            #print("connectivity" , connectivity)
            editMolec = Chem.RWMol(molec)
            double = editMolec.GetBondBetweenAtoms(C1, C2)
            double.SetBondType(Chem.BondType.SINGLE)
            molec = editMolec.GetMol()
        newSmiles = Chem.MolToSmiles(molec)
        newSMILES.append(newSmiles)
    return newSMILES
def cleanupSMILES(mainDir , savePath):
    txtFile = glob.glob(mainDir + "/*.txt")[0] #Should only be one 
    with open(txtFile, 'r') as file:
        rawSmilesList = [line.strip() for line in file]
        #print(rawSmilesList)

    if os.path.exists(savePath):
        with open(savePath , 'r') as file: 
            currentList = [line.strip() for line in file]
            currentSmile = currentList[-1]
        if currentSmile in rawSmilesList:
            j = rawSmilesList.index(currentSmile) + 1
            
        else: 
            print(currentSmile)
            sys.exit()
    else:
        j = 0

    for i in range ( j , len(rawSmilesList)):
        #print(i)
        cursedSMILES = []
        smiles = rawSmilesList[i]
        print(smiles)

        if smiles == '':
            #Ghost Smiles
            pass
        else:
            molec1 = Chem.MolFromSmiles(smiles)
            try:
                molec1.GetBonds()
                status , smiles = checkForFlaws(smiles)
                print(status)
                if status == 0:
                    pass
                else: 
                    #Smiles string is good, begin the atom mapping
                    doubleC = 0
                    CC = []
                    molec2 = Chem.MolFromSmiles(smiles)
                    molec2.GetBonds()
                    editMolec = Chem.RWMol(molec2)
                    for bond in editMolec.GetBonds():
                        at1 = bond.GetBeginAtom()
                        at2 = bond.GetEndAtom()
                        c1 = bond.GetBeginAtomIdx()
                        c2 = bond.GetEndAtomIdx()
                        if at1.GetAtomicNum() == 6 and at2.GetAtomicNum() == 6 and bond.GetBondType() == Chem.BondType.DOUBLE:
                            CC.append([c1 , c2])
                            doubleC +=1
                    if doubleC == 0:
                        continue
                    elif doubleC >= 2:
                        #More than 2 double bonds 
                        dbConfs = decompose(smiles , CC)

                        for i in range (len(dbConfs)):
                            conf1 = dbConfs[i]
                            saveSMILES(conf1 , savePath)

                    else:
                        saveSMILES(smiles, savePath)

            except AttributeError:
                #cursed Smiles save onto a seperate file 
                print("ERROR")
                cursedSMILES.append(smiles)
                continue
    return cursedSMILES

if __name__ == "__main__":
    mainDir = str(sys.argv[1]) #Gives directory that will contain all .txt files of interest: testReaxys_PostProcessing
    saveSuffix = str(sys.argv[2])

    savePath = mainDir + "/" + str(saveSuffix) + ".dat"

    cursedSMILES = cleanupSMILES(mainDir, savePath)
    #print(cursedSMILES)
    if not os.path.exists(mainDir + "/cursedSMILES.dat"):
        with open(mainDir + "/cursedSMILES.dat", "w") as file:
            for i in range (len(cursedSMILES)):
                cursed = cursedSMILES[i]
                file.write(cursed + "\n")

    createXLSX(savePath)


