import pandas as pd
import glob
import numpy as np
import os 
import sys
import json

#Sorts substrates into methods, assigning dicts to substrates that correspond to multiple methods and multiple Y values (Yield , %ee , etc)
#All .csv files must have consistent dataframe labeling 
def smilesToCanonical(smilesList):

    try:
        from rdkit import Chem
    except ImportError:
        raise ImportError("RDKit is required for this function. Please install it with 'pip install rdkit'")
    
    canonicalS = []
    
    for smiles in smilesList:
        try:
            # Convert SMILES to RDKit molecule
            mol = Chem.MolFromSmiles(smiles)
            
            if mol is None:
                # If conversion fails, add None to the results
                canonicalS.append(None)
                print(f"Warning: Could not parse SMILES: {smiles}")
            else:
                # Convert molecule to canonical SMILES
                canonical = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
                canonicalS.append(canonical)
        except Exception as e:
            # Handle any unexpected errors
            print(f"Error processing SMILES '{smiles}': {str(e)}")
            canonicalS.append(None)
    
    return canonicalS

def classifySMILES(chemistryLists):

    poolAssignments = {}
    
    for i, chemistry in enumerate(chemistryLists):
        list_value = 2**i #unique binary number based on chemistry
        
        for smiles in chemistry:
            if smiles in poolAssignments:
                # If we've seen this SMILES before
                poolAssignments[smiles] += list_value
            else:
                # First time seeing this SMILES
                poolAssignments[smiles] = list_value
    
    return poolAssignments

def createPool(mainDir, smileStr ):
    canonicalMain = []
    for dir in mainDir:
        df = pd.read_csv(dir)
        smilesList = list(df[smileStr])
        canonicalList = smilesToCanonical(smilesList)
        canonicalMain.append(canonicalList)

    return canonicalMain
if __name__ == "__main__":
    substrateDir = str(sys.argv[1])
    smileString = str(sys.argv[2])
    saveDir = str(sys.argv[3])

    chemistryDirs = glob.glob(substrateDir + "/*.csv" )
    canonicalChemistries = createPool(chemistryDirs, smileString)
    chemistriesDict = classifySMILES(canonicalChemistries)
    
    with open(saveDir + 'substratePool.json', 'w') as f:
        json.dump(chemistriesDict, f, indent=4) 


