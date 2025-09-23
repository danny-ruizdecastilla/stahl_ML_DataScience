#Danny Ruiz de Castilla 06.22.2025
import sys
import glob
import numpy as np
import os
import pandas as pd
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(parentDir)
from DFTWorkflow.AlkeneFeatures.alkeneSubstitution import  smiles_to_coords
from dimensionalityReduction.reactivityFeatures import boxGen
def write_xyz(elements, coords, xyzPath):
    numAtoms = len(elements)
    with open(xyzPath, "w") as f:
        f.write(f"{numAtoms}\n\n")
        for i in range (len(elements)):
            f.write(f"{elements[i]}  {coords[i][0]:.7f}  {coords[i][1]:.7f}  {coords[i][2]:.7f}\n")

def main(csvDir , outputDir):
    df = pd.read_csv(csvDir)
    colList = list(df.columns)
    colBox = boxGen(colList)
    while True:
        smileIdx = input(f"Here are the columns provided in the dataframe: {csvDir} Enter the number for the column corresponding to the SMILES strings {colBox}")
        try:
            smilesIdx = int(smileIdx)
            smilesCol = colList[smilesIdx]
            break
        except:
            print("Error: please enter an appropriate integer")
    while True:
        idIdx = input(f"Here are the columns provided in the dataframe: {csvDir} Enter the number for the column corresponding to the ID strings {colBox}")
        try:
            idIdx = int(idIdx)
            idCol = colList[idIdx]
            break
        except:
            print("Error: please enter an appropriate integer")
    for _, row in df.iterrows():
        smiles = row[smilesCol]
        id = row[idCol]
        elements , coords = smiles_to_coords(smiles)
        xyzDir = outputDir + "/" + str(id)
        xyzFile = xyzDir + "/" + str(id) + ".xyz"
        if not os.path.exists(xyzDir):
            os.makedirs(xyzDir)
        write_xyz(elements, coords, xyzFile)


if __name__ == "__main__":
    smilesCSV = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    main(smilesCSV , outputDir)