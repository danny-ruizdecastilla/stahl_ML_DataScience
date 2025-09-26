# Danny Ruiz de Castilla 06.22.2025
import sys
import glob
import numpy as np
import os
import pandas as pd
from pathlib import Path

parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))

from DFTWorkflow.AlkeneFeatures.alkeneSubstitution import smiles_to_coords
from dimensionalityReduction.reactivityFeatures import boxGen


def write_xyz(elements, coords, xyzPath):
    numAtoms = len(elements)
    with open(xyzPath, "w") as f:
        f.write(f"{numAtoms}\n\n")
        for i in range(len(elements)):
            f.write(f"{elements[i]}  {coords[i][0]:.7f}  {coords[i][1]:.7f}  {coords[i][2]:.7f}\n")


def main(csvDir, outputDir):
    df = pd.read_csv(csvDir)
    colList = list(df.columns)
    colBox = boxGen(colList)

    while True:
        smileIdx = input(
            f"Here are the columns provided in the dataframe: {csvDir} "
            f"Enter the number for the column corresponding to the SMILES strings {colBox}"
        )
        try:
            smilesIdx = int(smileIdx)
            smilesCol = colList[smilesIdx]
            break
        except:
            print("Error: please enter an appropriate integer")

    while True:
        idIdx = input(
            f"Here are the columns provided in the dataframe: {csvDir} "
            f"Enter the number for the column corresponding to the ID strings {colBox}"
        )
        try:
            idIdx = int(idIdx)
            idCol = colList[idIdx]
            break
        except:
            print("Error: please enter an appropriate integer")

    for _, row in df.iterrows():
        smiles = row[smilesCol]
        id = row[idCol]
        elements, coords = smiles_to_coords(smiles)

        xyzDir = Path(outputDir) / str(id)
        xyzFile = xyzDir / f"{id}.xyz"

        xyzDir.mkdir(parents=True, exist_ok=True)
        write_xyz(elements, coords, xyzFile)


if __name__ == "__main__":
    smilesCSV = Path(sys.argv[1])
    outputDir = Path(sys.argv[2])

    outputDir.mkdir(parents=True, exist_ok=True)
    main(smilesCSV, outputDir)
