#Danny Ruiz de Castilla 06.22.2025
import sys
import glob
import numpy as np
import os
import pandas as pd
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(parentDir)
from DFTWorkflow.alkeneSubstitution import  smiles_to_coords
def write_xyz(elements, coords, xyzPath):
    numAtoms = len(elements)
    with open(xyzPath, "w") as f:
        f.write(f"{numAtoms}\n\n")
        for i in range (len(elements)):
            f.write(f"{elements[i]}  {coords[i][0]:.7f}  {coords[i][1]:.7f}  {coords[i][2]:.7f}\n")

def main(csvDir , outputDir):
    df = pd.read_csv(csvDir)
    for _, row in df.iterrows():
        smiles = row["SMILES"]
        id = row["ID"]
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