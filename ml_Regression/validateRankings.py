#uses a fixed training sample to predict orderings of a test sample and select the values closed to the fixed prediction to evaluate an ordering
import os 
import sys
import glob
import pandas as pd 
import numpy as np
import chemdraw
from pathlib import Path
from rdkit import Chem
import chemdraw
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
def main():


if __name__ == "__main__":
    trainingCSV = str(sys.argv[1])
    testingCSV = str(sys.argv[2])
    binningNum = int(sys.argv[3])
    outputDir = str(sys.argv[4])
    trainingDF = pd.read_csv(trainingCSV)
    testingDF = pd.read_csv(testingCSV)
    main(trainingDF  , testingDF,binningNum  , outputDir  )