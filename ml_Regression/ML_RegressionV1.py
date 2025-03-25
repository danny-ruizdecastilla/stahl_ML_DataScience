import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
from rdkit.Chem.PandasTools import LoadSDF
import matplotlib.pyplot as plt
import random
import hdbscan
from itertools import combinations
from sklearn.ensemble import RandomForestRegressor
#from rdkit.Chem import Descriptors, MACCSkeys
from rdkit.Chem import AllChem
from hdbscan import HDBSCAN
from Code.stahl_ML-DataScience.DFTWorkflow.pitchingATent import transformations
import umap

#Danny Ruiz de Castilla 03.24.
#gives top X features responsible for random forest regression of given parameter


if __name__ == "__main__":

    datasetDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    regressionStr = str(sys.argv[3])#Usually will be "Yield"
    featureDataSets = glob.glob(datasetDir + "*.csv")
    print(featureDataSets)