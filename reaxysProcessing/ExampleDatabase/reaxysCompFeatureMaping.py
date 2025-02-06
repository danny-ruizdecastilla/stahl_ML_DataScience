import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
from rdkit.Chem.PandasTools import LoadSDF
import matplotlib.pyplot as plt
import random
from itertools import combinations
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from rdkit.Chem import Descriptors, MACCSkeys
from rdkit.Chem import AllChem
import datetime
from expt2_create_features import *
from expt2_feature_filtering import *
#Danny Ruiz de Castilla 02.01.2025

def createCSV(dataFrame , saveDir, saveStr: str):
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    dataFrame.to_csv(saveDir + "/" + saveStr + ".csv", index=True) 
    return dataFrame

def DimensionalityReduction(X, Y , labels, partition1 ):


    return dataframeMAST


if __name__ == "__main__":
#input an xlsx dataframe and depending on what features you want to calculate, export .csv of features, performs feature elimination, maps features, and creates clusters 
    mainDir = str(sys.argv[1])
    rdkitAllow = int(sys.argv[2])
    moldrAllow = int(sys.argv[3])
    maccAllow = int(sys.argv[4])
    
    featureTypes = []
    if rdkitAllow == 1:
        rdkitSave = "rdkit"
        featureTypes.append(rdkitSave)
    if moldrAllow == 1:
        moldrSave = "moldred"
        featureTypes.append(moldrSave)
    if maccAllow == 1:
        maccSave = "MACCSkeys"
        featureTypes.append(maccSave)
    #defines options for features you want 

    
    


    

    
    





    
