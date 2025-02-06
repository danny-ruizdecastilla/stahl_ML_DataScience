import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
from rdkit.Chem.PandasTools import LoadSDF
import matplotlib.pyplot as plt
import random
import networkx as nx
from networkx import Graph
from itertools import combinations
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from rdkit.Chem import Descriptors, MACCSkeys
from rdkit.Chem import AllChem
import datetime

from expt2_create_features import *
def createXLSX(dataFrame , saveDir, saveStr: str):
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)

    dataFrame.to_excel(saveDir + "/" + saveStr + ".xlsx", index=False) 
'''
def reducedDF(features, labels):
    X = features.to_numpy()
    from expt2_feature_filtering import *
    XFiltered , removedFeatures = spearmanr_correlation(X: np.ndarray, labels: list[str], threshold = 0.98, output  = False)
    return XFiltered , removedFeatures 

'''

if __name__ == "__main__":
    fileName = str(sys.argv[1])
    smilesStr = str(sys.argv[2])
    yieldStr = str(sys.argv[3])
    idStr = str(sys.argv[4])
    featureDir = str(sys.argv[5])
    smilesDF = pd.read_excel(fileName)

    smilesMAST = list(smilesDF[smilesStr])
    idMAST = list(smilesDF[idStr])

    rdkitFeatures = calc_rdkit_desc(smilesMAST, multiprocess_ = True)
    rdKitLabels = rdkitFeatures.columns.tolist()

    yieldsMAST = list(smilesDF[yieldStr])

    #rdkitReduced = 
    
    print(rdkitFeatures.shape)
    createXLSX(rdkitFeatures , featureDir , "rdkit")
    maccFeats = calc_maccs_keys(smilesMAST, multiprocess_ = True)
    createXLSX(maccFeats, featureDir , "macc")
    '''
    morganFingerprints = calc_morgan_keys(smilesMAST , multiprocess_ = True)
    createXLSX(morganFingerprints , "features/" , "morgan")
    '''

    

    
    





    
