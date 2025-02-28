import sys
#import glob
import numpy as np
import os
#from rdkit import Chem
import pandas as pd
#from rdkit.Chem.PandasTools import LoadSDF
import matplotlib.pyplot as plt
import random
import hdbscan
from itertools import combinations
from sklearn.cluster import KMeans, DBSCAN, AgglomerativeClustering
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
#from rdkit.Chem import Descriptors, MACCSkeys
#from rdkit.Chem import AllChem
from hdbscan import HDBSCAN
import umap
from expt2_feature_filtering import *
#Danny Ruiz de Castilla 02.28.2025

#pcaDict = {canonicalSmiles : [ xPCA , yPCA , yield]} for local Chemistries 
#smilesDict = {canonicalSmiles : [ xPCA , yPCA]}
def theTent(smilesDict  ,  )







class dfTransformer:
    def __init__(self, dataDir):
        """Initialize the class with a DataFrame from an Excel file."""
        self.df = pd.read_excel(dataDir)