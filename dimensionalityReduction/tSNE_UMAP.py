#Modular tSNE and UMAP dimensionality reduction code
import numpy as np
import pandas as pd
import umap 
import os 
import sys
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs

def generate_UMAP_tSNE(X):
    motifList = listInputs(f"Enter motif strings for features you do NOT want to include in the tSNE and UMAP:]\n")
    dropFeatures = [col for col in X.columns if any(motif in col for motif in motifList)]

    X = X.drop(columns=dropFeatures)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    tSNE_Hash = {}
    for perp in [5, 7,10,12,15,17,20]:
        tsne = TSNE(n_components=2, perplexity=perp, random_state=42)
        X_tsne = tsne.fit_transform(X_scaled)
        kmeans = KMeans(n_clusters=15, random_state=42)
        labels = kmeans.fit_predict(X_tsne)
        
        sil = silhouette_score(X_tsne, labels)
        tSNE_Hash[perp] = sil
    silScoreMin = max(list(tSNE_Hash.items()))
    for key , val in tSNE_Hash.items():
        if val == silScoreMin:
            tsne = TSNE(n_components=2, perplexity=int(key), random_state=0)
            tSNE_1 = X_tsne[:, 0] 
            tSNE_2 = X_tsne[:, 1]
    UMAP_Hash ={}
    for n_neighbors in [30,50,70,90,110,130,150]:
        for min_dist in [0.0, 0.1, 0.25 , 0.5]:
            reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, random_state=42)
            X_umap = reducer.fit_transform(X_scaled)
            
            kmeans = KMeans(n_clusters=15, random_state=42)
            labels = kmeans.fit_predict(X_umap)
            
            sil = silhouette_score(X_umap, labels)
            UMAP_Hash[sil] = [n_neighbors , min_dist]
    silScoreUMAP = max(list(UMAP_Hash.keys()))  
    values = UMAP_Hash[silScoreUMAP]
    reducer = umap.UMAP(n_neighbors=values[0], min_dist=values[1], random_state=42)
    X_umap = reducer.fit_transform(X_scaled)
    UMAP_1 = X_umap[:, 0]
    UMAP_2 = X_umap[:, 1]
    return UMAP_1, UMAP_2,tSNE_1,tSNE_2