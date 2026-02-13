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
from dimensionalityReduction.reactivityFeatures import boxGen

def generate_UMAP_tSNE(X):
    colList = list(X.columns)
    dfColBox = boxGen(colList)
    elimCols = listInputs(f"Here are all available columns for the dataframe {dfColBox}\n please enter the matching indexes corresponding to the columns you want to remove from the dataframe:\n")
    dropFeatures = []
    for idx in elimCols:
        column = colList[int(idx)]
        dropFeatures.append(column)

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
    perplexSilScoreMax = max(tSNE_Hash, key=tSNE_Hash.get)

    tsne = TSNE(n_components=2, perplexity=int(perplexSilScoreMax), random_state=0)
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

def main(df:np.array , outputDir:str , featureType:str):
    from pathlib import Path
    from dimensionalityReduction.kMeansClustering import kMeansCluster
    UMAP_1, UMAP_2,tSNE_1,tSNE_2 = generate_UMAP_tSNE(df)
    df[f"{featureType}_UMAP_1"] = UMAP_1
    df[f"{featureType}_UMAP_2"] = UMAP_2
    df[f"{featureType}_tSNE_1"] = tSNE_1
    df[f"{featureType}_tSNE_2"] = tSNE_2
    outputDir_tSNE = Path(outputDir) / "tSNE"
    outputDir_tSNE.mkdir(parents=True, exist_ok=True)
    dftSNE = kMeansCluster(df, f"{featureType}_tSNE_1", f"{featureType}_tSNE_1", 15, 45, "tSNE" , outputDir_tSNE)
    dftSNE.to_csv(outputDir_tSNE  / f"clustered_{featureType}.csv")
    outputDir_UMAP = Path(outputDir) / "UMAP"
    outputDir_UMAP.mkdir(parents=True, exist_ok=True)
    dftSNE = kMeansCluster(df, f"{featureType}_UMAP_1", f"{featureType}_UMAP_1", 15, 45 , "UMAP" , outputDir_tSNE)
    dftSNE.to_csv(outputDir_tSNE  / f"clustered_{featureType}.csv")
if __name__ == "__main__":
    datasetDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    featureType =  str(sys.argv[3])
    dataset = pd.read_csv(datasetDir)
    main(dataset, outputDir , featureType)