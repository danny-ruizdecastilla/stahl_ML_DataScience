import pandas as pd
import os 
import sys
import glob
import re
import chemdraw
import base64
import numpy as np
from pathlib import Path
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from dimensionalityReduction.reactivityFeatures import boxGen
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs

#Danny Ruiz de Castilla 05.02.25
def main(df , saveDir):
    dfMAST = pd.read_csv(df)
    cols = list(dfMAST.columns)
    boxCols = boxGen(cols)
    dropIdx = listInputs(f"Input the indices of all the columns you want to (temporarily) remove before doing a PCA\n{boxCols}\n")
    colsHash = {idx : dfMAST[cols[idx]] for idx in dropIdx}
    dropCols = [cols[idx] for idx in dropIdx]
    X = dfMAST.drop(columns=dropCols)
    np.random.seed(42)
    scaler = StandardScaler()
    scaledX = scaler.fit_transform(X)
    pca = PCA(n_components = X.shape[1] ,svd_solver="full" )    
    explainedVar = pca.explained_variance_ratio_
    if not os.path.exists(str(saveDir) + "/" + "explainedVarr.dat"):
        with open(str(saveDir) + "/" + "explainedVarr.dat", "w") as file:
            for i in range (len(explainedVar)):
                file.write(f"explained variance ratio: PC {i + 1} {explainedVar[i]:.6f}\n")

    top2 = np.argsort(explainedVar)[-2:][::-1] 
    #print(top2)
    top2PCA = pca.components_[top2]
    xPCA_ = scaledX @ top2PCA.T + 11
    #print(xPCA_)
    loadings = pd.DataFrame(pca.components_.T * np.sqrt(pca.explained_variance_), columns=[f'PC{i+1}' for i in range(pca.n_components_)] , index = X.columns)

    top_features_pc1 = loadings['PC1'].abs().sort_values(ascending=False)
    top_features_pc2 = loadings['PC2'].abs().sort_values(ascending=False)

    topFeatures = pd.DataFrame({
        "Feature": top_features_pc1.index,
        "PC1_Contribution": top_features_pc1.values,
        "PC2_Contribution": top_features_pc2.loc[top_features_pc1.index].values  # Align PC2 with PC1 sorting
    })

    filePath = "topFeatures.csv"
    if not os.path.exists(str(saveDir) + "/" + filePath):
        topFeatures.to_csv(str(saveDir) + "/" + filePath, sep="\t", index=False)
    dfPCA = pd.DataFrame(xPCA_, columns = ["PCA1" , "PCA2"] )
    pc1 = list(dfPCA["PCA1"]) 
    pc2 = list(dfPCA["PCA2"])
    dfMAST["pc1"] = pc1
    dfMAST["pc2"] = pc2
    for col , info in colsHash.items():
        dfMAST[col] = info
    
    dfMAST.to_csv(str(saveDir) + "/pcaDataframe.csv")


if __name__ == "__main__":
    dfMain = str(sys.argv[1])
    saveStr = str(sys.argv[2])
    dfPath = Path(dfMain)
    saveDir = dfPath.parent / f"pca_{saveStr}" 
    saveDir.mkdir(parents=True, exist_ok=True)
    main(dfMain , saveDir)

    