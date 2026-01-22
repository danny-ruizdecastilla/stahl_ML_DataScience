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
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def main(df, saveDir):
    df = Path(df)
    saveDir = Path(saveDir)
    dfMAST = pd.read_csv(df)

    cols = list(dfMAST.columns)
    boxCols = boxGen(cols)

    dropIdx = listInputs(
        "Input the indices of all the columns you want to (temporarily) "
        "remove before doing a PCA\n"
        f"{boxCols}\n"
    )
    dropIdx = [int(idx) for idx in dropIdx]

    colsHash = {idx: dfMAST[cols[idx]] for idx in dropIdx}
    dropCols = [cols[idx] for idx in dropIdx]

    X = dfMAST.drop(columns=dropCols)

    np.random.seed(42)

    scaler = StandardScaler()
    scaledX = scaler.fit_transform(X)

    pca = PCA(n_components=X.shape[1], svd_solver="full")
    xPCAFullRank = pca.fit_transform(scaledX)

    explainedVar = pca.explained_variance_ratio_

    explainedVarFile = saveDir / "explainedVarr.dat"
    if not explainedVarFile.exists():
        with explainedVarFile.open("w") as f:
            for i, var in enumerate(explainedVar, start=1):
                f.write(f"explained variance ratio: PC {i} {var:.6f}\n")

    top2 = np.argsort(explainedVar)[-2:][::-1]
    top2PCA = pca.components_[top2]

    xPCA_ = scaledX @ top2PCA.T + 11

    loadings = pd.DataFrame(
        pca.components_.T * np.sqrt(pca.explained_variance_),
        columns=[f"PC{i+1}" for i in range(pca.n_components_)],
        index=X.columns
    )

    top_features_pc1 = loadings["PC1"].abs().sort_values(ascending=False)
    top_features_pc2 = loadings["PC2"].abs().sort_values(ascending=False)

    topFeatures = pd.DataFrame({
        "Feature": top_features_pc1.index,
        "PC1_Contribution": top_features_pc1.values,
        "PC2_Contribution": top_features_pc2.loc[top_features_pc1.index].values
    })

    topFeaturesFile = saveDir / "topFeatures.csv"
    if not topFeaturesFile.exists():
        topFeatures.to_csv(topFeaturesFile, sep="\t", index=False)

    dfMAST["pc1"] = xPCA_[:, 0]
    dfMAST["pc2"] = xPCA_[:, 1]

    for col, info in colsHash.items():
        dfMAST[col] = info

    dfMAST.to_csv(saveDir / "pcaDataframe.csv", index=False)

if __name__ == "__main__":
    dfMain = str(sys.argv[1])
    saveStr = str(sys.argv[2])
    dfPath = Path(dfMain)
    dfPath = Path(dfMain).resolve()
    saveDir = dfPath.parent / f"pca_{saveStr}"
    saveDir.mkdir(parents=True, exist_ok=True)

    main(dfMain, saveDir)

    