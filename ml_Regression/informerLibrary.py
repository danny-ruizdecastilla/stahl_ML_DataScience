import pandas as pd
import os
import sys
import numpy as np
import chemdraw
from pathlib import Path
import random

parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))

#Danny Ruiz de Castilla
def getInformerLibraryFromCluster(dfMAST, smilesStr, informerInt,clusterStr , outputDir ):
    clusters = list(dfMAST[clusterStr])
    clusterCount = len(set(clusters))
    eachCluster = clusterCount // informerInt
    remain = clusterCount % informerInt

    informerLibrary = {}
    popType = {}
    for cluster in list(set(clusters)):
        df = dfMAST[dfMAST[clusterStr] == cluster]
        selectRows = df.sample(n=eachCluster)
        smilesList = []
        for row in selectRows.iterrows():
            smiles = row[smilesStr]
            smiles.append(smilesList)
        
        informerLibrary[cluster] = smilesList
    popType[cluster] = int(len(df))

    if remain != 0:
        sortedPop= dict(sorted(popType.items(), key=lambda x: x[1], reverse=True))
        count = 0
        for key, val in sortedPop.items():
            if count == remain:
                break
            df = dfMAST[dfMAST[clusterStr] == int(key)]
            smilesList = informerLibrary[key]
            selectSubstrates = list(df[smilesStr])
            while True:
                smiles = random.choice(selectSubstrates)
                if smiles not in smilesList:
                    smilesList.append(smiles)
                
                    break
                else:
                    continue
            informerLibrary[key] = smilesList

            count += 1
    smilesListMAST = []
    keyListMAST = []
    for key, val in informerLibrary.items():
        for smiles in val:
            smilesListMAST.append(smiles)
            keyListMAST.apppend("cluster: " + str(key))
    drawer = chemdraw.GridDrawer(smilesListMAST , title = keyListMAST)
    drawer.draw_png(outputDir + "/informerLibraray.png")
if __name__ == "__main__":
    dfMASTDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    smilesStr = "SMILES"
    informerInt = int(sys.argv[3])
    clusterStr = "Cluster"
    dfMAST = pd.read_csv(dfMASTDir)
    getInformerLibraryFromCluster(dfMAST, smilesStr, informerInt,clusterStr , outputDir )