import pandas as pd
import sys
import random
import chemdraw
#Danny Ruiz de Castilla

def getInformerLibraryFromCluster(dfMAST, smilesStr, informerInt, clusterStr, outputDir):
    clusters = list(dfMAST[clusterStr])
    clusterCount = len(set(clusters))
    eachCluster = int(informerInt/clusterCount)
    print(eachCluster)
    remain = informerInt % clusterCount 
    informerLibrary = {}
    popType = {}
    
    for cluster in list(set(clusters)):
        df = dfMAST[dfMAST[clusterStr] == cluster]
        selectRows = df.sample(n=eachCluster)
        smilesList = []
        for idx, row in selectRows.iterrows():
            smiles = row[smilesStr]
            smilesList.append(smiles)
        informerLibrary[cluster] = smilesList
        popType[cluster] = int(len(df))
    
    if remain != 0:
        sortedPop = dict(sorted(popType.items(), key=lambda x: x[1], reverse=True))
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
            keyListMAST.append("cluster: " + str(key))
    
    print(keyListMAST)
    drawer = chemdraw.GridDrawer(smilesListMAST)
    drawer.draw_png(outputDir + "/informerLibrary.png")

if __name__ == "__main__":
    dfMASTDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    smilesStr = "SMILES"
    informerInt = int(sys.argv[3])
    clusterStr = "Cluster"
    dfMAST = pd.read_csv(dfMASTDir)
    getInformerLibraryFromCluster(dfMAST, smilesStr, informerInt, clusterStr, outputDir)