#MLR workflow for combination of features
import pandas as pd
import numpy as np
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn import preprocessing
import sys
from pathlib import Path
from sklearn.preprocessing import StandardScaler
def combinations(pool, r):
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))

    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)
def rmseCalc(fit , actual):
    sum_ = np.sum((np.array(actual) - np.array(fit))**2)
    n = len(fit)
    sumSqr = sum_/n
    return np.sqrt(sumSqr)
def maeCalc(fit , actual):
    return np.sum(np.abs(np.array(actual) - np.array(fit)))/len(fit)
def MLR(df , yCol , paramNum):
    resultsHash = {}
    colList = tuple(df.columns)
    colsCombos = list(combinations(colList, paramNum))
    for cols in colsCombos:
        X = df[list(cols)]
        scaler = StandardScaler()
        X = scaler.fit_transform(X)
        rmseDump = []
        maeDump = []
        rmseTrain = []
        for i in range(5):
            mlr = LinearRegression()
            X_train, X_test, y_train, y_test = train_test_split(X, yCol, test_size=0.3, random_state=i)
            mlr.fit(X_train,y_train)
            y_fit = mlr.predict(X_test)
            y_fitT = mlr.predict(X_train)
            rmseTrain.append(rmseCalc(y_fitT , y_train))
            rmse = rmseCalc(y_fit , y_test)
            mae = maeCalc(y_fit , y_test)
            rmseDump.append(rmse)
            maeDump.append(mae)
        dataHash = {"rmseDump" : rmseDump , "maeDump" : maeDump , "rmseAvg" : np.mean(rmseDump) , "rmseTrainAvg" : np.mean(rmseTrain)}
        resultStr = ""
        for col in cols:
            resultStr += f"{col}&"
        resultsHash[resultStr] = dataHash
    
    sortedResults = dict(
        sorted(
            resultsHash.items(),
            key=lambda kv: kv[1]["rmseAvg"]
        )
    )
    return sortedResults
if __name__ == "__main__":
    dfStr = str(sys.argv[1])
    numParms = int(sys.argv[2])
    outStr = str(sys.argv[3])
    dfPath = Path(dfStr)
    outName = dfPath.stem
    outPath = Path(outStr)
    outPath.mkdir(parents=True, exist_ok=True)
    dfMain = pd.read_csv(dfStr)
    #yCol = dfMain["RelRate"]
    yCol2 = dfMain["k2"]
    yCol3 = dfMain["err"]
    dfMain = dfMain.drop(columns=[ "k2" , "err" , "SMILES" , "Canonicals" , "ID"])
    MLRHash = MLR(dfMain , yCol2 , numParms)
    MLRFile = outPath / f"{outName}_MLR.dat"
    with MLRFile.open("w") as f:
        for key , val in MLRHash.items():
            f.write(f"{key} : {val}\n\n")