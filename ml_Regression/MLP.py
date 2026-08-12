from sklearn.neural_network import MLPRegressor
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split,  KFold , cross_val_score
from sklearn.metrics import r2_score
from pathlib import Path
import sys
import pandas as pd 
import numpy as np
import shap
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
def rmseCalc(fit , actual):
    sum_ = np.sum((np.array(actual) - np.array(fit))**2)
    n = len(fit)
    sumSqr = sum_/n
    return np.sqrt(sumSqr)
def main( featDF , yStr , idStr , smileStr , outputDir , numSeeds):
    yVals = featDF[yStr]
    idVals = featDF[idStr]    
    smilesList = featDF[smilesStr]
    dfMAST = featDF.drop(columns = [yStr, idStr , smileStr])
    scaler = StandardScaler()

    dfMAST = pd.DataFrame(
        scaler.fit_transform(dfMAST),
        columns=dfMAST.columns,
        index=dfMAST.index
    )
    kfold = KFold(
        n_splits=5,
        shuffle=True,
        random_state=42
    )
    cvVarList = []
    meanCVList = []
    meanTestRMSEList = []
    testR2List = []
    for i in range(numSeeds):
        mlp = MLPRegressor(activation = 'relu' , alpha =  0.01, hidden_layer_sizes = (100, 50), learning_rate_init =  0.005)

        xTrain , xTest , yTrain , yTest = train_test_split(dfMAST, yVals, test_size=0.25, random_state=i)
        mlpCVRMSE = np.sqrt(-cross_val_score(mlp, xTrain, yTrain,cv=kfold,scoring="neg_mean_squared_error"))
        cvVar = np.var(np.array(mlpCVRMSE))
        cvVarList.append(cvVar)
        meanCV = np.mean(mlpCVRMSE)
        meanCVList.append(meanCV)
        mlp.fit(xTrain, yTrain)
        yPred = mlp.predict(xTest)
        testRMSE = rmseCalc(yPred, yTest)
        meanTestRMSEList.append(testRMSE)
        testR2 = r2_score(yTest , yPred)
        testR2List.append(testR2)
        '''
        shapExplainer = shap.KernelExplainer(mlp.predict_proba, xTrain)
        shapValues = shapExplainer.shap_values(xTest)
        shapDF = pd.DataFrame(shapValues,index=dfMAST.index,columns=dfMAST.columns)
        featureColsHash = {}
        for idx , feature in enumerate(shapDF.columns):
            shapVal = np.mean([np.abs(val) for val in list(shapDF[feature])])
            featureColsHash[feature] = shapVal
        sortedCols = sorted(featureColsHash.items(),key=lambda kv: kv[1])
        '''
    testRMSEVar = np.var(np.array(meanTestRMSEList))
    cvRMSEVar = np.var(np.array(meanCVList))
    print("Variance for all test splits: ", testRMSEVar)
    print("Variance for all avg. CV splits: ", cvRMSEVar)
    print("Mean R2: ", np.mean(testR2List))
    testPath = Path(outputDir) / "MLPmeanTest.dat"
    cvPath = Path(outputDir) / "MLPmeanCV.dat"
    cvVarPath = Path(outputDir) / "MLPcvVar.dat"
    testPath.write_text(",".join(str(item) for item in meanTestRMSEList))
    cvPath.write_text(",".join(str(item) for item in meanCVList))
    cvVarPath.write_text(",".join(str(item) for item in cvVarList))
    
        

if __name__ == "__main__":
    featureDFStr = str(sys.argv[1])
    yStr = str(sys.argv[2])
    idStr = str(sys.argv[3])
    smilesStr = str(sys.argv[4])
    outputDir = str(sys.argv[5])
    featureDF = pd.read_csv(featureDFStr)
    main(featureDF , yStr , idStr , smilesStr , outputDir , 200 )