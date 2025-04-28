import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from sklearn import ensemble
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_squared_error, root_mean_squared_error , make_scorer, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split , GridSearchCV, KFold , cross_val_score
from sklearn.utils.fixes import parse_version
from sklearn.ensemble import GradientBoostingRegressor
from scipy.spatial.distance import cdist
from statistics import multimode
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from DFTWorkflow.featureMaping import savePNG , createCSV
from DFTWorkflow.pitchingATent import compressData , locateNans , eliminateNans
#Danny Ruiz de Castilla 03.24.
#gives top X features responsible for random forest regression of given parameterf
def maskingColumns(featuresMAST , inputDF): 
    featureLabels = list(inputDF.columns)
    featureLabels = list(filter(None, featureLabels))  # Remove None values
    featuresMAST = list(filter(None, featuresMAST))  # Remove None values
    newFeatures = list(set(featureLabels)  & set(featuresMAST))
    newDF = inputDF[newFeatures]
    return newDF
def getOptimalParameters(parameterList):
    keyTypes = list(parameterList[0].keys())
    masterDict = {}
    for key in keyTypes:
        attributeList = []
        for dict_ in parameterList:
            attribute = dict_[str(key)]
            attributeList.append(attribute)
        listID = listTypes(attributeList)
        if not(str in listID):#if all of these are ints or floats 
            #get mathematical mode and set as dictionary parameter 
            mode1 = multimode(attributeList)
            masterDict[key] = mode1[-1]
        else:
            masterDict[key] = max(set(attributeList), key=attributeList.count)
    return masterDict
def gradientBoostingOptimization(X , Y ):
    gbrParameters = []
    gbrParmScores = []
    i = 0
    parameterDict  = { 'n_estimators': [50 , 100, 150 , 200 , 250] , 'max_depth' : [2,3,4,5,6] , 'criterion' : ['friedman_mse' , 'squared_error'] , 'loss' : ['squared_error' , 'absolute_error' , 'huber']}
    cvStrategy = KFold(n_splits=5 , shuffle = True , random_state = 42)
    gbr = GradientBoostingRegressor()
    while True:
        grid_search = GridSearchCV(gbr, parameterDict, scoring = 'neg_mean_squared_error' , cv = cvStrategy )
        grid_search.fit(X , Y)
        gbrParameters.append(grid_search.best_params_)
        gbrParmScores.append(grid_search.best_score_)
        i += 1 
        if i == 7:
            break
    optimalParameters = getOptimalParameters(gbrParameters)
    optimalParameters['random_state'] = 42

    return optimalParameters
def gradientBoostingRegression(optParms , X  , Y , saveDir , saveID , regressionStr):
#{'criterion': 'friedman_mse', 'loss': 'squared_error', 'max_depth': 2, 'n_estimators': 100}
    gbrOpt = GradientBoostingRegressor(**optParms)
    cvStrategy = KFold(n_splits=5 , shuffle = True , random_state = 42)
    X_train_CV, X_test_CV, y_train_CV, y_test_CV = train_test_split(X, Y, test_size=0.2, random_state=42)
    cvScores = cross_val_score(gbrOpt , X_train_CV , y_train_CV , cv = cvStrategy , scoring = 'r2')


    gbrOpt.fit(X_train_CV , y_train_CV)
    yPred = gbrOpt.predict(X_test_CV)
    r2Test = r2_score(y_test_CV , yPred)
    

    featureImportance = gbrOpt.feature_importances_
    sortedIdxs = np.argsort(featureImportance)[::-1]
    if not os.path.exists(saveDir + "regressionData/"):
        os.makedirs(saveDir + "regressionData/")
    with open(saveDir + "regressionData/"+ "r2Test" + str(saveID) + ".dat" , "w") as file:
        file.write(f"R^2 for {regressionStr} Test Set: {r2Test:.4f}\n")
    
    with open(saveDir + "regressionData/"+ "featureImportance" + str(saveID) + ".dat" , "w") as file:
        file.write("Feature,Importance\n")
        for idx in sortedIdxs:
            file.write(f"{X.columns[idx]},{featureImportance[idx]:.4f}\n")
class listTypes(list):
    def __contains__(self, typ):
        inst = len(self)
        count = 0
        for val in self:
            if isinstance(val , typ):
                count += 1 
        if count == inst:
            return True 
        elif count == 0:
            return False
        else:
            print(f"Error: Mixed type file {self}")
            sys.exit()
if __name__ == "__main__":

    datasetDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    regressionStr = str(sys.argv[3])#Usually will be "Yield"
    featureDataSets = glob.glob(datasetDir + "*.csv")
    dfMASTDir = str(sys.argv[4])
    initdataSets = sorted(featureDataSets)
    Xdataframe , smileList  , yieldList = compressData(initdataSets , regressionStr)
    nanDict = locateNans(Xdataframe)
    if len(nanDict) != 0:
        Xdataframe[regressionStr] = yieldList
        Xdataframe = eliminateNans(Xdataframe , nanDict)

    yields = Xdataframe[regressionStr].copy()
    Xdataframe = Xdataframe.drop(regressionStr, axis=1)
    mastFeatureSet = glob.glob(dfMASTDir + "*.csv")[0]   
    dfMAST = pd.read_csv(mastFeatureSet)
    featuresMAST = list(dfMAST.columns)
    X = maskingColumns(featuresMAST, Xdataframe)


    gbrOptParms = gradientBoostingOptimization(X , yields )
    k = 0
    while True:
        gradientBoostingRegression(gbrOptParms , X  , yields , outputDir , str(k), regressionStr)
        k += 1

        if k == 10: 
            break