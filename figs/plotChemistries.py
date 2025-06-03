import pandas as pd
import os 
import sys
import glob
import json
import numpy as np
import html
import plotly
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
#from figs.chemPlotlyV2 import createPNGDF,png64
#from figs.stericvselectroPCA import pcafeatureSplitter
from figs.plotSubstrates import safeStringHTML , plotSubstratesMain
#from DFTWorkflow.pitchingATent import compressData , locateNans , eliminateNans , convertCanonical , featureFiltering
#from DFTWorkflow.featureMaping import  createCSV
def checkSubstratePath(substrateDF):
    if not os.path.exists(substrateDF):
        print(f"\nWarning: Could not locate master substrate dataset at:\n   → {substrateDF}")
        print("It looks like it doesn't exist yet.")
        print("Would you like to build it from scratch?")
        while True:
            userInput = input("➡️  Enter 1 to build from scratch, or 2 to abort: ").strip()
            if userInput == "1":
                print("Proceeding to build the dataset from scratch...\n")
                return 1
            elif userInput == "2":
                print("Aborting the process\n")
                return 2
            else:
                print("⚠️  Invalid input. Please enter 1 (yes) or 2 (no).")
    else:
        print(f"Found existing master dataset at: {substrateDF}")
        return 0
def createAxisMotifs(axisNum):
    axisDict = {}
    for num in range(1,axisNum+1):
        motifList = listInputs(f"Enter motif names for axis {num}, Ex: [distance,Buried,angle,dihedral,Vbur]: ")
        naming = input(f"Do you want to name this axis? | Enter 1 for yes, or 2 to accept {num} as the name for axis  {num}:")
        while True:
            if naming == "1":
                name = input(f"Enter the name for axis {num}:")
                break
            elif naming == "2":
                name = str(num)
                break
            else:
                print("⚠️  Invalid input. Please enter 1 (yes) or 2 (no).")
        axisDict[name] = motifList
    return axisDict
def main(chemistryFiles, masterDF, chemistriesDict):

    


if __name__ == "__main__":
    chemistriesDir = str(sys.argv[1])
    substrateFile = str(sys.argv[2])
    initDataSets = glob.glob(chemistriesDir + "/*.csv")
    initdataSets = sorted(initDataSets)
    chemistryNames = [name.split("/")[-1].split(".")[0] for name in initdataSets]
    chemistriesDict = {}
    for name in chemistryNames:
        nameList = listInputs(f"For Chemistry: {name} Please type the name you want this chemistry to be represented by: ")
        name_ = nameList[0]
        chemistriesDict[name] = name_
    #print(chemistriesDict)
    dataOption = checkSubstratePath(substrateFile)
    if dataOption == 2:
        sys.close()
    elif dataOption == 1:
        chemistry = input("Enter the Chemistry Type: ").strip()
        datasetsDir = input("Enter the directory where the substrate Data resides: ")
        eliminatedPhrases = listInputs(f"Enter the dataframe eliminated phrases for {chemistry} chemistry: ")
        outputDir = input("Enter the output directory: ")
        figDir = input("Enter the figure Directory: ")
        axisMotifs = createAxisMotifs(2)        
        masterDF = plotSubstratesMain(datasetsDir,chemistry , figDir  , axisMotifs, eliminatedPhrases , outputDir)
    elif dataOption== 0:
        masaterDF = pd.read_csv(substrateFile)
    main(initDataSets , masterDF , chemistriesDict)