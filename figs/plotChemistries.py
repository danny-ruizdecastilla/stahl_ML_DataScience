import pandas as pd
import os 
import sys
import glob
import json
import numpy as np
import html
import plotly
import shutil
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from figs.chemPlotlyV2 import createPNGDF,png64
from figs.chemPlotlyV1 import convertCanonical
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
def htmlGenerator2(jsonDict , axisList , chemStr , outputDir):

    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Scatter Plot with Hover Images</title>
        <script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script>
        <script>
            const jsonDict = {{json.dumps(jsonDict)}};
            const groupedData = jsonDict;

            const dropdown = documnet,getElementById("groupDropDown"):
            for (const groupName in groupedData){{
            const option = documnet.createElement("option");
            option.value = groupName;
            option.textContent = groupName;
            dropdown.appendChild(option);
            
            }}
    """
def main(chemistryFiles, masterDF, chemistriesDict , chemistry, outputDir):
    cols = masterDF.columns
    if "base64" in cols:
        pass
    else:
        figDir = input("Enter the figure Directory where PNG's of the substrates will be stored (temporarily): ")
        masterDF = createPNGDF(masterDF ,"SMILES" , figDir)
        base64Col = []
        for img in list(masterDF["pngPath"]):

            base64 = png64(img)
            base64Col.append(base64)
            masterDF["base64"] = base64Col
            masterDF = masterDF.drop(columns=['pngPath'])
        shutil.rmtree(figDir)
    jsonDict = {}
    for file in chemistryFiles:
        chemName = file.split("/")[-1].split(".")[0]
        df = pd.DataFrame(columns=list(masterDF.columns) + ["Yield"])
        chemDF = pd.read_csv(file)  
        for _, row in chemDF.iterrows():
            smiles = row["SMILES"]
            canonical = convertCanonical(smiles)
            matches = masterDF[masterDF["canonicalSMILES"] == canonical]
            if not matches.empty:
                rowMatch = matches.head(1)
                rowMatch["Yield"] = float(row["Yields"])
                df = pd.concat([df, rowMatch], ignore_index=True)
            else:
                continue
        jsonChem = chemDF.to_json(orient="records" , force_ascii=False )
        jsonChem = jsonChem.replace('\\/', '/')
        dfName = str(chemistriesDict[chemName])
        jsonDict[dfName] = jsonChem
    axisList = masterDF.select_dtypes(include='number').columns
    htmlGenerator2(jsonDict , axisList , chemistry )


if __name__ == "__main__":
    chemistriesDir = str(sys.argv[1])
    substrateFile = str(sys.argv[2])
    chemistry = str(sys.argv[3])
    outputDir = str(sys.argv[4])
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
        datasetsDir = input("Enter the directory where the substrate Data resides: ")
        eliminatedPhrases = listInputs(f"Enter the dataframe eliminated phrases for {chemistry} chemistry: ")
        figDir = input("Enter the figure Directory: ")
        axisMotifs = createAxisMotifs(2)        
        masterDF = plotSubstratesMain(datasetsDir,chemistry , figDir  , axisMotifs, eliminatedPhrases , outputDir)
    elif dataOption== 0:
        masaterDF = pd.read_csv(substrateFile , encoding='utf-8')
    main(initDataSets , masterDF , chemistriesDict, chemistry, outputDir)