#Master Feature Extractor
import os 
import sys
import glob
import tkinter
from pathlib import Path
import pandas as pd 
import numpy as np
import chemdraw
import base64
from io import BytesIO
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from DFTWorkflow.cleanLogs import basicTerm
from dimensionalityReduction.reactivityFeatures import boxGen

def C13Assign(molecHash , subMolec): #for subMolecules that experience symmetry C13 shifts are used to assign C1 and C2 
    #molecHash = {id : {smiles: smilesStr , subMolecId : {nmrHash} } , }
    config = chemdraw.Config()
    config.atom_numbers.show = True
    drawer = chemdraw.Drawer(subMolec, title=subMolec, config= config)
    buf = BytesIO()
    fig = drawer.draw()
    fig.write_image(buf, format="png")  # needs kaleido installed
    buf.seek(0)
    subMolec64 = base64.b64encode(buf.read()).decode("utf-8")
def configureSubstructure(smilesMAST , idMAST , filteredLogs , symmetryCheck):


def compartmentalization(logDir , outputDir , substrateFile):
    substrateDF = pd.read_csv(substrateFile)
    cols = list(substrateDF.columns())
    boxCols = boxGen(cols)
    while True:
        colId = input(f"Here are the columns for the file {substrateFile}\n {boxCols}\n Enter the ID number for the column corresponding to the SMILES strings").strip()
        try:
            smilesId = int(colId)
            smilesCol = cols[smilesId]
            smilesMAST = substrateDF[smilesCol]
            break
        except:
            print("Try again, enter an integer")
    while True:
        colId = input(f"Here are the columns for the file {substrateFile}\n {boxCols}\n Enter the ID number for the column corresponding to the substrate ID strings").strip()
        try:
            substrateId = int(colId)
            subId = cols[substrateId]
            idMAST = substrateDF[subId]
            break
        except:
            print("Try again, enter an integer")
    


    logPaths = Path(logDir) 
    logFiles = logPaths.glob('*.log')
    fileSplit = input(f"{logFiles[0]} Enter the string iteral that seperates the common name with the conf. type : ")

    filteredLogs = []
    for log in logFiles:
        termError = basicTerm(log , "Error termination" , "Normal termination")
        if not termError:
            filteredLogs.append(log)

    subMolecHash , substrateIDHash = configureSubstructure(smilesMAST , idMAST , filteredLogs , True)

if __name__ == "__main__":
    logDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    substrateCSV = str(sys.argv[3])