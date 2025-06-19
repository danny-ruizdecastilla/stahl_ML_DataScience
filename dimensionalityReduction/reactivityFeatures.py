import os 
import sys
import glob
import shutil
import re
import pandas as pd 
import numpy as np
import plotly 
import plotly.graph_objects as go
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
def selectChemistries(mainDir ):
    csvPool = glob.glob(mainDir + "/*.csv")
    csvOptions = [csv.split("/")[-1] for csv in csvPool]
    prompt1 = '''Here are the dataframe options from this directory:
            {csvOptions[i]} === [{i}]
            '''
    while True:
        
if __name__ == "__main__":
    dataDir = str(sys.argv[1])
    chemList = selectChemistries(dataDir)

