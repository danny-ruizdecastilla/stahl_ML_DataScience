#Master Feature Extractor
import os 
import sys
import glob
from pathlib import Path
import pandas as pd 
import numpy as np
import plotly 
import plotly.graph_objects as go
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))




def compartmentalization(logDir , outputDir):
    


if __name__ == "__main__":
    logDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
