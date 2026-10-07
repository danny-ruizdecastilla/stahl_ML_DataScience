import os 
import sys
import glob
import pandas as pd 
import numpy as np
import plotly
from pathlib import Path
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from DFTWorkflow.dftFeatureExtractor import getBoltzmannWeightsGauss
def splitGather(logDir , outputDir , outputSplit , energyStr):
    # outputSplit splits all logs into the corresponding molecular conformers 
    # logDir is the path of all log files 
    # output Dir is where we will place all the interactive boltzmann diagrams
    allLogs = list(Path(logDir).glob("*.log"))
    molecs = set(main for main in [log.stem.split(outputSplit)[0] for log in allLogs])
    for mol in molecs:
        molLogs = [log for log in allLogs if log.stem.split(outputSplit)[0] == mol]
        boltzDF = getBoltzmannWeightsGauss(molLogs, 298, "electronic" , energyStr) #cols : 'logID', 'E_Ha', 'rE_kCal', 'boltzWeights'
        weights = boltzDF['boltzWeights'].values
        names = boltzDF['logID'].values
        fig = plotly.graph_objects.Figure(data=[plotly.graph_objects.Bar(x=names, y=weights)])
        fig.update_layout(title=f"Boltzmann Distribution for {mol}",
                          xaxis_title="Conformers",
                          yaxis_title="Boltzmann Weights",
                          yaxis=dict(range=[0, 1]))
        fig.write_html(os.path.join(outputDir, f"{mol}_boltzmann_distribution.html"))
