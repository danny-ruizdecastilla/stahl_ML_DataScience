import os 
import sys
import glob
import pandas as pd 
import numpy as np
import plotly
from pathlib import Path
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
def extractEnergiesSimple(logFile ,energyStr):
    patterns = {
        'electronic':r'SCF Done', 
        'gibbs': r'Sum of electronic and thermal Free Energies',
        'enthalpy': r'Sum of electronic and thermal Enthalpies',
        'zpe': r'Sum of electronic and zero-point Energies'
    }
    try:
        energyType = patterns[energyStr]
    except:
        raise ValueError(f"Unknown energy type: {energyStr}. "
                f"Available types: {list(patterns.keys())}")
    energyLevels = []
    with open(logFile, 'r') as f:
        for i, line in enumerate(f):
            if energyType in line:
                energyLevel = float(line.strip().split("=")[-1].split("A.U.")[0].strip())
                energyLevels.append(energyLevel)
    if len(energyLevels) > 0:
        return energyLevels[-1]
    else:
        return "Poison"
def getBoltzmannWeightsGauss(logDirs, temperature, energyStr , logEnergyStr):
    R = 1.987204e-3
    HARTREE_TO_KCAL = 627.5094740631
    results = []
    
    for log in logDirs:
        if not log.exists():
            print(f"Warning: File {str(log)} not found. Skipping.")
            continue        
        try:
            energy = extractEnergiesSimple(log,  energyStr)
            if energy != "Poison":
                results.append({'logID': Path(log.name).stem,  'E_Ha': energy})
            else:
                print(f"Warning: Could not extract {energyStr} energy from {log}")     
        except Exception as e:
            print(f"Error processing {str(log)}: {str(e)}")
            continue
    if len(results) == 0:
        raise ValueError(f"No energies could be extracted from any log files for the log types {logDirs[-1]}")
    df = pd.DataFrame(results)
    df['E_kCal'] = df['E_Ha'] * HARTREE_TO_KCAL
    minE = df['E_kCal'].min()
    df['rE_kCal'] = df['E_kCal'] - minE
    df['boltzmannFacts'] = np.exp(-df['rE_kCal'] / (R * temperature))
    normFacts = df['boltzmannFacts'].sum()
    df['boltzWeights'] = df['boltzmannFacts'] / normFacts

    df = df.sort_values('E_kCal').reset_index(drop=True)
    
    return df[['logID', 'E_Ha', 'rE_kCal', 'boltzWeights']]
def splitGather(logDir , outputDir , outputSplit , energyStr):
    # outputSplit splits all logs into the corresponding molecular conformers 
    # logDir is the path of all log files 
    # output Dir is where we will place all the interactive boltzmann diagrams
    allLogs = list(Path(logDir).glob("*.log"))
    molecs = set(main for main in [log.stem.split(outputSplit)[0] for log in allLogs])
    for mol in molecs:
        molLogs = [log for log in allLogs if log.stem.split(outputSplit)[0] == mol]
        print(f"Processing {mol} with {molLogs} conformers")
        boltzDF = getBoltzmannWeightsGauss(molLogs, 298, "electronic" , energyStr) #cols : 'logID', 'E_Ha', 'rE_kCal', 'boltzWeights'
        weights = boltzDF['boltzWeights'].values
        names = boltzDF['logID'].values
        fig = plotly.graph_objects.Figure(data=[plotly.graph_objects.Bar(x=names, y=weights)])
        fig.update_layout(title=f"Boltzmann Distribution for {mol}",
                          xaxis_title="Conformers",
                          yaxis_title="Boltzmann Weights",
                          yaxis=dict(range=[0, 1]))
        fig.write_html(os.path.join(outputDir, f"{mol}_boltzmann_distribution.html"))
if __name__ == "__main__":
    logDir = sys.argv[1]  # Directory containing log files
    outputDir = sys.argv[2]  # Directory to save the output HTML files
    #make outputDir if it doesn't exist
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    outputSplit = sys.argv[3]  # String to split the log file names
    energyStr = sys.argv[4]  # Energy type to use for Boltzmann weights (e.g., "electronic")
    
    splitGather(logDir, outputDir, outputSplit, energyStr)