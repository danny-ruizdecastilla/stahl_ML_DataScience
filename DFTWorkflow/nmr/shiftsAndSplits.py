#Daniel Ruiz de Castilla 07.12.2026
#Takes in log files and calculates their isotropic shifts and coupling constants 
import os 
import sys
import pandas as pd 
import numpy as np
from pathlib import Path
import chemdraw
import plotly 
from rdkit import Chem
from networkx import Graph
import re
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
from DFTWorkflow.dftFeatureExtractorMAST import regexLocateInLog, getBoltzmannWeightsGauss , getAtomCoordsRobust , extractShiftsByIdx 
def getProtons(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    protonHash = {}
    atomCount = 0
    for atom in mol.GetAtoms():
        atomCount +=1 
        if atom.GetAtomicNum() == 1:
            continue
        protonHash[atom.GetIdx()] = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 1]

    return protonHash , mol , atomCount
def extractJCoupling(logFile , extract , atoms , downSelectIdx = None):
    lowerInd = regexLocateInLog(logFile , extract[0], 0)
    upperInd = regexLocateInLog(logFile , extract[1], 4 )
    if upperInd == "Poison" or lowerInd == "Poison":
        raise ValueError(f"The Log file {logFile} did not terminate properly")
    JCoupling = np.zeros((atoms, atoms)) #empty zero by zero matrix for j coupling
    with open(logFile , "r") as f:
        for idx, line in enumerate(f):
            if idx > int(lowerInd)  and idx < int(upperInd)-1:
                if "D" in line:
                    line = re.sub(r'([0-9])([Dd])([+-]\d+)', r'\1E\3', line)
                    parts = line.split()
                    row = int(parts[0]) - 1
                    values = [float(x) for x in parts[1:]]
                    colIdx = 0
                    for coupling in values:
                        JCoupling[row][newCols[colIdx]] = coupling
                        colIdx +=1
                else:
                    newCols = [int(shift) - 1 for shift in line.split()]
    JCoupling = JCoupling + JCoupling.T - np.diag(np.diag(JCoupling))
    if downSelectIdx is not None:
        downSelectIdx = np.asarray(downSelectIdx, dtype=int)
        JCoupling = JCoupling[np.ix_(downSelectIdx, downSelectIdx)]
    return JCoupling

def protonNMRData(smiles , comFiles , idStr , logEnergyStr , m , b , mHz ):
    protonHash , molec , totAtoms= getProtons(smiles)
    boltzmannDF = getBoltzmannWeightsGauss(comFiles, 298, "electronic" , logEnergyStr)
    protonIdx = [p + 1 for protons in protonHash.values() for p in protons]
    print(protonIdx)
    protonShiftHash  = {proton : [] for proton in protonIdx}
    jCouplingList = []
    for name in list(boltzmannDF["logID"]):
        fileStr = f"{name}.log"
        conformer  = next((f for f in comFiles if fileStr in f.name), None)
        protonShifts = extractShiftsByIdx(conformer, "SCF GIAO Magnetic shielding tensor (ppm):", "Eigenvalues:  ", "earliest", 
                                          "latest", protonIdx, m , b)
        for proton , shifts in protonShifts.items():
            protonShiftHash[proton].append(shifts[1])
        jCoupling = extractJCoupling(conformer, [ "Total nuclear spin-spin coupling J" , "End of Minotr F.D. properties file"] , 
                                     totAtoms , [p - 1 for p in protonIdx])
        jCouplingList.append(jCoupling)
    weights = boltzmannDF["boltzWeights"]
    avgDelta = [((protonShiftHash[proton] * weights).sum() / weights.sum() * mHz) for proton in protonIdx]
    print(protonShifts.keys())
    jStack = np.stack(jCouplingList)   # shape = (nConf, N, N)
    jAvg = np.tensordot(weights, jStack, axes=(0, 0)) / weights.sum()
    
    return avgDelta , jAvg

        
