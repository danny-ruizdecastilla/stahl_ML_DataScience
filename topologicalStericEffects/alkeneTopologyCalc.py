#Based on DOI: 10.1021/ci034266b 
#Modified to fit alkene structures by Danny Ruiz de Castilla
import numpy as np 
import argparse
import pandas as pd 
import sys 
from pathlib import Path
from networkx import Graph
from rdkit import Chem
from rdkit.Chem import rdmolops
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from breadthFirstSearch.alkeneBreadthFirstSearchMAST import initAlkene , breadthFirstSearch
radii = {
    "H":  [32,  None, None],
    "B":  [85,  78,  73],
    "C":  [75,  67,  60],
    "N":  [71,  60,  54],
    "O":  [63,  57,  53],
    "F":  [64,  59,  53],
    "Si": [116, 107, 102],
    "P":  [111, 102, 94],
    "S":  [103, 94,  95],
    "Cl": [99,  95,  93],
    "Br": [114, 109, 110],
    "I":  [133, 129, 125],
    "Sn": [140, 130, 132]
}
def calcTopology(atomHash , molec):
    v_rc = 0
    mainAtom = atomHash[0][0]
    atomMain = molec.GetAtomWithIdx(mainAtom)
    
    for period , contacts in atomHash.items():
        for idx in contacts:
            path = rdmolops.GetShortestPath(molec, mainAtom, idx)
            if len(path) > 2:
                for id in path[1:-1]:
                    atom = molec.GetAtomWithIdx(id)
            else:
                

    return v_rc
def main(smilesList , saveDir , period , addHs):
    maxTopology = []
    minTopology = []
    for smiles in smilesList:
        CC , molec = initAlkene(smilesList , addHs)
        c1Hash = breadthFirstSearch(molec, period , {0 : [CC[0]]} , [CC[1]])


        c2Hash = breadthFirstSearch(molec, period , {0 : [CC[1]]} , [CC[0]])


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="alkene Topology Calculator")

    parser.add_argument("--smiles_csv", type=str, required=True,
                        help="Path to CSV file containing SMILES")

    parser.add_argument("--smiles_col", type=str, required=True,
                        help="Column inside .csv to call the smiles list")

    parser.add_argument("--save_dir", type=str, required=True,
                        help="Directory to save outputs")

    parser.add_argument("--period", type=int, required=True,
                        help="BFS expansion depth")

    parser.add_argument("--addHs", action="store_true",
                        help="Include hydrogens if flag is present")
    args = parser.parse_args()
    smilesList = pd.read_csv(args.smiles_col)
    savePath = Path(args.save_dir)
    savePath.mkdir(parents=True, exist_ok=True)
    main(smilesList , savePath ,args.period, args.addHs )
    print(args)

