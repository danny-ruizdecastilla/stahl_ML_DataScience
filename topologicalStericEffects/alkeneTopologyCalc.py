#Based on DOI: 10.1021/ci034266b 
#Modified to fit alkene structures by Danny Ruiz de Castilla
import numpy as np 
import argparse
import pandas as pd 
import sys 
from pathlib import Path
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
    "Sn": [140, 130, 132],
    "Se" : [116 , 107 , 107]
}
def calcTopology(atomHash , molec):
    print(atomHash)
    TSEI = 0
    for period , nestedPaths in atomHash.items():
        for path in nestedPaths:
            print(path)
            if len(path) > 1:
                lastIdx = path[-1]
                denominator = 0
                for i in range(len(path)-1):
                    startIdx = path[i]
                    endIdx = path[i+1]
                    bond = molec.GetBondBetweenAtoms(startIdx, endIdx)
                    bondType = str(bond.GetBondType())
                    print(bondType)
                    if bondType == "SINGLE":
                        bondIdx = 0
                    elif bondType == "AROMATIC" or bondType == "DOUBLE":
                        bondIdx = 1
                    else:
                        bondIdx = 2
                    startAtom = molec.GetAtomWithIdx(startIdx)
                    endAtom = molec.GetAtomWithIdx(endIdx)
                    startSymbol = startAtom.GetSymbol() 
                    endSymbol = endAtom.GetSymbol()
                    startRadii = radii[startSymbol][bondIdx]
                    endRadii = radii[endSymbol][bondIdx]
                    denominator += startRadii
                    denominator += endRadii
                    if endIdx == lastIdx:
                        numerator = endRadii
                v_rc = numerator**3 / denominator**3 
                TSEI += v_rc

    return TSEI
def main(smilesList , saveDir , period , addHs):
    maxTopology = []
    minTopology = []
    for smiles in smilesList:
        CC , molec = initAlkene(smiles , addHs)
        c1Hash = breadthFirstSearch(molec, period , {0 : [[CC[0]]]} , [CC[1]])
        c1TSEL = calcTopology(c1Hash, molec)
        c2Hash = breadthFirstSearch(molec, period , {0 : [[CC[1]]]} , [CC[0]])
        c2TSEL = calcTopology(c2Hash , molec)
        tsel = [c1TSEL , c2TSEL]
        maxTopology.append(max(tsel))
        minTopology.append(min(tsel))
    deltaTopology = np.abs(np.array(maxTopology) - np.array(minTopology))
    meanTopology = (np.array(maxTopology) + np.array(minTopology))/2
    topologyDF = pd.DataFrame()
    topologyDF["SMILES"] = smilesList
    topologyDF[f"maxTopology_{period}"] = maxTopology
    topologyDF[f"minTopology_{period}"] = minTopology
    topologyDF[f"deltaTopology_{period}"] = deltaTopology
    topologyDF[f"meanTopology_{period}"] = meanTopology
    topologyDF.to_csv(saveDir / f"topology_atPeriod_{period}.csv" )
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
    dfMAST = pd.read_csv(args.smiles_csv)
    smilesList = dfMAST[args.smiles_col]
    savePath = Path(args.save_dir)
    savePath.mkdir(parents=True, exist_ok=True)
    print(args)
    main(smilesList , savePath ,args.period, args.addHs )


