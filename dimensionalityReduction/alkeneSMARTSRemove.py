#Remove functional groups based on robustness screen results
from rdkit import Chem
import chemdraw
import pandas as pd
import sys
from pathlib import Path
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from figs.chemPlotlyV1 import dat2List
def saveToChemdraw(accepted:list , outputDir , conditionStr , status):
    if len(accepted) <= 100:
        acceptedChems = chemdraw.GridDrawer(accepted)
        acceptedChems.draw_png(Path(outputDir) / f"{conditionStr}_{status}.png" )
    else:
        pop = 100
        acceptedBins = [accepted[i:i+pop] for i in range(0, len(accepted), pop)]
        for i in range (len(acceptedBins)):
            bin_ = acceptedBins[i]
            acceptedChems = chemdraw.GridDrawer(bin_)
            condition = conditionStr + "_" +  str(i)
            acceptedChems.draw_png(Path(outputDir) / f"{condition}_{status}.png" )
def main(dfMAST , outputDir , smartsSubstructures , saveStr ):
    structureList = dfMAST["SMILES"]
    
    acceptedHash = {}
    acceptedHash["initAlkenes"] = structureList
    for smarts in smartsSubstructures:
        acceptedList = []
        rejectedList = []
        subStruct = Chem.MolFromSmarts(smarts)
        smartsString = input(f"Input the save string for this smarts substructure {smarts}: ")
        structureStr = next(reversed(acceptedHash))
        structures = acceptedHash[structureStr]
        for smiles in structures:
            alkene = Chem.MolFromSmiles(smiles)
            indeces = alkene.GetSubstructMatches(subStruct)
            if len(indeces) == 0:
                #No substructure exists 
                acceptedList.append(smiles)
            else:
                rejectedList.append(smiles)
        acceptedHash[smartsString] = acceptedList
        saveToChemdraw(acceptedList , outputDir , smartsString , "accepted")
        saveToChemdraw(rejectedList , outputDir , smartsString , "rejected")

    finalstr = next(reversed(acceptedHash))
    finalStructures = acceptedHash[finalstr]
    dfFinal = dfMAST[dfMAST["SMILES"].isin(finalStructures)].copy()
    dfFinal.to_csv(Path(outputDir) / f"{saveStr}_removedFunctionals.csv")
    
if __name__ == "__main__":
    chemistryCSV = str(sys.argv[1])
    dfMAST = pd.read_csv(chemistryCSV)
    chemdrawDir = str(sys.argv[2])
    pathdir = Path(chemdrawDir)
    pathdir.mkdir(parents=True, exist_ok=True)
    smartsDat = str(sys.argv[3])
    smartsList = dat2List(smartsDat , ",")
    dfpath = Path(chemistryCSV)
    name = dfpath.name.split(".")[0]
    main(dfMAST,chemdrawDir, smartsList , name)
