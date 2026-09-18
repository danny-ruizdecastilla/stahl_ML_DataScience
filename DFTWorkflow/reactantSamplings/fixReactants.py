#Danny Ruiz de Castilla | 09.14.2026
#Filters out conformers based on fixed bond lengths 
import sys
from pathlib import Path
import numpy as np
def newBoxGen(indexList , inputList , menuStr):
    """
    Generates a formatted table-of-contents style menu string 
    from a list of indices and their corresponding inputs.
    """
    box_lines = [
        "┌────────────────────────┐",
        f"│{menuStr}│",
        "├────────────────────────┤"
    ]
    
    # Pair up indices and inputs safely
    for idx, item in zip(indexList, inputList):
        # Format each line to fit nicely inside the box
        box_lines.append(f"│  [{idx}] {str(item):<16} │")
        
    box_lines.append("└────────────────────────┘")
    
    return "\n".join(box_lines)
def read_conformers(confFile):
    conformerHash = {}
    moleculeHash = None
    creatingMolecule = False
    with open(confFile , "r", encoding="utf-8") as f:
        for idx , line in enumerate(f):
            if line.strip().isdigit():
                #New conformer 
                conformer = {}
                if moleculeHash is None:
                    atomTot = int(line.strip())
                    
                    creatingMolecule = True
                    moleculeHash = {}
                else:
                    creatingMolecule = False 
                atomCount = 1
            elif len(line.strip().split(".")) == 2:
                #New Energy 
                energy = float(line.strip())

            else:
                #Just collect coords 
                options = line.strip().split("      ")
                atom = options[0]
                x = float(options[1].strip())
                y = float(options[2].strip())
                z = float(options[3].strip())
                conformer[f"{atom}_{atomCount}"] = [x,y,z]
                if creatingMolecule:
                    #Create atom indices in moleculeHash 
                    moleculeHash[atomCount] = atom 
                if atomCount == atomTot:
                    #Close out the conformer
                    conformerHash[energy] = conformer
                atomCount += 1 
    return conformerHash , moleculeHash 
def main(conformerStr , outputPath  ):
    conformerFile = Path(conformerStr)
    conformerHash , moleculeHash = read_conformers(conformerFile)
    indices = moleculeHash.keys() 
    atoms = moleculeHash.values()
    inputStr = newBoxGen(indices , atoms  , "Select the atom index pairs you want to constrain")
    indexPairs = [] 
    while True:
        indexprs = input(inputStr)
        if indexprs == "break":
            break
        tol = input(f"For pair {indexprs} enter the distance threshold for screening in Angstroms:\n")
        atomChoices = [atom for atom in indexprs.split(",")]
        print("atomChoices" , atomChoices)
        atomChoices.append(float(tol.strip()))
        indexPairs.append(atomChoices)
    conf = 0
    for energy, coords in conformerHash.items():
        #print(coords.keys())
        conf += 1
        valid = True
        for pair in indexPairs:
            atom1 = np.array(coords[pair[0]])
            atom2 = np.array(coords[pair[1]])
            thresh = pair[2]
            dist = atom1 - atom2 
            
            if np.linalg.norm(dist) > thresh:
                valid = False
                break  # Break out of the inner loop

            #save coords 
            totAtoms = len(coords.keys())
            with open(outputPath / f"{conf}.xyz", "w", encoding="utf-8") as f:
                f.write(f"{totAtoms}\n")
                f.write(f"Created in DFTWorkflow\n")
                for atom, coordinates in coords.items():
                    atomStr = atom.split("_")[0]
                    f.write(f"{atomStr:<4} {coordinates[0]} {coordinates[1]} {coordinates[2]}\n")
        if not valid:
            continue
if __name__ == "__main__":
    confPath = Path(sys.argv[1])
    outputPath = Path(sys.argv[2])
    outputPath.mkdir(parents=True, exist_ok=True)
    main(confPath , outputPath  )



        