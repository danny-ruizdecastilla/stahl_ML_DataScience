#The goal of this set of scripts is to place oxygen atoms in the correct positions for the FeCl_oxo system.
#Danny Ruiz de Castilla | 10.07.2026
"""
3 choices: 
axial attachement on the opposite side of the Cl
axial attachement replacing the Cl and then bonding the Cl to the O
axial attachement on the opposite side and then another O on the same side as the Cl
"""
from itertools import combinations, permutations
from pathlib import Path
import sys
import numpy as np
import math
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
from DFTWorkflow.dftFeatureExtractorMAST import getAtomCoordsRobust #logFile , xyzStr , commaSplit:int , locationIdx
oxoFeDouble = 1.64 #Angstroms
oxoFeSingle = 1.91 
oxoClSingle = 1.73
oxoClAngle = 123 #degrees
def saveAtomCoords(atomHash , outputFile):
    with open(outputFile , "w") as f:
        f.write(f"{len(atomHash)}\nSavedFromFeCl_oxo.py\n")
        for idx , atom in atomHash.items():
            f.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]}\n")
def angle(a, b, c):
    ba = [a[i] - b[i] for i in range(len(a))]
    bc = [c[i] - b[i] for i in range(len(c))]

    dot = sum(x * y for x, y in zip(ba, bc))
    mag = math.hypot(*ba) * math.hypot(*bc)

    if mag == 0:
        raise ValueError("Cannot calculate an angle with coincident points.")

    cos_theta = max(-1.0, min(1.0, dot / mag))

    return math.degrees(math.acos(cos_theta))
def main(logDir , outputDir , commaSplit , locationIdx , xyzStr):
    logDir = Path(logDir)
    outputDir = Path(outputDir)
    while True:
        oxoTransferType = input(f"Please enter the oxo transfer type:\n1. Axial attachment on the opposite side of the Cl\n2. Axial attachment replacing the Cl and then bonding the Cl to the O\n3. Axial attachment on the opposite side and then another O on the same side as the Cl\n") #Present each type with a description of the oxo transfer type.
        if oxoTransferType in ["1", "2", "3"]:
            break
        else:
            print("Invalid input. Please enter 1, 2, or 3.")
    if not outputDir.exists():
        outputDir.mkdir(parents=True, exist_ok=True)
    for logFile in logDir.glob("*.log"):
        parent = logFile.stem.split(".")[0] #output str 
        atomHash = getAtomCoordsRobust(logFile , xyzStr , commaSplit , locationIdx)
        feIdx = []
        clIdx = []
        pyroleIdx = []
        for idx , atom in atomHash.items():
            if atom[0] == "Fe":
                feIdx.append(idx)
            elif atom[0] == "Cl":
                clIdx.append(idx)
            elif atom[0] == "N":
                pyroleIdx.append(idx)
        #Get all N N N angles and find the one that is closest to 90 degrees. This will make our x and y for change of basis
        allAngles = []
        combos = list(permutations(pyroleIdx, 3))

        for combo in combos:
            #Calculate the angle between the three atoms
            a = np.array(atomHash[combo[0]][1:4])
            b = np.array(atomHash[combo[1]][1:4])
            c = np.array(atomHash[combo[2]][1:4])
            ang = angle(a , b, c)
            allAngles.append((ang , combo[0] , combo[1] , combo[2]))
        closestAngle = min(allAngles, key=lambda x: abs(x[0] - 90))
        A = np.array(atomHash[closestAngle[1]][1:4], dtype=float)
        B = np.array(atomHash[closestAngle[2]][1:4], dtype=float)
        C = np.array(atomHash[closestAngle[3]][1:4], dtype=float)

        xAxis = A - B
        xAxis /= np.linalg.norm(xAxis)

        yAxis = B - C
        yAxis /= np.linalg.norm(yAxis)

        # Gram-Schmidt orthogonalization
        yAxis -= np.dot(yAxis, xAxis) * xAxis
        yAxis /= np.linalg.norm(yAxis)

        zAxis = np.cross(xAxis, yAxis)
        zAxis /= np.linalg.norm(zAxis)

        basisMatrix = np.column_stack([xAxis, yAxis, zAxis])

        assert np.allclose(basisMatrix.T @ basisMatrix, np.eye(3), atol=1e-10)

        feClVec = np.array(atomHash[clIdx[0]][1:4]) - np.array(atomHash[feIdx[0]][1:4])
        feClVec = feClVec / np.linalg.norm(feClVec)
        if np.dot(feClVec, zAxis) < 0:
            zAxis = -zAxis
        center = np.array(atomHash[feIdx[0]][1:4])
        #Change of basis to the new coordinate system
        for idx , atom in atomHash.items():
            atomCoords = np.array(atom[1:4]) - center
            newX = np.dot(atomCoords, xAxis)
            newY = np.dot(atomCoords, yAxis)
            newZ = np.dot(atomCoords, zAxis)
            atomHash[idx][1:4] = [newX , newY , newZ] 
        if oxoTransferType == "1":
            #Axial attachment on the opposite side of the Cl
            newO = [0 , 0 , -oxoFeDouble]
            atomHash[len(atomHash)] = ["O" , newO[0] , newO[1] , newO[2]]
            saveAtomCoords(atomHash , outputDir / f"{parent}_oxoType1.xyz")
        elif oxoTransferType == "2":
            #Axial attachment replacing the Cl and then bonding the Cl to the O
            newO = [0 , 0 , -oxoFeDouble]
            atomHash[len(atomHash)] = ["O" , newO[0] , newO[1] , newO[2]]
            #Move Cl to be bonded to O
            newCl = [0 , 0 , -oxoFeDouble - oxoClSingle]
            
            atomHash[clIdx[0]][1:4] = newCl
            saveAtomCoords(atomHash , outputDir / f"{parent}_oxoType2.xyz")
        elif oxoTransferType == "3":
            #Axial attachment on the opposite side and then another O on the same side as the Cl
            newO1 = [0 , 0 , -oxoFeDouble]
            atomHash[len(atomHash)] = ["O" , newO1[0] , newO1[1] , newO1[2]]
            #Calculate position for second O based on Cl position and angle
            clPos = np.array(atomHash[clIdx[0]][1:4])
            angleRad = np.radians(oxoClAngle)
            distance = oxoClSingle
            newO2X = clPos[0] + distance * np.sin(angleRad)
            newO2Y = clPos[1]
            newO2Z = clPos[2] + distance * np.cos(angleRad)
            atomHash[len(atomHash)+ 1 ] = ["O" , newO2X , newO2Y , newO2Z]
            saveAtomCoords(atomHash , outputDir / f"{parent}_oxoType3.xyz")
                



                    



