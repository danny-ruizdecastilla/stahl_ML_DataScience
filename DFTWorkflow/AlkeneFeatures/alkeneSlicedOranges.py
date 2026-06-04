#Vburr for alkenes with quadrants and octants 
#Danny Ruiz de Castilla 05.25.2026
#Based on doi.org/10.1021/acsatal.2c00976
import numpy as np
import sys
from pathlib import Path
from networkx import Graph
from rdkit import Chem
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
from breadthFirstSearch.radialBasedCorrelation import getCC
from figs.stericVisuals import stericDrawer
class alkeneSlicedOranges:
    def __init__(self, C1Hash , C2Hash, radius):
        self.C1Idx = int(C1Hash["idx"]) #Index of alkene Carbon 1
        self.C2Idx = int(C2Hash["idx"]) #Index of alkene Carbon 2

        self.C1Coords = C1Hash["0"]
        self.C2Coords = C2Hash["0"]

        #XAxis
        XVec = self.C1Coords - self.C2Coords
        XVec /= np.linalg.norm(XVec)

        # vectors from carbons to neighbors
        c1Contacts = C1Hash["1"]
        c2Contacts = C2Hash["1"]

        c1Vec1 = c1Contacts[0] - self.C1Coords
        c1Vec2 = c1Contacts[1] - self.C1Coords

        c2Vec1 = c2Contacts[0] - self.C2Coords
        c2Vec2 = c2Contacts[1] - self.C2Coords

        # plane normals
        c1Orth = np.cross(c1Vec1, c1Vec2)
        c2Orth = np.cross(c2Vec1, c2Vec2)

        # normalize
        c1Orth /= np.linalg.norm(c1Orth)
        c2Orth /= np.linalg.norm(c2Orth)

        # align directions
        if np.dot(c1Orth, c2Orth) < 0:
            c2Orth *= -1

        # scale
        c1OrthScaled = radius * c1Orth
        c2OrthScaled = radius * c2Orth

        #YAxis
        yTargetC1 = self.C1Coords + c1OrthScaled 
        yTargetC2 = self.C2Coords + c2OrthScaled

        yVecC1 = (yTargetC1 - self.C1Coords) / np.linalg.norm(yTargetC1 - self.C1Coords)
        yVecC2 = (yTargetC2 - self.C2Coords) / np.linalg.norm(yTargetC2 - self.C2Coords)
        #Gram Schmidt for y axis 
        yVecC1 = yVecC1 - np.dot(yVecC1, XVec) * XVec
        yVecC1 /= np.linalg.norm(yVecC1)
        yVecC2 = yVecC2 - np.dot(yVecC2, XVec) * XVec
        yVecC2 /= np.linalg.norm(yVecC2)

        self.radius = radius

        #z axis 
        zVecC1 = np.cross(XVec, yVecC1)
        zVecC1 /= np.linalg.norm(zVecC1)
        zVecC2 = np.cross(XVec, yVecC2)
        zVecC2 /= np.linalg.norm(zVecC2)

        #change of basis 
        self.basisMatrixC1 = np.column_stack([XVec, yVecC1, zVecC1])
        self.basisMatrixC2 = np.column_stack([XVec, yVecC2, zVecC2])
        #print(self.basisMatrix)
        self.basisMatrixInvC1 = self.basisMatrixC1.T
        self.basisMatrixInvC2 = self.basisMatrixC2.T
        #print(self.basisMatrixInv@self.basisMatrix)

        assert np.allclose(self.basisMatrixInvC1 @ self.basisMatrixC1, np.eye(3), rtol=1e-5)  # should be identity
        assert np.allclose(self.basisMatrixInvC2 @ self.basisMatrixC2, np.eye(3), rtol=1e-5)
    def getAtoms(self, atomList , idxList , symbolList , bondHash, makeFig ):
        #x1Atoms 
        atomsC1 = np.array(atomList)
        translatedC1 = atomsC1 - self.C1Coords
        newAtomsC1 = (self.basisMatrixInvC1 @ translatedC1.T).T

        acceptedAtomsC1 = []
        acceptedIdxC1 = []
        acceptedSymbolsC1 = []
        for i in range(len(newAtomsC1)):
            atom = newAtomsC1[i]
            rad = np.linalg.norm(atom)
            if rad <= self.radius:
                acceptedAtomsC1.append(atom)
                acceptedIdxC1.append(idxList[i])
                acceptedSymbolsC1.append(symbolList[i])
        self.atomCoordsC1 = np.array(acceptedAtomsC1)
        self.atomIndexC1 = acceptedIdxC1 
        self.atomSymbolsC1 = acceptedSymbolsC1

        #x2Atoms 
        atomsC2 = np.array(atomList)
        translatedC2 = atomsC2 - self.C2Coords
        newAtomsC2 = (self.basisMatrixInvC2 @ translatedC2.T).T
        acceptedAtomsC2 = []
        acceptedIdxC2 = []
        acceptedSymbolsC2 = []
        for i in range(len(newAtomsC2)):
            atom = newAtomsC2[i]
            rad = np.linalg.norm(atom)
            if rad <= self.radius:
                acceptedAtomsC2.append(atom)
                acceptedIdxC2.append(idxList[i])
                acceptedSymbolsC2.append(symbolList[i])
        self.atomCoordsC2 = np.array(acceptedAtomsC2)
        self.atomIndexC2 = acceptedIdxC2 
        self.atomSymbolsC2 = acceptedSymbolsC2
        print(acceptedIdxC1)
        print(acceptedIdxC2)
        acceptedIdx = list(set(acceptedIdxC2 + acceptedIdxC1))
        print(acceptedIdx)
        if makeFig:
            pt = Chem.GetPeriodicTable()
            radiiList = []
            for symbol in symbolList:
                atomRadii = pt.GetRvdw(pt.GetAtomicNumber(symbol))
                radiiList.append(atomRadii)
            figHash = {"atomCoords" : atomList , "atomSymbol" : symbolList , "radii" : radiiList}
            stericFigure = stericDrawer(figHash , 30 , highlighted = acceptedIdx)
            for i in range(len(atomList)):
                stericFigure.drawAtom(i)
            for _ , bdTable in bondHash.items():
                bdAtoms =  bdTable["idxList"]
                bondType = str(bdTable["bondType"])
                stericFigure.drawBond( bdAtoms[0], bdAtoms[1], bondType)
            self.Fig = stericFigure
    def slicedOranges(self, quad , oct , makeFig ):
        sphereVol = (4/3)*np.pi*self.radius**3
        if not hasattr(self, "atomCoordsC1"):
            raise RuntimeError("getAtoms() must be called before getBurriedVolume()")
        from rdkit import Chem
        def getAtomsVol(symbolList):
            atomVol = 0
            pt = Chem.GetPeriodicTable()
            for atom in symbolList:
                atomRadii = pt.GetRvdw(pt.GetAtomicNumber(atom))
                vol = (4/3) * np.pi * atomRadii**3
                atomVol += vol
            return atomVol
        returnHash = {}
        c1Tot = getAtomsVol(self.atomSymbolsC1)
        c2Tot = getAtomsVol(self.atomSymbolsC2)
        c1Bur = (c1Tot / sphereVol) * 100
        c2Bur = (c2Tot / sphereVol ) * 100
        if c1Bur >= c2Bur:
            flip = True
        else:
            flip = False
        def sliceOrange(coordinates , symbols, hashMap , label):
            quadLabelsC1 = (coordinates[:, 0] > 0).astype(int) + 2 * (coordinates[:, 2] > 0).astype(int)
            c1Idx = [np.where(quadLabelsC1 == i)[0] for i in range(4)]
            quadNum = 0
            for q in c1Idx:
                atoms = [symbols[i] for i in q]
                atomVol = getAtomsVol(atoms)
                hashMap[f"C{label}_quad_{quadNum}"] = float(atomVol / (sphereVol/4)) 
                quadNum += 1
            return hashMap
        if quad: #split along z and x axis 
            #C1 is min burriedVol
            if flip:
                returnHash = sliceOrange(self.atomCoordsC1 , self.atomSymbolsC1, returnHash , "2")
                returnHash = sliceOrange(self.atomCoordsC2 , self.atomSymbolsC2, returnHash , "1")
                if makeFig:
                    c1Hash = {"atom" : self.C2Coords ,  "inverseCoB" : self.basisMatrixInvC2 , "radius" : self.radius , "shapeType" : "slicedOrange"}
                    c2Hash = {"atom" : self.C1Coords ,  "inverseCoB" : self.basisMatrixInvC1 ,"radius" : self.radius , "shapeType" : "slicedOrange"}
                    slicedFig = self.Fig.drawShapes(c1Hash)  
                    #slicedFig = self.Fig.drawShapes(c2Hash, fig = slicedFig)
                    slicedFig.write_html("slicedOranges.html")
            else:
                returnHash = sliceOrange(self.atomCoordsC1 , self.atomSymbolsC1, returnHash , "1")
                returnHash = sliceOrange(self.atomCoordsC2 , self.atomSymbolsC2, returnHash , "2")
                if makeFig:
                    c1Hash = {"atom" : self.C1Coords ,  "inverseCoB" : self.basisMatrixInvC1 , "radius" : self.radius , "shapeType" : "slicedOrange"}
                    c2Hash = {"atom" : self.C2Coords ,  "inverseCoB" : self.basisMatrixInvC2 ,"radius" : self.radius , "shapeType" : "slicedOrange"}
                    slicedFig = self.Fig.drawShapes(c1Hash)  
                    #slicedFig = self.Fig.drawShapes(c2Hash, fig = slicedFig)
                    slicedFig.write_html("slicedOranges.html")
        if oct:
            #C1
            octLabelsC1 = (
                (self.atomCoordsC1[:, 0] > 0).astype(int)
                + 2 * (self.atomCoordsC1[:, 1] > 0).astype(int)
                + 4 * (self.atomCoordsC1[:, 2] > 0).astype(int)
            )
            c1Idx = [np.where(octLabelsC1 == i)[0] for i in range(4)]
            octNum = 0
            for o in c1Idx:
                atoms = [self.atomSymbolsC1[i] for i in o]
                atomVol = getAtomsVol(atoms)
                returnHash[f"C1_oct_{octNum}"] = float(atomVol / (sphereVol/8))
                octNum += 1
            #C2
            octLabelsC2 = (
                (self.atomCoordsC2[:, 0] > 0).astype(int)
                + 2 * (self.atomCoordsC2[:, 1] > 0).astype(int)
                + 4 * (self.atomCoordsC2[:, 2] > 0).astype(int)
            )
            c2Idx = [np.where(octLabelsC2 == i)[0] for i in range(4)]
            octNum = 0
            for o in c2Idx:
                atoms = [self.atomSymbolsC2[i] for i in o]
                atomVol = getAtomsVol(atoms)
                returnHash[f"C2_oct_{octNum}"] = float(atomVol / (sphereVol/8))
                octNum += 1
        return returnHash
def main(logFile , smilesStr, radius , linkIdx):
    from DFTWorkflow.dftFeatureExtractorMAST import getAtomCoordsRobust
    cc , molec = getCC(smilesStr)
    #print(cc)
    C1 = cc[1] + 1
    C2 = cc[0] + 1 
    g = Graph()
    bondHash = {}
    idx = 0
    for bond in molec.GetBonds():

        start, end = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        bondStr = bond.GetBondType()
        #print(bondStr)
        bondHash[idx] = {"idxList" : [start , end ] , "bondType" : str(bondStr)}
        idx +=1 
        g.add_edge(start, end)
    #print(bondHash)
    CminNeighbors = list(g.neighbors(int(C1-1)))
    CminNeighbors.remove(int(C2-1))
    CmaxNeighbors = list(g.neighbors(int(C2-1)))
    CmaxNeighbors.remove(int(C1-1))
    CmaxContacts = []
    CminContacts = []
    coordHash = getAtomCoordsRobust(logFile , "GINC-COMPUTE" , linkIdx  , 1)
    #print(coordHash.keys())
    atoms = []
    idxList = []
    symbols = []
    for atomIdx , coords in coordHash.items():
        idxList.append(atomIdx)
        atoms.append(coords[2:5])
        symbols.append(coords[0])
        atomIdx +=1
        if atomIdx == C1:
            c1_coords = np.array(coords[2:5])
            c1Idx = atomIdx
        elif atomIdx == C2:
            c2_coords = coords[2:5]
            c2Idx = atomIdx
        if (atomIdx-1) in CminNeighbors:
            crds = np.array(coords[2:5])
            CminContacts.append(crds) 
        if (atomIdx-1) in CmaxNeighbors:
            crds = np.array(coords[2:5])
            CmaxContacts.append(crds) 
    C1Hash = {"0" : c1_coords , "1" : CminContacts , "idx" : c1Idx}
    C2Hash = {"0" : c2_coords , "1" : CmaxContacts , "idx" : c2Idx}
    alkeneOranges = alkeneSlicedOranges(C1Hash , C2Hash , radius)
    alkeneOranges.getAtoms(atoms , idxList , symbols , bondHash, True)
    orangeSlices = alkeneOranges.slicedOranges(True , False , True)
    print(orangeSlices)
if __name__ == "__main__":


    logFile = str(sys.argv[1])
    smilesStr = str(sys.argv[2])
    radius = float(sys.argv[3])
    linkIdx = int(sys.argv[4])
    main(logFile , smilesStr , radius , linkIdx)