#The steric doughnut | Danny Ruiz de Castilla
#Not really a doughnut
import plotly 
import numpy as np
import plotly.graph_objects as go
import sys
from pathlib import Path
from networkx import Graph
from rdkit import Chem
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
from breadthFirstSearch.radialBasedCorrelation import getCC
from figs.stericVisuals import stericDrawer
def vdw_radius(symbol):
    pt = Chem.GetPeriodicTable()
    return pt.GetRvdw(pt.GetAtomicNumber(symbol))
def create_hollow_cylinder(inner_radius, outer_radius, height, resolution=50):
    if inner_radius >= outer_radius:
        raise ValueError("Inner radius must be smaller than outer radius")
    
    # Create theta values for circular cross-sections
    theta = np.linspace(0, 2 * np.pi, resolution)
    
    # Create z values for height
    z = np.linspace(0, height, resolution)
    
    # Create meshgrids
    theta_grid, z_grid = np.meshgrid(theta, z)
    
    # Outer wall
    x_outer = outer_radius * np.cos(theta_grid)
    y_outer = outer_radius * np.sin(theta_grid)
    z_outer = z_grid
    
    # Inner wall
    x_inner = inner_radius * np.cos(theta_grid)
    y_inner = inner_radius * np.sin(theta_grid)
    z_inner = z_grid
    
    # Top ring (z = height)
    r_ring = np.linspace(inner_radius, outer_radius, resolution)
    theta_ring, r_grid = np.meshgrid(theta, r_ring)
    
    x_top = r_grid * np.cos(theta_ring)
    y_top = r_grid * np.sin(theta_ring)
    z_top = np.full_like(x_top, height)
    
    # Bottom ring (z = 0)
    x_bottom = r_grid * np.cos(theta_ring)
    y_bottom = r_grid * np.sin(theta_ring)
    z_bottom = np.zeros_like(x_bottom)
    
    return {
        'outer_wall': (x_outer, y_outer, z_outer),
        'inner_wall': (x_inner, y_inner, z_inner),
        'top_ring': (x_top, y_top, z_top),
        'bottom_ring': (x_bottom, y_bottom, z_bottom)
    }

def check_atoms_in_cylinder(atomHash, inner_radius, outer_radius, height, center):
    #coordHash structure : {key : [atomName, atomInd, x,y,z]}
    atoms = []
    idxList = []
    for _ , val in atomHash.items():
        coords = [float(a) for a in val[2:5]]
        idx = val[1]
        idxList.append(idx)
        atoms.append(coords)
    atoms = np.array(atoms)
    
    if atoms.ndim != 2 or atoms.shape[1] != 3:
        raise ValueError("Atoms array must be of shape (N, 3)")
    
    cx, cy, cz = center
    
    xList = atoms[:, 0] - cx
    yList = atoms[:, 1] - cy
    zList = atoms[:, 2] - cz
    
    # Calculate radial distance from cylinder axis (z-axis)
    r = np.sqrt(xList**2 + yList**2)

    inside_mask = (r >= inner_radius) & (r <= outer_radius) & (zList >= 0) & (zList <= height)
    #print("mask" , inside_mask)
    insideID = []
    inside_indices = np.where(inside_mask)[0]
    
    for i in inside_indices:
        atomID = idxList[i]
        insideID.append(atomID)
    
    return insideID

def get_atoms_summary(inside_atoms_dict, atoms):
    """
    Print summary of atoms inside and outside the cylinder.
    """
    total = len(atoms)
    inside = len(inside_atoms_dict)
    outside = total - inside
    
    #print(f"Total atoms: {total}")
    #print(f"Atoms inside cylinder: {inside}")
    #print(f"Atoms outside cylinder: {outside}")
    #print(f"Percentage inside: {100*inside/total:.1f}%")
    
    return inside, outside

def plot_cylinder_with_atoms(atoms, inside_atoms_dict, inner_radius, outer_radius, 
                             height, center, resolution=50):
    surfaces = create_hollow_cylinder(inner_radius, outer_radius, height, resolution)
    
    # Translate surfaces to cylinder center
    cx, cy, cz = center
    
    fig = go.Figure()
    
    # Plot outer wall (blue)
    x_outer, y_outer, z_outer = surfaces['outer_wall']
    fig.add_trace(go.Surface(
        x=x_outer + cx, y=y_outer + cy, z=z_outer + cz,
        opacity=0.3,
        colorscale=[[0, 'steelblue'], [1, 'steelblue']],
        showscale=False,
        name='Outer Wall',
        hoverinfo='skip'
    ))
    
    # Plot inner wall (lightblue)
    x_inner, y_inner, z_inner = surfaces['inner_wall']
    fig.add_trace(go.Surface(
        x=x_inner + cx, y=y_inner + cy, z=z_inner + cz,
        opacity=0.3,
        colorscale=[[0, 'lightblue'], [1, 'lightblue']],
        showscale=False,
        name='Inner Wall',
        hoverinfo='skip'
    ))
    
    # Plot top ring (coral)
    x_top, y_top, z_top = surfaces['top_ring']
    fig.add_trace(go.Surface(
        x=x_top + cx, y=y_top + cy, z=z_top + cz,
        opacity=0.4,
        colorscale=[[0, 'coral'], [1, 'coral']],
        showscale=False,
        name='Top',
        hoverinfo='skip'
    ))
    
    # Plot bottom ring (coral)
    x_bottom, y_bottom, z_bottom = surfaces['bottom_ring']
    fig.add_trace(go.Surface(
        x=x_bottom + cx, y=y_bottom + cy, z=z_bottom + cz,
        opacity=0.4,
        colorscale=[[0, 'coral'], [1, 'coral']],
        showscale=False,
        name='Bottom',
        hoverinfo='skip'
    ))
    
    # Separate atoms into inside and outside
    atoms = np.array(atoms)
    inside_indices = set(inside_atoms_dict.keys())
    outside_mask = np.array([i not in inside_indices for i in range(len(atoms))])
    
    # Plot atoms outside (red)
    if np.any(outside_mask):
        outside_atoms = atoms[outside_mask]
        fig.add_trace(go.Scatter3d(
            x=outside_atoms[:, 0],
            y=outside_atoms[:, 1],
            z=outside_atoms[:, 2],
            mode='markers',
            marker=dict(size=5, color='red', opacity=0.7),
            name='Outside',
            text=[f'Atom {i}' for i in range(len(atoms)) if i not in inside_indices],
            hovertemplate='<b>%{text}</b><br>x: %{x:.2f}<br>y: %{y:.2f}<br>z: %{z:.2f}<extra></extra>'
        ))
    
    # Plot atoms inside (green)
    if inside_atoms_dict:
        inside_coords = np.array(list(inside_atoms_dict.values()))        # scale
        inside_idx_list = list(inside_atoms_dict.keys())
        fig.add_trace(go.Scatter3d(
            x=inside_coords[:, 0],
            y=inside_coords[:, 1],
            z=inside_coords[:, 2],
            mode='markers',
            marker=dict(size=5, color='green', opacity=0.9),
            name='Inside',
            text=[f'Atom {i}' for i in inside_idx_list],
            hovertemplate='<b>%{text}</b><br>x: %{x:.2f}<br>y: %{y:.2f}<br>z: %{z:.2f}<extra></extra>'
        ))
    
    # Update layout
    max_range = max(outer_radius, height/2) * 1.2
    fig.update_layout(
        title=f'Hollow Cylinder with Atoms<br>Inner R={inner_radius}, Outer R={outer_radius}, Height={height}',
        scene=dict(
            xaxis=dict(range=[-max_range, max_range], title='X'),
            yaxis=dict(range=[-max_range, max_range], title='Y'),
            zaxis=dict(range=[min(0, cz), max(height + cz, cz + height)], title='Z'),
            aspectmode='cube'
        ),
        width=900,
        height=700,
        showlegend=True
    )
    
    fig.show()

def getBurriedDougnut(atomHash, innerR, outerR, h, center ):
    cylinderHash = create_hollow_cylinder(innerR, outerR, h, resolution=50)
    insideList = check_atoms_in_cylinder(atomHash, innerR, outerR, h, center)
    #print("insiderList" , insideList)
    volDoughnut = np.pi*outerR**2*h - np.pi*innerR**2*h
    #print("volDoughnut" , volDoughnut)
    atomVol = 0
    for _ , atoms in atomHash.items():
        atomID = str(atoms[0])
        if atomID in insideList:
            atomRadii = vdw_radius(atomID)
            vol = (4/3) * np.pi * atomRadii**3
            atomVol += vol
    percentBurriedVol = (atomVol/volDoughnut) * 100
    return percentBurriedVol

class alkeneSemiCylinders:
    def __init__(self, C1Hash, C2Hash,radius ,padding ):
        self.padding = padding
        self.C1Idx = int(C1Hash["idx"]) #Index of alkene Carbon 1
        self.C2Idx = int(C2Hash["idx"]) #Index of alkene Carbon 2
        self.C1 =  C1Hash["0"]
        self.C2 =  C2Hash["0"]

        #create new x vector 
        center = [(float(a) + float(b)) * 0.5 for a, b in zip(list(C1Hash["0"]), list(C2Hash["0"]))]
        self.centerPoint = center 
        xVec = self.C1 - self.centerPoint #unit x axis starts at center between C1 and C2 and goes towards C1
        xVec /= np.linalg.norm(xVec)
        self.xVec = xVec
        alkeneVec = self.C1- self.C2 #vector between C1 and C2
        self.alkeneVec = alkeneVec 
        self.height = np.linalg.norm(alkeneVec)#alkeneBondLength 
        # vectors from carbons to neighbors
        c1Contacts = C1Hash["1"]
        c2Contacts = C2Hash["1"]
        c1Vec1 = c1Contacts[0] - C1Hash["0"]
        self.c1Vec1 = c1Vec1
        c1Vec2 = c1Contacts[1] - C1Hash["0"]
        self.c1Vec2 = c1Vec2
        c2Vec1 = c2Contacts[0] - C2Hash["0"]
        self.c2Vec1 = c2Vec1
        c2Vec2 = c2Contacts[1] - C2Hash["0"]
        self.c2Vec2 = c2Vec2

        # cap heights 
        alkeneUnit = alkeneVec / self.height
        c1Proj1 = np.dot(c1Vec1, alkeneUnit)
        c1Proj2 = np.dot(c1Vec2, alkeneUnit)
        c1CapHeight = max([np.linalg.norm(c1Proj1) , np.linalg.norm(c1Proj2)])
        self.C1CapHeight = c1CapHeight + 0.45
        #print(self.C1CapHeight)

        c2Proj1 = np.dot(c2Vec1, alkeneUnit)
        c2Proj2 = np.dot(c2Vec2, alkeneUnit)
        c2CapHeight = max([np.linalg.norm(c2Proj1) , np.linalg.norm(c2Proj2)])
        self.C2CapHeight = c2CapHeight + 0.45
        #print(self.C2CapHeight)

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
        bisector = (c1OrthScaled + c2OrthScaled ) #bisector between the 2 orthonormal vectors of the alkene planes 
        bisector /= np.linalg.norm(bisector)
        yTarget = self.centerPoint + bisector
        yVec = (yTarget - self.centerPoint) / np.linalg.norm(yTarget - self.centerPoint)
        bisector *= radius
        self.bisector = bisector
        self.radius = radius

        #Gram Schmidt for y axis 
        yVec = yVec - np.dot(yVec, xVec) * xVec
        yVec /= np.linalg.norm(yVec)
        self.yVec = yVec
        #z axis 
        zVec = np.cross(self.xVec, self.yVec)
        zVec /= np.linalg.norm(zVec)
        self.zVec = zVec

        #change of basis 
        self.basisMatrix = np.column_stack([self.xVec, self.yVec, self.zVec])
        #print(self.basisMatrix)
        self.basisMatrixInv = self.basisMatrix.T
        #print(self.basisMatrixInv@self.basisMatrix)

        assert np.allclose(self.basisMatrixInv @ self.basisMatrix, np.eye(3), rtol=1e-5)  # should be identity

    def getAtoms(self ,  atomHash ,bondHash ,  makeFig):
        atoms = []
        idxList = []
        symbolList = [] 
        for _ , val in atomHash.items(): #break down the hash table into atoms their respective indices and the symbols connected to each atom 
            coords = [float(a) for a in val[2:5]]
            idxList.append(_)
            atoms.append(coords)
            symbolList.append(val[0])
        atoms = np.array(atoms)
        translated = atoms - self.centerPoint
        newAtoms = (self.basisMatrixInv @ translated.T).T

        xMin = -1 * (self.height/2 + self.padding) #alkene bond range 
        xMax = self.height/2 + self.padding #a "padded 0.15 Angstroms is added to resolve edge cases sitting on the bond"

        newXList = newAtoms[:,0]
        newYList = newAtoms[:,1]
        newZList = newAtoms[:,2]
        r = np.sqrt(newYList**2 + newZList**2)
        #print(newXList)
        self.xMin = xMin
        self.xMax = xMax
        newC1 = newAtoms[self.C1Idx - 1]
        newC2 = newAtoms[self.C2Idx - 1]
        self.paramC1 = [xMin , newC1[1], newC1[2]]
        self.paramC2 = [xMax , newC2[1], newC2[2]]

        inside_mask = (newXList >= xMin) & (newXList <= xMax) & (r <= self.radius)
        acceptedAtoms = []
        acceptedIdx = []
        acceptedSymbols = []
        c1c2Idx = [self.C1Idx -1 , self.C2Idx -1 ]

        for i in range (len(inside_mask)):
            idx = idxList[i]
            if inside_mask[i] and idx not in c1c2Idx:
                acceptedAtoms.append(newAtoms[i])
                acceptedIdx.append(idx)
                acceptedSymbols.append(symbolList[i])
        self.acceptedAtoms = acceptedAtoms
        self.acceptedIdx = acceptedIdx
        self.acceptedSymbols = acceptedSymbols 

        #create caps 
        c1Cap = float(self.C1CapHeight + xMax)
        c2Cap = float(-1* self.C2CapHeight + xMin)
        c1End = [c1Cap , newC1[1] , newC1[2] ]
        c2End = [c2Cap , newC2[1] , newC2[2] ]
        c1Int = [xMax , self.radius , newC1[2]]
        c2Int = [xMin , self.radius , newC2[2]]
        self.C1Range = [xMax , c1Cap]
        self.C2Range = [xMin , c2Cap]

        def createCaps(atomList , idxList , symbolList ,localMin, xInt1 , lowerBound , upperBound ):
            print(xInt1[1])
            a = (xInt1[0] -  localMin[0]) / (xInt1[1] - localMin[1])**2 
            xAtoms = np.array(atomList[:,0])
            expRadii = np.sqrt((xAtoms - localMin[0] )/a) + localMin[1] 
            #print("k" , localMin[0])
            funcList = [a , localMin[1] , localMin[0]]
            maskedList = []
            maskedIdx = []
            maskedSymbol = []
            ySpace = np.linspace(localMin[1] , xInt1[0] , 10000)
            xSpace = (ySpace - localMin[0]) / a
            vol = np.pi * np.trapz(xSpace , ySpace)
            for i in range(len(expRadii)):
                radii = expRadii[i]
                atom = atomList[i]
                x = atom[0]
                y = atom[1]
                z = atom[2]
                actRad = np.sqrt(y**2 + z**2)
                if x >= lowerBound and x <= upperBound and actRad <= radii:
                    maskedList.append(atom)
                    maskedIdx.append(idxList[i])
                    maskedSymbol.append(symbolList[i])
            return maskedList , maskedIdx, maskedSymbol , np.abs(vol) , funcList
        c1CapAtoms , c1CapIdx , c1CapSymbols , c1Vol , c1FuncList = createCaps(newAtoms , idxList , symbolList , c1End, c1Int ,c1Int[0] , c1End[0] )
        self.C1FuncList = c1FuncList
        c2CapAtoms , c2CapIdx , c2CapSymbols , c2Vol , c2FuncList = createCaps(newAtoms , idxList , symbolList , c2End, c2Int , c2End[0] , c2Int[0])
        self.C2FuncList = c2FuncList
        acceptedIdx.extend(c1CapIdx)
        acceptedIdx.extend(c2CapIdx)
        self.c1CapAtoms = c1CapAtoms 
        self.c2CapAtoms = c2CapAtoms
        self.c1CapSymbols = c1CapSymbols 
        self.c2CapSymbols = c2CapSymbols
        self.c1CapVol = c1Vol 
        self.c2CapVol = c2Vol
        if makeFig:
            pt = Chem.GetPeriodicTable()
            radiiList = []
            for symbol in symbolList:
                atomRadii = pt.GetRvdw(pt.GetAtomicNumber(symbol))
                radiiList.append(atomRadii)
            figHash = {"atomCoords" : newAtoms , "atomSymbol" : symbolList , "radii" : radiiList}
            stericFigure = stericDrawer(figHash , 30 , highlighted = acceptedIdx)
            for i in range(len(newAtoms)):
                stericFigure.drawAtom(i)
            for _ , bdTable in bondHash.items():
                bdAtoms =  bdTable["idxList"]
                bondType = str(bdTable["bondType"])
                stericFigure.drawBond( bdAtoms[0], bdAtoms[1], bondType)
            self.Fig = stericFigure
    def getBurriedVolume(self, semiCylinders , saveFig):
        if not hasattr(self, "acceptedAtoms"):
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
        volCylinder = np.pi * self.radius**2 * self.height
        if semiCylinders:
            atoms = self.acceptedSymbols
            atomCoords = self.acceptedAtoms

            for i in range(len(self.c1CapAtoms)):
                coords = self.c1CapAtoms[i]
                symbol = self.c1CapSymbols[i]
                atoms.append(symbol)
                atomCoords.append(coords)
            for i in range(len(self.c2CapAtoms)):
                coords = self.c2CapAtoms[i]
                symbol = self.c2CapSymbols[i]
                atoms.append(symbol)
                atomCoords.append(coords)

            atoms = np.array(atoms)
            atomCoords = np.array(atomCoords)

            c1Atoms = np.asarray(self.c1CapSymbols)
            c2Atoms = np.asarray(self.c2CapSymbols)
            c1CapVol = getAtomsVol(c1Atoms)
            c2CapVol = getAtomsVol(c2Atoms)
            c1CapBurr = (c1CapVol / self.c1CapVol)*100  
            c2CapBurr = (c2CapVol / self.c2CapVol)*100
            capBurr = [c1CapBurr , c2CapBurr]     

            if len(atoms) == 0:
                return [0 ,0 , 0, 0 , max(capBurr) , min(capBurr) , (c1CapVol + c2CapVol) / (self.c1CapVol / self.c2CapVol)*100]
            else:
                #pi is the z axis and orth is the y axis
                posIdx_Pi= np.where(atomCoords[:,2] > 0)[0]
                negIdx_Pi = np.where(atomCoords[:,2] < 0)[0]
                posIdx_Orth= np.where(atomCoords[:,1] > 0)[0]
                negIdx_Orth = np.where(atomCoords[:,1] < 0)[0]

                posAtoms_Pi = atoms[posIdx_Pi]
                negAtoms_Pi = atoms[negIdx_Pi]
                posAtoms_Orth = atoms[posIdx_Orth]
                negAtoms_Orth = atoms[negIdx_Orth]

                posVol_Pi = getAtomsVol(posAtoms_Pi)
                negVol_Pi = getAtomsVol(negAtoms_Pi)
                posVol_Orth = getAtomsVol(posAtoms_Orth)
                negVol_Orth = getAtomsVol(negAtoms_Orth)

                posBuried_Pi = (posVol_Pi / ((volCylinder / 2)  + (self.c1CapVol / 2 ) + (self.c2CapVol/2))) * 100
                negBuried_Pi = (negVol_Pi / ((volCylinder / 2)  + (self.c1CapVol / 2 ) + (self.c2CapVol/2))) * 100
                posBuried_Orth = (posVol_Orth / ((volCylinder / 2)  + (self.c1CapVol / 2 ) + (self.c2CapVol/2))) * 100
                negBuried_Orth = (negVol_Orth / ((volCylinder / 2)  + (self.c1CapVol / 2 ) + (self.c2CapVol/2))) * 100
                newXVec = np.array(self.paramC1) - np.array(self.paramC2)
                newYVec = np.array([0,1,0])*self.radius
                newZVec = np.array([0,0,1])*self.radius      
                if saveFig:
                    piPlaneHash = {"shapeType" : "plane" , "C1" : self.paramC1 , "C2" : self.paramC2,
                                    "vector" : newYVec, "reflect" : True , "color" : "blue" }
                    orthPlaneHash = {"shapeType" : "plane" , "C1" : self.paramC1 , "C2" : self.paramC2 , 
                                    "vector" : newZVec, "reflect" : True , "color" : "yellow" }
                    
                    c1SemiHash_Plus = {"shapeType" :  "semiCircle" , "origin" : self.paramC1 , "vector" : newYVec , 
                                        "orthogonal" : newXVec  , "color" : "green"}
                    c1SemiHash_Neg = {"shapeType" :  "semiCircle" , "origin" : self.paramC1 , "vector" : newYVec, 
                                        "orthogonal" : -1 * np.array(newXVec)  , "color" : "red"}
                    c2SemiHash_Plus = {"shapeType" :  "semiCircle" , "origin" : self.paramC2 , "vector" : newYVec, 
                                        "orthogonal" : newXVec , "color" : "green"}
                    c2SemiHash_Neg = {"shapeType" :  "semiCircle" , "origin" : self.paramC2 , "vector" : newYVec , 
                                        "orthogonal" :  -1 * np.array(newXVec)  , "color" : "red"}
                    negPiSemiHash = {"shapeType" :  "semiArc" , "origin" : [0,0,0], "vector" : newYVec , "color" : "red",
                                        "orthogonal" : -1 * np.array(newXVec) , "start" : self.paramC1 , "end" : self.paramC2 } 
                    posPiSemiHash = {"shapeType" :  "semiArc" , "origin" : [0,0,0], "vector" : newYVec, "color" : "green",
                                        "orthogonal" :  np.array(newXVec) , "start" : self.paramC2 , "end" : self.paramC1 } 
                    
                    c1SemiHash_Plus_orth = {"shapeType" :  "semiCircle" , "origin" : self.paramC1 , "vector" :newZVec , 
                                        "orthogonal" : newXVec , "color" : "green"}
                    c1SemiHash_Neg_orth = {"shapeType" :  "semiCircle" , "origin" : self.paramC1 , "vector" :newZVec , 
                                        "orthogonal" : -1 * np.array(newXVec)  , "color" : "red"}
                    c2SemiHash_Plus_orth = {"shapeType" :  "semiCircle" , "origin" : self.paramC2 , "vector" : newZVec , 
                                        "orthogonal" : newXVec  , "color" : "green"}
                    c2SemiHash_Neg_orth = {"shapeType" :  "semiCircle" , "origin" : self.paramC2 , "vector" : newZVec , 
                                        "orthogonal" :  -1 * np.array(newXVec)  , "color" : "red"} 
                    negOrthSemiHash = {"shapeType" :  "semiArc" , "origin" : [0,0,0], "vector" : newZVec  , "color" : "red",
                                        "orthogonal" : -1 * np.array(newXVec) ,"alkeneVec" : -1 * np.array(newXVec) , "start" : self.paramC1, "end" : self.paramC2  } 
                    posOrthSemiHash = {"shapeType" :  "semiArc" , "origin" : [0,0,0] , "vector" : newZVec  , "color" : "green",
                                        "orthogonal" :  np.array(newXVec) ,"alkeneVec" : np.array(newXVec) , "start" : self.paramC1 , "end" : self.paramC2  }    

                    orthFig = self.Fig.drawShapes(orthPlaneHash)  
                    orthFig = self.Fig.drawShapes(c1SemiHash_Plus_orth , fig = orthFig)  
                    orthFig = self.Fig.drawShapes(c1SemiHash_Neg_orth , fig = orthFig)  
                    orthFig = self.Fig.drawShapes(c2SemiHash_Plus_orth ,fig =  orthFig)  
                    orthFig = self.Fig.drawShapes(c2SemiHash_Neg_orth , fig = orthFig) 
                    orthFig = self.Fig.drawShapes(negOrthSemiHash , fig = orthFig) 
                    orthFig = self.Fig.drawShapes(posOrthSemiHash , fig = orthFig)
                    orthFig = self.Fig.drawShapes({ "shapeType" : "cap", "yRange" : self.C1Range , "funcVariables" : self.C1FuncList} , fig = orthFig) 
                    orthFig = self.Fig.drawShapes({"shapeType" : "cap","yRange" : self.C2Range , "funcVariables" : self.C2FuncList} , fig = orthFig) 
                    orthFig.write_html("orthFig.html")

                    piFig = self.Fig.drawShapes(piPlaneHash)  
                    piFig = self.Fig.drawShapes(c1SemiHash_Plus , fig = piFig)  
                    piFig = self.Fig.drawShapes(c1SemiHash_Neg ,fig =  piFig)  
                    piFig = self.Fig.drawShapes(c2SemiHash_Plus , fig = piFig)  
                    piFig = self.Fig.drawShapes(c2SemiHash_Neg ,fig =  piFig) 
                    piFig = self.Fig.drawShapes(negPiSemiHash , fig = piFig) 
                    piFig = self.Fig.drawShapes(posPiSemiHash , fig = piFig)     
                    piFig = self.Fig.drawShapes({"shapeType" : "cap","yRange" : self.C1Range , "funcVariables" : self.C1FuncList} , fig = piFig) 
                    piFig = self.Fig.drawShapes({"shapeType" : "cap","yRange" : self.C2Range , "funcVariables" : self.C2FuncList} , fig = piFig)   
                    piFig.write_html("piFig.html")
                    
                    c1Burr = {"shapeType" :"Sphere" , "radius" : self.radius , "coordinates" : self.paramC1 , "color" : "green"}
                    vBurrFig = self.Fig.drawShapes(c1Burr) 
                    c2Burr = {"shapeType" :"Sphere" , "radius" : self.radius , "coordinates" : self.paramC2 , "color" : "red"}
                    vBurrFig = self.Fig.drawShapes(c2Burr , fig = vBurrFig) 
                    vBurrFig.write_html("vBurrFig.html")
                    '''
                    from plotly_gif import GIF, three_d_scatter_rotate
                    gifVburr = GIF()
                    three_d_scatter_rotate(gifVburr, vBurrFig)
                    gifOrth = GIF()
                    three_d_scatter_rotate(gifOrth, orthFig)
                    gifPi = GIF()
                    three_d_scatter_rotate(gifPi, piFig)
                    
                    gifPi.create_gif()
                    gifVburr.create_gif()
                    gifOrth.create_gif()
                    '''
                totBurrVol = ((posVol_Pi + negVol_Pi)/(volCylinder + self.c1CapVol + self.c2CapVol)) * 100 
                return [max(posBuried_Pi, negBuried_Pi), min(posBuried_Pi, negBuried_Pi) , 
                        max(posBuried_Orth, negBuried_Orth), min(posBuried_Orth, negBuried_Orth),
                        max(capBurr) , min(capBurr) , totBurrVol]    
        else:
            c1Atoms = np.asarray(self.c1CapSymbols)
            c2Atoms = np.asarray(self.c2CapSymbols)
            c1CapVol = getAtomsVol(c1Atoms)
            c2CapVol = getAtomsVol(c2Atoms)
            atoms = list(self.acceptedSymbols)
            atomVol = getAtomsVol(atoms)
            percentBurriedVol = ((atomVol + c1CapVol + c2CapVol)/(volCylinder + self.c1CapVol + self.c2CapVol)) * 100 
            self.BuriedCylinder = percentBurriedVol
            return [percentBurriedVol]
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
    for atomIdx , coords in coordHash.items():
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
    mainCylinder = alkeneSemiCylinders(C1Hash , C2Hash , radius , 0.15)
    print(C1Hash , C2Hash)
    mainCylinder.getAtoms(coordHash ,bondHash, True)
    maxSemi_Pi , minSemi_Pi , maxSemi_Orth, minSemi_Orth , maxCap , minCap ,totalCap = mainCylinder.getBurriedVolume(True , True)
    print(maxSemi_Pi , minSemi_Pi , maxSemi_Orth, minSemi_Orth , maxCap , minCap , totalCap)

if __name__ == "__main__":


    logFile = str(sys.argv[1])
    smilesStr = str(sys.argv[2])
    radius = float(sys.argv[3])
    linkIdx = int(sys.argv[4])
    main(logFile , smilesStr , radius , linkIdx)
    