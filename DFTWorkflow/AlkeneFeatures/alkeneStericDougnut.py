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
    print("mask" , inside_mask)
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
    
    print(f"Total atoms: {total}")
    print(f"Atoms inside cylinder: {inside}")
    print(f"Atoms outside cylinder: {outside}")
    print(f"Percentage inside: {100*inside/total:.1f}%")
    
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
    print("insiderList" , insideList)
    volDoughnut = np.pi*outerR**2*h - np.pi*innerR**2*h
    print("volDoughnut" , volDoughnut)
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
    def __init__(self, C1Hash, C2Hash,radius  ):
        self.C1Idx = int(C1Hash["idx"]) #Index of alkene Carbon 1
        self.C2Idx = int(C2Hash["idx"]) #Index of alkene Carbon 2
        self.C1 =  C1Hash["0"]
        self.C2 =  C2Hash["0"]

        alkeneVec = C1Hash["0"] - C2Hash["0"] #vector between C1 and C2
        self.alkeneVec = alkeneVec 
        self.height = np.linalg.norm(C1Hash["0"] - C2Hash["0"])#alkeneBondLength 
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

        # plane normals
        c1Orth = np.cross(c1Vec1, c1Vec2)
        c2Orth = np.cross(c2Vec1, c2Vec2)

        # normalize
        c1Orth /= np.linalg.norm(c1Orth)
        c2Orth /= np.linalg.norm(c2Orth)

        # ALIGN DIRECTIONS
        if np.dot(c1Orth, c2Orth) < 0:
            c2Orth *= -1

        # scale
        c1OrthScaled = radius * c1Orth
        c2OrthScaled = radius * c2Orth
        center = [(float(a) + float(b)) * 0.5 for a, b in zip(list(C1Hash["0"]), list(C2Hash["0"]))]
        self.centerPoint = center 
        bisector = (c1OrthScaled + c2OrthScaled )
        bisector /= np.linalg.norm(c1OrthScaled + c2OrthScaled )
        bisector *= radius
        self.bisector = bisector
        self.radius = radius

    def getAtoms(self ,  atomHash ,bondHash ,  makeFig):
        atoms = []
        idxList = []
        symbolList = [] 
        for _ , val in atomHash.items():
            coords = [float(a) for a in val[2:5]]
            idx = val[1]
            if idx in [self.C1Idx,self.C2Idx]:
                pass
            else:
                idxList.append(idx)
                atoms.append(coords)
                symbolList.append(val[0])
        if makeFig:
            pt = Chem.GetPeriodicTable()
            radiiList = []
            for symbol in symbolList:
                atomRadii = pt.GetRvdw(pt.GetAtomicNumber(symbol))
                radiiList.append(atomRadii)
            figHash = {"atomCoords" : atoms , "atomSymbol" : symbolList , "radii" : radiiList}
            stericFigure = stericDrawer(figHash , 30)
            for i in range(len(atoms)):
                stericFigure.drawAtom(i)
            for _ , bdTable in bondHash.items():
                bdAtoms =  bdTable["idxList"]
                bondType = str(bdTable["bondType"])
                stericFigure.drawBond( bdAtoms[0], bdAtoms[1], bondType)

            self.Fig = stericFigure
        atoms = np.array(atoms)
        #print(atoms)
        cx, cy, cz = self.centerPoint
        xList = atoms[:, 0] - cx
        yList = atoms[:, 1] - cy
        zList = atoms[:, 2] - cz

        # Calculate radial distance from cylinder axis (z-axis)
        r = np.sqrt(xList**2 + yList**2)

        inside_mask = (r <= self.radius)  & (zList <= self.height)
        acceptedAtoms = []
        acceptedIdx = []
        acceptedSymbols = []
        for i in range (len(inside_mask)):
            if inside_mask[i]:
                acceptedAtoms.append(atoms[i])
                acceptedIdx.append(idxList[i])
                acceptedSymbols.append(symbolList[i])
        self.acceptedAtoms = acceptedAtoms
        self.acceptedIdx = acceptedIdx
        self.acceptedSymbols = acceptedSymbols 

    def getBurriedVolume(self, semiSpheres , saveFig):
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
        if semiSpheres:
            tol = 1e-8
            atoms = np.asarray(self.acceptedSymbols)
            atomCoords = np.asarray(self.acceptedAtoms)

            planeNormal = np.cross(self.alkeneVec, self.bisector)
            n = planeNormal / np.linalg.norm(planeNormal)
            n *= self.radius 
            signedDist_Pi = np.dot(atomCoords - self.centerPoint, n)
            signedDist_Orth = np.dot(atomCoords - self.centerPoint, self.bisector)

            posIdx_Pi= np.where(signedDist_Pi > tol)[0]
            negIdx_Pi = np.where(signedDist_Pi < -tol)[0]
            posIdx_Orth= np.where(signedDist_Orth > tol)[0]
            negIdx_Orth = np.where(signedDist_Orth < -tol)[0]

            posAtoms_Pi = atoms[posIdx_Pi]
            negAtoms_Pi = atoms[negIdx_Pi]
            posAtoms_Orth = atoms[posIdx_Orth]
            negAtoms_Orth = atoms[negIdx_Orth]

            posVol_Pi = getAtomsVol(posAtoms_Pi)
            negVol_Pi = getAtomsVol(negAtoms_Pi)
            posVol_Orth = getAtomsVol(posAtoms_Orth)
            negVol_Orth = getAtomsVol(negAtoms_Orth)

            posBuried_Pi = (posVol_Pi / (volCylinder / 2)) * 100
            negBuried_Pi = (negVol_Pi / (volCylinder / 2)) * 100
            posBuried_Orth = (posVol_Orth / (volCylinder / 2)) * 100
            negBuried_Orth = (negVol_Orth / (volCylinder / 2)) * 100
            if saveFig:
                piPlaneHash = {"shapeType" : "plane" , "C1" : self.C1 , "C2" : self.C2 ,
                                "vector" : self.bisector, "reflect" : True , "color" : "blue" }
                orthPlaneHash = {"shapeType" : "plane" , "C1" : self.C1 , "C2" : self.C2 , 
                                 "vector" : n, "reflect" : True , "color" : "yellow" }
                
                c1SemiHash_Plus = {"shapeType" :  "semiCircle" , "origin" : self.C1 , "vector" : self.bisector , 
                                    "orthogonal" : self.alkeneVec  , "color" : "green"}
                c1SemiHash_Neg = {"shapeType" :  "semiCircle" , "origin" : self.C1 , "vector" : self.bisector , 
                                    "orthogonal" : -1 * np.array(self.alkeneVec)  , "color" : "red"}
                c2SemiHash_Plus = {"shapeType" :  "semiCircle" , "origin" : self.C2 , "vector" : self.bisector , 
                                    "orthogonal" : self.alkeneVec  , "color" : "green"}
                c2SemiHash_Neg = {"shapeType" :  "semiCircle" , "origin" : self.C2 , "vector" : self.bisector , 
                                    "orthogonal" :  -1 * np.array(self.alkeneVec)  , "color" : "red"}
                negPiSemiHash = {"shapeType" :  "semiArc" , "origin" : self.centerPoint , "radius" : self.radius , 
                                    "length" : float(np.linalg.norm(self.alkeneVec)) , "vector" : self.bisector , "color" : "red",
                                    "orthogonal" : -1 * np.array(self.alkeneVec) ,"alkeneVec" : -1 * np.array(self.alkeneVec)  } 
                posPiSemiHash = {"shapeType" :  "semiArc" , "origin" : self.centerPoint , "radius" : self.radius , 
                                    "length" : float(np.linalg.norm(self.alkeneVec)) , "vector" : self.bisector , "color" : "green",
                                    "orthogonal" :  np.array(self.alkeneVec) ,"alkeneVec" : np.array(self.alkeneVec)  } 
                
                c1SemiHash_Plus_orth = {"shapeType" :  "semiCircle" , "origin" : self.C1 , "vector" : n, 
                                    "orthogonal" : self.alkeneVec  , "color" : "green"}
                c1SemiHash_Neg_orth = {"shapeType" :  "semiCircle" , "origin" : self.C1 , "vector" : n, 
                                    "orthogonal" : -1 * np.array(self.alkeneVec)  , "color" : "red"}
                c2SemiHash_Plus_orth = {"shapeType" :  "semiCircle" , "origin" : self.C2 , "vector" : n , 
                                    "orthogonal" : self.alkeneVec  , "color" : "green"}
                c2SemiHash_Neg_orth = {"shapeType" :  "semiCircle" , "origin" : self.C2 , "vector" : n , 
                                    "orthogonal" :  -1 * np.array(self.alkeneVec)  , "color" : "red"} 
                negOrthSemiHash = {"shapeType" :  "semiArc" , "origin" : self.centerPoint , "radius" : self.radius , 
                                    "length" : float(np.linalg.norm(self.alkeneVec)) , "vector" : n , "color" : "red",
                                    "orthogonal" : -1 * np.array(self.alkeneVec) ,"alkeneVec" : -1 * np.array(self.alkeneVec)  } 
                posOrthSemiHash = {"shapeType" :  "semiArc" , "origin" : self.centerPoint , "radius" : self.radius , 
                                    "length" : float(np.linalg.norm(self.alkeneVec)) , "vector" : n , "color" : "green",
                                    "orthogonal" :  np.array(self.alkeneVec) ,"alkeneVec" : np.array(self.alkeneVec)  }    

                orthFig = self.Fig.drawShapes(orthPlaneHash)  
                orthFig = self.Fig.drawShapes(c1SemiHash_Plus_orth , fig = orthFig)  
                orthFig = self.Fig.drawShapes(c1SemiHash_Neg_orth , fig = orthFig)  
                orthFig = self.Fig.drawShapes(c2SemiHash_Plus_orth ,fig =  orthFig)  
                orthFig = self.Fig.drawShapes(c2SemiHash_Neg_orth , fig = orthFig) 
                orthFig = self.Fig.drawShapes(negOrthSemiHash , fig = orthFig) 
                orthFig = self.Fig.drawShapes(posOrthSemiHash , fig = orthFig) 
                orthFig.write_html("orthFig.html")

                piFig = self.Fig.drawShapes(piPlaneHash)  
                piFig = self.Fig.drawShapes(c1SemiHash_Plus , fig = piFig)  
                piFig = self.Fig.drawShapes(c1SemiHash_Neg ,fig =  piFig)  
                piFig = self.Fig.drawShapes(c2SemiHash_Plus , fig = piFig)  
                piFig = self.Fig.drawShapes(c2SemiHash_Neg ,fig =  piFig) 
                piFig = self.Fig.drawShapes(negPiSemiHash , fig = piFig) 
                piFig = self.Fig.drawShapes(posPiSemiHash , fig = piFig)       
                piFig.write_html("piFig.html")
                
                c1Burr = {"shapeType" :"Sphere" , "radius" : 3.0 , "coordinates" : self.C1 , "color" : "green"}
                vBurrFig = self.Fig.drawShapes(c1Burr) 
                c2Burr = {"shapeType" :"Sphere" , "radius" : 3.0 , "coordinates" : self.C2 , "color" : "red"}
                vBurrFig = self.Fig.drawShapes(c2Burr , fig = vBurrFig) 
                vBurrFig.write_html("vBurrFig.html")
            return [max(posBuried_Pi, negBuried_Pi), min(posBuried_Pi, negBuried_Pi) , 
                    max(posBuried_Orth, negBuried_Orth), min(posBuried_Orth, negBuried_Orth)]    
        else:
            atoms = list(self.acceptedSymbols)
            atomVol = getAtomsVol(atoms)
            percentBurriedVol = (atomVol/volCylinder) * 100 
            self.BuriedCylinder = percentBurriedVol
            return [percentBurriedVol]
def main(logFile , smilesStr, radius , linkIdx):
    from DFTWorkflow.AlkeneFeatures.alkeneFeatureExtractorMAST import getAtomCoordsRobust
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
    #print(coordHash)
    idx = 0
    for _ , coords in coordHash.items():
        idx +=1
        if idx == C1:
            c1_coords = np.array(coords[2:5])
            c1Idx = idx
        elif idx == C2:
            c2_coords = coords[2:5]
            c2Idx = idx 
        if (idx-1) in CminNeighbors:
            crds = np.array(coords[2:5])
            CminContacts.append(crds) 
        if (idx-1) in CmaxNeighbors:
            crds = np.array(coords[2:5])
            CmaxContacts.append(crds) 
    C1Hash = {"0" : c1_coords , "1" : CminContacts , "idx" : c1Idx}
    C2Hash = {"0" : c2_coords , "1" : CmaxContacts , "idx" : c2Idx}
    mainCylinder = alkeneSemiCylinders(C1Hash , C2Hash , radius)
    mainCylinder.getAtoms(coordHash ,bondHash, True)
    maxSemi_Pi , minSemi_Pi , maxSemi_Orth, minSemi_Orth = mainCylinder.getBurriedVolume(True , True)

if __name__ == "__main__":


    logFile = str(sys.argv[1])
    smilesStr = str(sys.argv[2])
    radius = float(sys.argv[3])
    linkIdx = int(sys.argv[4])
    main(logFile , smilesStr , radius , linkIdx)
    