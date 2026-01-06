#The steric doughnut | Danny Ruiz de Castilla
#Not really a doughnut
import numpy as np
import plotly.graph_objects as go
import sys
from pathlib import Path
from rdkit import Chem
parentDir = Path(__file__).resolve().parents[2]
sys.path.append(str(parentDir))
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
        self.C1Idx = int(C1Hash["idx"])
        self.C2Idx = int(C2Hash["idx"])
        #print("228" , C1Hash["0"])
        alkeneVec = C1Hash["0"] - C2Hash["0"]
        self.alkeneVec = alkeneVec 
        self.height = np.linalg.norm(C1Hash["0"] - C2Hash["0"])#alkeneBondLength 
        # vectors from carbons to neighbors
        c1Contacts = C1Hash["1"]
        c2Contacts = C2Hash["1"]
        c1Vec1 = c1Contacts[0] - C1Hash["0"]
        c1Vec2 = c1Contacts[1] - C1Hash["0"]
        c2Vec1 = c2Contacts[0] - C2Hash["0"]
        c2Vec2 = c2Contacts[1] - C2Hash["0"]

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
        self.bisector = c1OrthScaled + c2OrthScaled
        self.radius = radius

    def getAtoms(self ,  atomHash):
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
        atoms = np.array(atoms)
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

    def getBurriedVolume(self, semiSpheres):
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
        if semiSpheres:
            tol = 1e-8
            volCylinder = np.pi * self.radius**2 * self.height

            atoms = np.asarray(self.acceptedSymbols)
            atomCoords = np.asarray(self.acceptedAtoms)

            planeNormal = np.cross(self.alkeneVec, self.bisector)
            n = planeNormal / np.linalg.norm(planeNormal)

            signedDist = np.dot(atomCoords - self.centerPoint, n)

            pos_idx = np.where(signedDist > tol)[0]
            neg_idx = np.where(signedDist < -tol)[0]

            posAtoms = atoms[pos_idx]
            negAtoms = atoms[neg_idx]

            posVol = getAtomsVol(posAtoms)
            negVol = getAtomsVol(negAtoms)

            posBuried = (posVol / (volCylinder / 2)) * 100
            negBuried = (negVol / (volCylinder / 2)) * 100

            return [max(posBuried, negBuried), min(posBuried, negBuried)]
        else:
            atoms = list(self.acceptedSymbols)
            atomVol = getAtomsVol(atoms)
            volCylinder = np.pi*self.radius**2*self.height
            percentBurriedVol = (atomVol/volCylinder) * 100 
            self.BuriedCylinder = percentBurriedVol
            return [percentBurriedVol]
if __name__ == "__main__":
    # Define cylinder parameters
    inner_r = 2.0
    outer_r = 5.0
    h = 8.0
    cylinder_center = (0, 0, 0)
    
    # Create sample atoms (random coordinates)
    np.random.seed(42)
    n_atoms = 100
    atoms = np.random.uniform(-7, 7, size=(n_atoms, 3))
    atoms[:, 2] = np.random.uniform(-2, 10, size=n_atoms)  # z coordinates
    
    print("Sample atom coordinates (first 5):")
    print(atoms[:5])
    print()
    
    # Check which atoms are inside the cylinder
    inside_atoms = check_atoms_in_cylinder(atoms, inner_r, outer_r, h, cylinder_center)
    
    print("Atoms inside the cylinder:")
    for idx, coord in inside_atoms.items():
        print(f"  Atom {idx}: {coord}")
    print()
    
    # Get summary
    get_atoms_summary(inside_atoms, atoms)
    print()
    
    # Visualize with Plotly (interactive!)
    plot_cylinder_with_atoms(atoms, inside_atoms, inner_r, outer_r, h, 
                            cylinder_center, resolution=50)