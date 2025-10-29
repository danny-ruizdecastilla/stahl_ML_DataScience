#The steric doughnut | Danny Ruiz de Castilla
#Not really a doughnut
import numpy as np
import plotly.graph_objects as go

def create_hollow_cylinder(inner_radius, outer_radius, height, resolution=50):
    """
    Create a 3D hollow cylinder (doughnut shape).
    
    Parameters:
    -----------
    inner_radius : float
        Inner radius of the hollow cylinder
    outer_radius : float
        Outer radius of the hollow cylinder
    height : float
        Height of the cylinder
    resolution : int
        Number of points for mesh resolution (default: 50)
    
    Returns:
    --------
    dict containing the surfaces: outer_wall, inner_wall, top_ring, bottom_ring
    """
    
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

def check_atoms_in_cylinder(atoms, inner_radius, outer_radius, height, center=(0, 0, 0)):
    """
    Check which atoms are inside the hollow cylinder and return them in a hash table.
    
    Parameters:
    -----------
    atoms : numpy array or list
        Array of shape (N, 3) containing [x, y, z] coordinates of atoms
    inner_radius : float
        Inner radius of the hollow cylinder
    outer_radius : float
        Outer radius of the hollow cylinder
    height : float
        Height of the cylinder
    center : tuple
        Center coordinates (cx, cy, cz) of the cylinder base (default: (0, 0, 0))
    
    Returns:
    --------
    dict : Hash table with atom indices as keys and coordinates as values
        Example: {0: [x, y, z], 5: [x, y, z], ...}
    """
    atoms = np.array(atoms)
    
    if atoms.ndim != 2 or atoms.shape[1] != 3:
        raise ValueError("Atoms array must be of shape (N, 3)")
    
    cx, cy, cz = center
    
    # Translate atoms relative to cylinder center
    x = atoms[:, 0] - cx
    y = atoms[:, 1] - cy
    z = atoms[:, 2] - cz
    
    # Calculate radial distance from cylinder axis (z-axis)
    r = np.sqrt(x**2 + y**2)
    
    # Check conditions for being inside the hollow cylinder:
    # 1. Radial distance between inner and outer radius
    # 2. Height between 0 and cylinder height
    inside_mask = (r >= inner_radius) & (r <= outer_radius) & (z >= 0) & (z <= height)
    
    # Create hash table with indices of atoms inside
    inside_atoms = {}
    inside_indices = np.where(inside_mask)[0]
    
    for idx in inside_indices:
        inside_atoms[int(idx)] = atoms[idx].tolist()
    
    return inside_atoms

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
                             height, center=(0, 0, 0), resolution=50):
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
        inside_coords = np.array(list(inside_atoms_dict.values()))
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

# Example usage
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