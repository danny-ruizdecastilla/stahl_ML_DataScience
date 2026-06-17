import plotly.graph_objects as go
import plotly 
import pandas as pd
import os 
import sys
import glob
import base64
import numpy as np
import json
from pathlib import Path
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from figs.repositoryGraphVisualization import str_to_rgb
atomColorHashRGB = { "H" : "rgba(220,220,220,0.7)" , "C" : "rgba(128,128,128,0.7)"
                    , "O" : "rgba(197,5,12,0.7)" , "N" : "rgba(001,031,091,0.7)"
                     , "S" : "rgba(254,242,80,0.7)" , "Si" : "rgba(170,089,043,0.7)"
                      , "P" : "rgba(245,128,037,0.7)" , "F" : "rgba(138,206,0,0.7)"
                       , "Br" :  "rgba(165,028,048,0.7)" , "Cl" :"rgba(255,206,207,0.7)" 
                       , "B" : "rgba(255,192,203)"
                         }
class stericDrawer:
    def __init__(self, atomsHash , resolution, **kwargs):
        self.atomList = atomsHash["atomCoords"]
        self.atomSymbols = atomsHash["atomSymbol"]
        self.atomRadii = atomsHash["radii"]
        if kwargs.get("highlighted") is None:
            pass
        else:
            self.highlighted = kwargs["highlighted"]
        fig = go.Figure()
        fig.update_layout(
        paper_bgcolor="white",
        plot_bgcolor="white",
        scene=dict(
            xaxis=dict(
                visible=False,
                showbackground=False,
                showgrid=False,
                zeroline=False
            ),
            yaxis=dict(
                visible=False,
                showbackground=False,
                showgrid=False,
                zeroline=False
            ),
            zaxis=dict(
                visible=False,
                showbackground=False,
                showgrid=False,
                zeroline=False
            ),
            bgcolor="white"
        ),
        margin=dict(l=0, r=0, t=0, b=0)
    )
        self.Figure = fig
        self.resolution = resolution
    def drawAtom(self ,idx,):
      theta = np.linspace(0, 2 * np.pi, self.resolution)  # azimuthal angle
      phi = np.linspace(0, np.pi, self.resolution)        # polar angle
      theta_grid, phi_grid = np.meshgrid(theta, phi)
      atomStr = self.atomSymbols[idx]
      if hasattr(self , "highlighted") and idx in self.highlighted:
        atomColor = "rgba(100,31.8,0,0.8)"
        #print("highlighted")
      else:
         atomColor = atomColorHashRGB[atomStr]
      radius = float(0.5* self.atomRadii[idx])
      atomCoord = self.atomList[idx]
      x_outer = radius * np.sin(phi_grid) * np.cos(theta_grid)
      y_outer = radius * np.sin(phi_grid) * np.sin(theta_grid)
      z_outer = radius * np.cos(phi_grid)
      self.Figure.add_trace(go.Surface(
          x=x_outer + atomCoord[0], y=y_outer + atomCoord[1], z=z_outer + atomCoord[2],
          opacity=0.75,colorscale=[[0, atomColor], [1, atomColor]],showscale=False,name=f'Atom_{idx}',
          hovertemplate=f'<b>Atom_{idx}</b><br>' +
                  'x: %{x:.2f}<br>' +
                  'y: %{y:.2f}<br>' +
                  'z: %{z:.2f}<br>' +
                  '<extra></extra>',))
    def drawBond(self, idx1, idx2, bondStr):
        atom1 = self.atomList[idx1]
        atom2 = self.atomList[idx2]
        if bondStr == "SINGLE":
          w = 5
        elif bondStr == "DOUBLE":
          w = 15
        elif bondStr == "TRIPLE":
          w = 30    
        elif bondStr == "AROMATIC":
          w = 10
        self.Figure.add_trace(go.Scatter3d(
            x=[atom1[0], atom2[0]],
            y=[atom1[1], atom2[1]],
            z=[atom1[2], atom2[2]],
            mode='lines',
            line=dict(
                color='rgba(0,0,0,0.5)',
                width=w
            ),
            name=f'Bond_{idx1}-{idx2}',
            hovertemplate=f'<b>{bondStr}</b><br>' +
                          f'Atom {idx1} ↔ Atom {idx2}<br>' +
                          '<extra></extra>',
            showlegend=False
        ))
    def drawShapes(self , shapeHash , **kwargs):
      fig = kwargs.get("fig", self.Figure)
      typeStr = shapeHash["shapeType"]
      if typeStr == "line":
        origin = shapeHash["origin"]
        vector = shapeHash["vector"]
        fig.add_trace(go.Scatter3d(
            x=[origin[0], vector[0]],
            y=[origin[1], vector[1]],
            z=[origin[2], vector[2]],
            mode='lines',
            line=dict(
                color=shapeHash["color"],
                width=5
            ),
            name=f'line_{shapeHash["name"]}',
            hovertemplate=f'<b>line_{shapeHash["name"]}</b><br>' +
                          '<extra></extra>',
            showlegend=False))
      elif typeStr == "semiCircle":
        newFig = go.Figure(fig.to_dict())
        resolution = self.resolution 
        origin = shapeHash["origin"]
        vector = shapeHash["vector"]
        orthogonal = shapeHash["orthogonal"]
        colorStr = shapeHash["color"]
        rgb = (*str_to_rgb(colorStr) , 0.35)
        r,g,b,a = rgb
        rgba_str = f'rgba({r},{g},{b},{a})'

        R = np.linalg.norm(vector)
        u_hat = vector / R

        v_hat = np.cross(u_hat , orthogonal)
        v_hat /= np.linalg.norm(v_hat)

        r = np.linspace(0,R , resolution)
        theta = np.linspace(0,np.pi,resolution)
        Rg,Tg = np.meshgrid(r, theta)

        X = origin[0] + Rg*np.cos(Tg) * u_hat[0] + Rg*np.sin(Tg) * v_hat[0]
        Y = origin[1] + Rg*np.cos(Tg) * u_hat[1] + Rg*np.sin(Tg) * v_hat[1]
        Z = origin[2] + Rg*np.cos(Tg) * u_hat[2] + Rg*np.sin(Tg) * v_hat[2]

        newFig.add_trace(go.Surface(
          x = X,
          y = Y,
          z = Z,
          showscale = False,
          colorscale = [[0 , rgba_str] , [1 , rgba_str]],
          name = "Semicircle Face"
        ))
        return newFig
      elif typeStr == "Sphere":
        newFig = go.Figure(fig.to_dict())
        theta = np.linspace(0, 2 * np.pi, self.resolution)  # azimuthal angle
        phi = np.linspace(0, np.pi, self.resolution)        # polar angle
        theta_grid, phi_grid = np.meshgrid(theta, phi)
        radius = shapeHash["radius"]
        atomCoord = shapeHash["coordinates"]
        colorStr = shapeHash["color"]
        rgb = (*str_to_rgb(colorStr), 0.35)
        r, g, b, a = rgb
        rgba_str = f'rgba({r}, {g}, {b}, {a})'
        x_outer = radius * np.sin(phi_grid) * np.cos(theta_grid)
        y_outer = radius * np.sin(phi_grid) * np.sin(theta_grid)
        z_outer = radius * np.cos(phi_grid)
        newFig.add_trace(go.Surface(
            x=x_outer + atomCoord[0], y=y_outer + atomCoord[1], z=z_outer + atomCoord[2],
            opacity=0.6,colorscale=[[0, rgba_str], [1, rgba_str]],showscale=False,name=f"Vburr_{radius}",
            hovertemplate=f'<b>Vburr_{radius}</b><br>' +
                    'x: %{x:.2f}<br>' +
                    'y: %{y:.2f}<br>' +
                    'z: %{z:.2f}<br>' +
                    '<extra></extra>',))
        return newFig
      elif typeStr == "slicedOrange":
        newFig = go.Figure(fig.to_dict())
        rad = shapeHash["radius"]
        atomCoord = shapeHash["atom"]
        inverseMatrix = shapeHash["inverseCoB"]
        flipped = shapeHash["flipped"]
        theta = np.linspace(0, 2*np.pi, self.resolution)
        phi = np.linspace(0, np.pi, self.resolution)
        theta_grid, phi_grid = np.meshgrid(theta, phi)
        x = rad * np.sin(phi_grid) * np.cos(theta_grid)
        y = rad * np.sin(phi_grid) * np.sin(theta_grid)
        z = rad * np.cos(phi_grid)
        if flipped:
           region = ((x <= 0).astype(int)+ 2 * (z <= 0).astype(int))
        else:
            region = ((x >= 0).astype(int)+ 2 * (z >= 0).astype(int)) #gives 4 distinct assignments 

        pts_local = np.stack(
            [
                x.ravel(),
                y.ravel(),
                z.ravel()
            ],
            axis=1
        )

        pts_global = (
            pts_local @ inverseMatrix
            + atomCoord
        )

        xGlob = pts_global[:, 0].reshape(x.shape)
        yGlob = pts_global[:, 1].reshape(y.shape)
        zGlob = pts_global[:, 2].reshape(z.shape)

        colorscale = [
            [0.00, "red"],
            [0.25, "red"],

            [0.25, "blue"],
            [0.50, "blue"],

            [0.50, "green"],
            [0.75, "green"],

            [0.75, "orange"],
            [1.00, "orange"],
        ]

        newFig.add_trace(
            go.Surface(
                x=xGlob,
                y=yGlob,
                z=zGlob,
                surfacecolor=region,
                cmin=0,
                cmax=3,
                colorscale=colorscale,
                showscale=False,
                opacity=0.3,
                name=f"PartitionedSphere_{rad}",
                hovertemplate=
                    "<b>Region %{surfacecolor}</b><br>"
                    "x: %{x:.2f}<br>"
                    "y: %{y:.2f}<br>"
                    "z: %{z:.2f}<br>"
                    "<extra></extra>"
            )
        )
        return newFig
      elif typeStr == "plane":
        newFig = go.Figure(fig.to_dict())
        C1 = shapeHash["C1"]
        C2 = shapeHash["C2"]
        vector = shapeHash["vector"]
        resolution = self.resolution
        reflect = shapeHash["reflect"]
        colorStr = shapeHash["color"]
        rgb = (*str_to_rgb(colorStr), 0.8)
        r, g, b, a = rgb
        rgba_str = f'rgba({r}, {g}, {b}, {a})'

        if reflect:
            P00 = C1 + vector
            P10 = C2 + vector
            P01 = C1 - vector
            P11 = C2 - vector
        else:
            P00 = C1
            P10 = C2
            P01 = C1 + vector
            P11 = C2 + vector

        # Parameter space
        s = np.linspace(0, 1, resolution)
        t = np.linspace(0, 1, resolution)
        S, T = np.meshgrid(s, t)

        # Bilinear interpolation
        X = ((1 - S) * (1 - T) * P00[0] + S * (1 - T) * P10[0] + (1 - S) * T * P01[0] + S * T * P11[0])
        Y = ((1 - S) * (1 - T) * P00[1] + S * (1 - T) * P10[1] + (1 - S) * T * P01[1] + S * T * P11[1])
        Z = ((1 - S) * (1 - T) * P00[2] + S * (1 - T) * P10[2] + (1 - S) * T * P01[2] + S * T * P11[2])
        newFig.add_trace(go.Surface(
            x=X,
            y=Y,
            z=Z,
            opacity=0.3,
            colorscale=[[0, rgba_str], [1, rgba_str]],
            showscale=False,
            hoverinfo='skip'))
        return newFig
      elif typeStr == "cap":
        newFig = go.Figure(fig.to_dict())

        resolution = self.resolution

        # y-range for the cap surface
        yMin  = shapeHash["yRange"][0]
        yMax  = shapeHash["yRange"][1]

        a = shapeHash["funcVariables"][0]
        h = shapeHash["funcVariables"][1]
        k = shapeHash["funcVariables"][2]
        # sample along y
        y_vals = np.linspace(yMin, yMax, resolution)

        theta = np.linspace(0, 2*np.pi, resolution)

        X = []
        Y = []
        Z = []

        for y in y_vals:

            r = np.sqrt((y - k )/a) + h


            # rotate around x-axis
            x_ring = np.full_like(theta, y)
            y_ring = r * np.cos(theta)
            z_ring = r * np.sin(theta)

            X.append(x_ring)
            Y.append(y_ring)
            Z.append(z_ring)

        X = np.array(X)
        Y = np.array(Y)
        Z = np.array(Z)

        newFig.add_trace(go.Surface(
            x=X,
            y=Y,
            z=Z,
            opacity=0.3,
            showscale=False,
            name="cap"
        ))
        return newFig

      elif typeStr == "semiArc":
        newFig = go.Figure(fig.to_dict())
        resolution = self.resolution 
        origin = np.asarray(shapeHash["origin"], dtype=float)
        vector = shapeHash["vector"] #vector where the arcs span
        colorStr = shapeHash["color"]
        point1 = np.array(shapeHash["start"]) #coordinates of the lower bound of the arcs 
        point2 = np.array(shapeHash["end"]) #coordiantes of the upper bound of the arcs
        orthogonal = shapeHash["orthogonal"]
        rgb = (*str_to_rgb(colorStr) , 0.35)
        r,g,b,a = rgb
        rgba_str = f'rgba({r},{g},{b},{a})'

        R = np.linalg.norm(vector)
        u_hat = vector / R

        v_hat = np.cross(u_hat , orthogonal)
        v_hat /= np.linalg.norm(v_hat)

        theta = np.linspace(0, np.pi, resolution)

        x_local = R * np.cos(theta)
        y_local = R * np.sin(theta)

        points = (
            point1[None, :]
            + x_local[:, None] * u_hat[None, :]
            + y_local[:, None] * v_hat[None, :]
        )

        arc_grid = points[:, None, :]

        extrude_vec = np.array(point2) - np.array(point1)
        L = np.linalg.norm(extrude_vec)

        extrude_hat = extrude_vec / L

        s = np.linspace(0, L, resolution)
        s_grid = s[None, :, None]

        surface = arc_grid + s_grid * extrude_hat[None, None, :]

        X = surface[:, :, 0]
        Y = surface[:, :, 1]
        Z = surface[:, :, 2]
        newFig.add_surface(
            x=X,
            y=Y,
            z=Z,
            showscale=False,
            opacity=a,
            colorscale=[[0, rgba_str], [1, rgba_str]],
            lighting=dict(
                ambient=0.6,
                diffuse=0.6,
                specular=0.1,
                roughness=0.9
            )
        )
        return newFig
    def saveHTML(self , figDir  , figName):
       from pathlib import Path
       newOutputDir = Path(figDir)
       newOutputDir.mkdir(parents=True, exist_ok=True)
       newFig = go.Figure(self.Figure.to_dict())
       newFig.write_html(newOutputDir / f"{figName}.html")


