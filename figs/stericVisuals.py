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
                         }
class stericDrawer:
    def __init__(self, atomsHash , resolution):
        self.atomList = atomsHash["atomCoords"]
        self.atomSymbols = atomsHash["atomSymbol"]
        self.atomRadii = atomsHash["radii"]
        fig = go.Figure()
        self.Figure = fig
        self.resolution = resolution
    def drawAtom(self ,idx,):
      theta = np.linspace(0, 2 * np.pi, self.resolution)  # azimuthal angle
      phi = np.linspace(0, np.pi, self.resolution)        # polar angle
      theta_grid, phi_grid = np.meshgrid(theta, phi)
      atomStr = self.atomSymbols[idx]
      atomColor = atomColorHashRGB[atomStr]
      radius = float(0.5* self.atomRadii[idx])
      atomCoord = self.atomList[idx]
      x_outer = radius * np.sin(phi_grid) * np.cos(theta_grid)
      y_outer = radius * np.sin(phi_grid) * np.sin(theta_grid)
      z_outer = radius * np.cos(phi_grid)
      self.Figure.add_trace(go.Surface(
          x=x_outer + atomCoord[0], y=y_outer + atomCoord[1], z=z_outer + atomCoord[2],
          opacity=0.3,colorscale=[[0, atomColor], [1, atomColor]],showscale=False,name=f'Atom_{idx}',
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
          bondColor = "grey"
        elif bondStr == "DOUBLE":
          bondColor = "red"
          w = 15
        elif bondStr == "TRIPLE":
          bondColor = "grey"
          w = 30    
        self.Figure.add_trace(go.Scatter3d(
            x=[atom1[0], atom2[0]],
            y=[atom1[1], atom2[1]],
            z=[atom1[2], atom2[2]],
            mode='lines',
            line=dict(
                color=bondColor,
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
        rgb = (*str_to_rgb(colorStr), 0.8)
        r, g, b, a = rgb
        rgba_str = f'rgba({r}, {g}, {b}, {a})'

        R = np.linalg.norm(vector)
        u_hat = vector / R

        v_hat = np.cross(u_hat, orthogonal)
        v_hat /= np.linalg.norm(v_hat)

        r = np.linspace(0, R, resolution)
        theta = np.linspace(0, np.pi, resolution)
        Rg, Tg = np.meshgrid(r, theta)

        X = origin[0] + Rg * np.cos(Tg) * u_hat[0] + Rg * np.sin(Tg) * v_hat[0]
        Y = origin[1] + Rg * np.cos(Tg) * u_hat[1] + Rg * np.sin(Tg) * v_hat[1]
        Z = origin[2] + Rg * np.cos(Tg) * u_hat[2] + Rg * np.sin(Tg) * v_hat[2]

        newFig.add_trace(go.Surface(
        x=X,
        y=Y,
        z=Z,
        opacity=0.5,
        showscale=False,
        colorscale=[[0, rgba_str], [1, rgba_str]],
        name='Semicircle Surface'))
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
      #elif typeStr == "semiArc":
         
       
        