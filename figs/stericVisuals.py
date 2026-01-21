from dash import Dash, dcc, html, Input, Output, no_update, callback
import plotly.graph_objects as go
import plotly 
import pandas as pd
import os 
import sys
import glob
import chemdraw
import base64
import numpy as np
import json
from pathlib import Path
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
def string_to_hex_color(s):
  """Convert any string to a consistent hex color."""
  hash_value = hash(s)
  color_int = abs(hash_value) % 0xFFFFFF
  return f'#{color_int:06x}'
class stericDrawer:
    def __init__(self, atomsHash):
        self.atomList = atomsHash["atomCoords"]
        self.atomSymbols = atomsHash["atomSymbol"]
        self.atomRadii = atomsHash["radii"]
        fig = go.Figure()
        self.Figure = fig
    def drawAtom(self ,idx,  resolution):
      theta = np.linspace(0, 2 * np.pi, resolution)  # azimuthal angle
      phi = np.linspace(0, np.pi, resolution)        # polar angle
      theta_grid, phi_grid = np.meshgrid(theta, phi)
      atomStr = self.atomSymbols[idx]
      atomColor = string_to_hex_color(atomStr)
      radius = float(0.01 * self.atomRadii[idx])
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
        if bondStr == "single":
          w = 5
        elif bondStr == "double":
           w = 8
        elif bondStr == "triple":
           w = 12
        bondColor = string_to_hex_color(bondStr)
        
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
    def drawShapes(self , shapeHash ):
      typeStr = shapeHash["shapeType"]
      if typeStr == "line":
        origin = shapeHash["origin"]
        vector = shapeHash["vector"]
        self.Figure.add_trace(go.Scatter3d(
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
      elif typeStr == "semiCylinder":
        radius = shapeHash["radius"]
        origin = shapeHash["origin"]
        resolution = shapeHash["resolution"]
        height = shapeHash["length"]
        theta = np.linspace(0, np.pi, resolution)
        z = np.linspace(0, height, resolution)
        theta_grid, z_grid = np.meshgrid(theta, z)
        x_outer = radius * np.cos(theta_grid)
        y_outer = radius * np.sin(theta_grid)
        z_outer = z_grid
        self.Figure.add_trace(go.Surface(
        x=x_outer + origin[0], y=y_outer + origin[1], z=z_outer + origin[2],
        opacity=0.3,
        colorscale=[[0, 'steelblue'], [1, 'steelblue']],
        showscale=False,
        name='Outer Wall',
        hoverinfo='skip'))

