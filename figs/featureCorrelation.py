#creates the visual matrix of features used to visualize feature spearman and pearson correlation 
import numpy as np
import os 
import sys
import glob
from dash import Dash, html
import plotly.figure_factory as ff
import plotly.graph_objects as go
import pandas as pd

def openPlotlyTemplate():#for use when plotting correlation matrix's , variable typesettings based on length of features \
    template = go.layout.Template()
    template.layout.font = dict(family="Arial", size=18, color="black")
    template.layout.plot_bgcolor = "white"
    template.layout.xaxis.linewidth = 5
    template.layout.xaxis.linecolor = "black"
    template.layout.xaxis.showgrid = False
    template.layout.xaxis.tickangle = 45 
    template.layout.yaxis.linewidth = 5
    template.layout.yaxis.linecolor = "black"
    template.layout.yaxis.showgrid = False
    return template
def diagonalChecker(matrix):
    isDiagonal =True
    for i in range(len(matrix)):        
        for j in range(len(matrix[0])):  
            val1 = matrix[i][j]
            val2 = matrix[j][i]
            if val1 != val2:
                isDiagonal = False
    return isDiagonal 
def correlationGenerator(corrDF, corrStr , savePath , template: str = None):
    if template is None:
        plotTemplate = openPlotlyTemplate
    if savePath is None:
        save_path = False
    else:
        save_path = True

    isSymmetric = diagonalChecker(corrDF.values)
    if not isSymmetric:
        print("Error: Input correlation matrix is not symmetric")
        return None
    fig = go.Figure(
        data=go.Heatmap(
            z=corrDF.values,
            x=corrDF.columns,
            y=corrDF.columns,
            colorscale="RdBu", 
            colorbar=dict(title=corrStr, titleside="right")
        ),
        layout=dict(template=plotTemplate())
    )
    annotations = []
    for i, row in enumerate(corrDF.values):
        for j, val in enumerate(row):
            annotations.append(
                dict(
                    x=corrDF.columns[j],
                    y=corrDF.index[i],
                    text=str(round(val, 2)),
                    showarrow=False,
                    font=dict(color="black" if abs(val) < 0.6 else "white") # contrast for readability
                )
            )
    
    fig.update_layout(
        title=f"{corrStr} Correlation Matrix",
        xaxis=dict(tickangle=45),
        yaxis=dict(autorange="reversed"),
        annotations=annotations
    )
    if save_path:
        fig.write_html(savePath)
        print(f"✅ Correlation matrix saved to {savePath}")
    
    return fig