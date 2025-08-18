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
def correlationGenerator(matrx , corrStr , plotTemplate):
    isSymmetric = diagonalChecker(matrx)
    if isSymmetric:
        fig = go.Figure(layout=dict(template=plotTemplate()))
        

    else:
        print("Error in correlation Generator. Input correlation matrix is not symmetric")