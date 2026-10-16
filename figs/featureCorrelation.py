#creates the visual matrix of features used to visualize feature spearman and pearson correlation 
import numpy as np
import os 
import sys
import glob
from dash import Dash, html
import plotly.figure_factory as ff
import plotly.graph_objects as go
import pandas as pd
from scipy import stats
from pathlib import Path
def openPlotlyTemplate():#for use when plotting correlation matrix's , variable typesettings based on length of features \
    template = go.layout.Template()
    template.layout.font = dict(family="Arial", size=8, color="black")
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
    mask = np.triu(np.ones(corrDF.shape, dtype=bool), k=1)
    z_masked = corrDF.mask(mask)
    if template is None:
        plotTemplate = openPlotlyTemplate
    if savePath is None:
        save_path = False
    else:
        save_path = True

    fig = go.Figure(
    data=go.Heatmap(
        z=z_masked.to_numpy(),
        x=corrDF.columns,
        y=corrDF.columns,
        colorscale="RdBu",
        colorbar=dict(title=corrStr)
    ),
        layout=dict(template=plotTemplate())
    )
    annotations = []
    for i, row in enumerate(corrDF.values):
        for j, val in enumerate(row):
            if i >= j:
                font_color = "black" if abs(val) < 0.6 else "white"
                text_val = str(round(val, 2))
            else:
                text_val = ""
                font_color = "rgba(0,0,0,0)"  

            annotations.append(
                dict(
                    x=corrDF.columns[j],
                    y=corrDF.index[i],
                    text=text_val,
                    showarrow=False,
                    font=dict(color=font_color)
                )
            )
    
    fig.update_layout(
        title=f"{corrStr} Correlation Matrix",
        xaxis=dict(
            tickangle=45,
            constrain="domain"
        ),
        yaxis=dict(
            autorange="reversed"
        ), 
        annotations=annotations
    )
    if save_path: 
        fig.write_html(savePath) 
        print(f"✅ Correlation matrix saved to {savePath}")
    
    return fig
def main(df , corrType , savePath):
    dfMain = pd.read_csv(df , encoding="utf-8") 
    featureLabels = dfMain.columns.tolist()
    if corrType == "Spearman":
        result = stats.spearmanr(dfMain.values, axis=0)
        correlationMatrix = result.correlation
        spearmanDF = pd.DataFrame(correlationMatrix, index=featureLabels, columns=featureLabels)
        correlationGenerator(spearmanDF, "Spearman", Path(savePath) / "spearman.html", None)
    elif corrType == "Pearson":
        correlationMatrix = np.corrcoef(dfMain.values, rowvar=False)
        pearsonDF = pd.DataFrame(correlationMatrix, index=featureLabels, columns=featureLabels)
        correlationGenerator(pearsonDF, "Pearson", Path(savePath) / "pearson.html", None )
    else:
        raise ValueError("Invalid correlation type. Choose 'pearson' or 'spearman'.")

if __name__ == "__main__":
    corrDF = sys.argv[1]
    corrStr = sys.argv[2]
    savePath = sys.argv[3]
    if not os.path.exists(os.path.dirname(savePath)):
        os.makedirs(os.path.dirname(savePath))
    main(corrDF, corrStr , savePath)