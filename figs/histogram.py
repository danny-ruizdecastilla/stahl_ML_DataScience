#Histogram maker
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
import sys
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))

PENN_BLUE = "#011F5B"


def histogramTemplate():
    template = go.layout.Template()
    template.layout.font = dict(family="Arial", size=18, color="black")
    template.layout.plot_bgcolor = "white"
    template.layout.width, template.layout.height = 1200, 600

    template.layout.xaxis.showline = True
    template.layout.xaxis.linewidth = 5
    template.layout.xaxis.linecolor = "black"
    template.layout.xaxis.showgrid = False
    template.layout.xaxis.mirror = True

    template.layout.yaxis.tickprefix = "<b>"
    template.layout.yaxis.ticksuffix = "</b>"
    template.layout.yaxis.showline = True
    template.layout.yaxis.linewidth = 5
    template.layout.yaxis.linecolor = "black"
    template.layout.yaxis.ticks = "outside"
    template.layout.yaxis.tickwidth = 4
    template.layout.yaxis.showgrid = False
    template.layout.yaxis.mirror = True
    return template


def makeHistogram(df, colStr, saveDir, saveStr):
    fig = px.histogram(
        df,
        x=colStr,
        nbins=10,                           # 10 bins
        color_discrete_sequence=[PENN_BLUE],
    )
    # white outlines keep adjacent bars distinct
    fig.update_traces(marker_line_color="white", marker_line_width=2)
    fig.update_layout(template=histogramTemplate())

    tbePath = Path(saveDir) / f"{saveStr}.html"
    fig.write_html(str(tbePath))

if __name__ == "__main__":
    dfDir = input("Enter the source of the dataframe: ").strip()
    df = pd.read_csv(dfDir)
    colStr = input("Enter the name of the column from the dataframe you want to make a histogram for: ").strip()
    saveDir = input("Enter the directory where you want to save the histogram: ").strip()
    savePath = Path(saveDir)
    savePath.mkdir(parents=True, exist_ok=True)
    saveStr = input("Enter the name for the saved histogram file (without extension): ").strip()
    makeHistogram(df, colStr, savePath, saveStr)
