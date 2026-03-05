#violin plot generator 
#Daniel Ruiz de Castilla | 03.04.2026
import plotly.express as px
import plotly.graph_objects as go
from dash import html
from pathlib import Path
import pandas as pd
import sys
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
def violinTemplate(): 
    template = go.layout.Template()
    template.layout.font = dict(family="Arial", size=18, color="black")
    template.layout.plot_bgcolor = "white"
    template.layout.width, template.layout.height = 1200, 600
    template.layout.xaxis.showgrid = False
    template.layout.xaxis.mirror = True
    template.layout.yaxis.tickprefix = "<b>"
    template.layout.yaxis.ticksuffix = "<b>"
    template.layout.yaxis.showline = True
    template.layout.yaxis.linewidth = 5
    template.layout.yaxis.linecolor = "black"
    template.layout.yaxis.ticks = "outside"
    template.layout.yaxis.tickwidth = 4
    template.layout.yaxis.showgrid = False
    template.layout.yaxis.mirror = True
    return template
def askingForAddition(askStr):
    response = input(f"{askStr}").strip().lower()

    if response in {"yes", "y"}:
        return True
    elif response in {"no", "n"}:
        return False
    else:
        print("Please enter a valid 'yes' or 'no' response.")
        return askingForAddition(askStr)   # recursive call
def checkCol(df):
    colStr = input("Enter the name of the column from the dataframe you want to add: ").strip()

    if colStr in df.columns:
        return colStr
    else:
        print(f"'{colStr}' is not a valid column.")
        print("Available columns:")
        print(", ".join(df.columns))
        return checkCol(df)

def compareSources(fig , saveDir , saveStr):
    addSources = askingForAddition("Do you want to add another source to your violin plot? (yes/no): ")

    if addSources:
        dfDir = input("Enter the source of the dataframe: ").strip()
        df = pd.read_csv(dfDir)
        col = checkCol(df)
        fig.add_trace(go.Violin(
            y=df[col].dropna(),
            name=col,
            box_visible=True,
            meanline_visible=True
        ))
        
        fig.update_layout(template=violinTemplate())
        
        return compareSources(fig ,saveDir , saveStr)
    else:
        fig.write_html(saveDir + f"/{saveStr}_Violin.html") 
if __name__ == "__main__":
    saveDir = str(sys.argv[1])
    saveStr = str(sys.argv[2])
    savePath = Path(saveDir)
    savePath.mkdir()
    fig = go.Figure()
    compareSources(fig , saveDir , saveStr)