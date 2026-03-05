#violin plot generator 
#Daniel Ruiz de Castilla | 03.04.2026
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path
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

def make_violin_plots(df, columns, box=True, meanline=True):
    """
    Create violin plots from selected DataFrame columns.
    
    Parameters:
        df (pd.DataFrame): source dataframe
        columns (list): list of column names to plot
        box (bool): show mini boxplot inside violin
        meanline (bool): show mean line
    """
    
    fig = go.Figure()
    
    for col in columns:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in DataFrame")
        
        fig.add_trace(go.Violin(
            y=df[col].dropna(),
            name=col,
            box_visible=box,
            meanline_visible=meanline
        ))
    
    fig.update_layout(template=violinTemplate())
    fig.show()
    
    return fig