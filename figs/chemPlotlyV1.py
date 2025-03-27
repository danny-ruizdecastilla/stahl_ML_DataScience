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
def dat2DF(fileDir, delimiter=None):
    with open(fileDir, 'r') as file:
        headers = file.readline().strip().split(delimiter)
    df = pd.read_csv(fileDir, delimiter=delimiter, skiprows=1, names=headers)

    return df
def createPNGDF(df , smileStr , saveDir  ):
    smileList = list(df[smileStr])
    pngPaths = []
    for i , smile in enumerate(smileList):
        drawer = chemdraw.Drawer(smile, title=smile)
        pngPath = saveDir + "/" + str(i) + ".png" #has to be entire path from root 
        drawer.draw_img(pngPath)
        pngPaths.append(pngPath)

    df["pngPath"] = pngPaths
    return df 
def png64(imagePath):
    with open(imagePath, "rb") as img_file:
        return "data:image/png;base64," + base64.b64encode(img_file.read()).decode()

def interactiveFigGenerator(mainDF , backgroundDF , partition):
    symbols = np.where(mainDF["Yield"] > partition, "circle", "circle-open")

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=mainDF["PC1"], 
        y=mainDF["PC2"], 
        mode="markers", 
        name="Dataset 1",
        marker=dict(symbol=symbols, size=10, color='blue'),
    ))

    fig.add_trace(go.Scatter(
        x=backgroundDF["PC1"], 
        y=backgroundDF["PC2"], 
        mode="markers", 
        name="Dataset 2",
        marker=dict(color='grey' , size = 6), 
        opacity=0.1
    ))


    fig.update_traces(hoverinfo="none", hovertemplate=None)

    fig.update_layout(
        xaxis=dict(title='PC 1' , scaleanchor="y"),
        yaxis=dict(title='PC 2'),  
        plot_bgcolor='rgba(255,255,255,0.1)' , width = 900 , length = 600  , margin=dict(l=60, r=60, t=50, b=60)
    )



    return fig 

def create_dash_app(fig, df , mainStr , background):

    app = Dash(__name__)
    
    app.layout = html.Div([
        dcc.Graph(id=mainStr, figure=fig, clear_on_unhover=True),
        dcc.Tooltip(id="graph-tooltip"),
    ])
    
    @callback(
        Output("graph-tooltip", "show"),
        Output("graph-tooltip", "bbox"),
        Output("graph-tooltip", "children"),
        Input(mainStr, "hoverData")
    )
    def display_hover(hoverData):
        if hoverData is None:
            return False, no_update, no_update
        
        # Get the first hovered point
        pt = hoverData["points"][0]
        bbox = pt["bbox"]
        trace_index = pt["curveNumber"]
        if trace_index == 0:
            pointID = pt["pointNumber"]
            
            df_row = df.iloc[pointID]
            
            img_file = df_row["pngPath"]
            name = df_row['SMILES']

            img_src = png64(img_file)
            
            # Create tooltip children
            children = [
                html.Div([
                    html.Img(src=img_src, style={"width": "100%"}),
                    html.H2(f"{name}", style={"color": "darkblue", "overflow-wrap": "break-word"}),
                ], style={'width': '200px', 'white-space': 'normal'})
            ]
            return True, bbox, children

        else: 
            return False, no_update, no_update
        
    
    return app
if __name__ == "__main__":
    pathDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
        print("created Dir")
    partitionVal = float(sys.argv[3])
    chemistryStr = str(sys.argv[4])
    substrateFiles = glob.glob(pathDir  + "*.dat")
    #print(substrateFiles)
    for file in substrateFiles:
        fileName = file.split("/")[-1]
        if "Grey" in fileName:
            datasetStr1 = "BackgroundSubstrates" #no need to create images for these 

            backgroundDF = dat2DF(file , ",")
        else:
            datasetStr2 = chemistryStr
            chemDF = dat2DF(file , ",")
    chemDF = createPNGDF(chemDF , "SMILES" , outputDir)
    fig = interactiveFigGenerator(chemDF , backgroundDF , partitionVal)
    app = create_dash_app(fig , chemDF  , datasetStr2 , datasetStr1)
    app.run(debug=True)
