from dash import Dash, html, dcc, callback, Output, Input, callback_context
import plotly.graph_objects as go
import plotly 
import pandas as pd
import os 
import sys
import glob
import chemdraw
import base64
import numpy as np
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from DFTWorkflow.featureMaping import createCSV
def dat2DF(fileDir, delimiter=None):
    with open(fileDir, 'r') as file:
        headers = file.readline().strip().split(delimiter)
    df = pd.read_csv(fileDir, delimiter=delimiter, skiprows=1, names=headers)

    return df
def createPNGDF(df, smileStr, saveDir):
    if not os.path.exists(saveDir):
        os.makedirs(saveDir)
    smileList = list(df[smileStr])
    pngPaths = []
    
    for idx, smile in zip(df.index, smileList):
        drawer = chemdraw.Drawer(smile, title=smile)
        pngPath = f"{saveDir}/{idx}.png" 
        drawer.draw_img(pngPath)
        pngPaths.append(pngPath)

    df["pngPath"] = pngPaths
    return df
def png64(imagePath):
    with open(imagePath, "rb") as img_file:
        return "data:image/png;base64," + base64.b64encode(img_file.read()).decode()

def interactiveFigGenerator(mainDF , backgroundDF , partition , dataStr1 , dataStr2):
    symbols = np.where(mainDF["Yield"] >= partition, "circle", "circle-open")
    colors = np.where(mainDF["Yield"] >= partition, "blue", "red")
    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=mainDF["PC1"], 
        y=mainDF["PC2"], 
        mode="markers", 
        name= str(dataStr2),
        marker=dict(symbol=symbols, size=10, color=colors),
    ))

    fig.add_trace(go.Scatter(
        x=backgroundDF["PC1"], 
        y=backgroundDF["PC2"], 
        mode="markers", 
        name= str(dataStr1) ,
        marker=dict(color='grey' , size = 8), 
        opacity=0.4
    ))


    fig.update_traces(hoverinfo="none", hovertemplate=None)

    fig.update_layout(
        xaxis=dict(title='PC 1', scaleanchor="y"),  # Keeps x and y scales equal
        yaxis=dict(title='PC 2'),
        plot_bgcolor='rgba(255,255,255,0.1)',  # Light background transparency
        width=600,  
        height=600,  
        margin=dict(l=60, r=60, t=50, b=60),  
        legend=dict(
            x=0.6,  
            y=0.98,  
            bgcolor="rgba(255,255,255,0.7)",  
            bordercolor="black",
            borderwidth=1
        ),
        title=dict(
            text=f"PCA for {dataStr2} Alkenes at {partition}% Yield",  # Fixed f-string
            font=dict(size=18, color="black"),  
            x=0.5,  # Center the title
            y=0.95  
        )
    )
    return fig 

def create_dash_app(figDict, dfDict, strList):
    app = Dash(__name__)
    
    app.layout = html.Div([
        html.Div([
            dcc.Graph(id=fig_id, figure=fig, clear_on_unhover=True)
            for fig_id, fig in figDict.items()
        ]),
        dcc.Tooltip(id="graph-tooltip"),
    ])
    
    @callback(
        Output("graph-tooltip", "show"),
        Output("graph-tooltip", "bbox"),
        Output("graph-tooltip", "children"),
        *[Input(fig_id, "hoverData") for fig_id in figDict.keys()]
    )
    def display_tooltip(*hover_data):
        ctx = callback_context
        if not ctx.triggered:
            return False, {}, []
            
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
        
        for i, fig_id in enumerate(figDict.keys()):
            if fig_id == triggered_id and hover_data[i] is not None:
                hover = hover_data[i]
                df_id = list(dfDict.keys())[i]
                
                pt = hover["points"][0]
                bbox = pt["bbox"]
                trace_index = pt["curveNumber"]
                
                # Make sure we're only processing points from the first trace (blue circles)
                if trace_index == 0:
                    # Get the x,y coordinates of the point
                    x_val = pt["x"]
                    y_val = pt["y"]
                    
                    # Find the corresponding row in the dataframe by matching coordinates
                    if df_id in dfDict:
                        # Find the row with matching coordinates
                        mask = (dfDict[df_id]["PC1"] == x_val) & (dfDict[df_id]["PC2"] == y_val)
                        if any(mask):
                            df_row = dfDict[df_id][mask].iloc[0]  # Take the first match
                            img_file = df_row["pngPath"]
                            
                            # Verify the image file exists
                            if os.path.exists(img_file):
                                name ="Yield:" + str(df_row['Yield'])
                                try:
                                    img_src = png64(img_file)
                                    children = [
                                        html.Div([
                                            html.Img(src=img_src, style={"width": "100%"}),
                                            html.H2(f"{name}", style={"color": "darkblue", "overflow-wrap": "break-word"}),
                                        ], style={'width': '200px', 'white-space': 'normal'})
                                    ]
                                    return True, bbox, children
                                except Exception as e:
                                    # Display error if image can't be loaded
                                    children = [html.Div(f"Error loading image: {str(e)}")]
                                    return True, bbox, children
                            else:
                                # Image file doesn't exist
                                children = [html.Div(f"Image not found: {img_file}")]
                                return True, bbox, children
        
        return False, {}, []
    
    return app
def getPathDirs(mainDir , string1 , string2 ):
    fileDict = {}
    allDirs = glob.glob(mainDir + "/*/")
    print(allDirs)
    for dir_ in allDirs:
        if len(glob.glob(dir_ + "*" +str(string1) +  ".dat")) == 1 and len(glob.glob(dir_ + "*" +str(string2) +  ".dat")) == 1:
            substrates = glob.glob(dir_ + "*" +str(string1) +  ".dat")[0]
            greyedOut = glob.glob(dir_ + "*" +str(string2) +  ".dat")[0]
            keyString = str(dir_.split("/")[-2])
            fileDict[keyString] = [ substrates ,   greyedOut ]
        elif len(glob.glob(dir_ + "*" +str(string1) +  ".dat")) == 0 or len(glob.glob(dir_ + "*" +str(string2) +  ".dat")) == 0:
            print(f"Directory {dir_} is missing an input file")
            
        else:
            print(f"Error: Too many Input Files for Plotly map")
            quit()
    return fileDict
if __name__ == "__main__":
    pcaDir = str(sys.argv[1])
    backgroundStr = str(sys.argv[2])
    substrateStr = str(sys.argv[3])
    pathDict = getPathDirs(pcaDir , backgroundStr , substrateStr)
    partitionVal = float(sys.argv[4])
    chemistryStr = str(sys.argv[5]) #Alkene
    chemDict = {}
    dfDict = {}
    saveStrList = []
    for i , key in enumerate(list(pathDict.keys())):
        print("i" , i)
        saveDir = pcaDir + "/" + key 
        if not os.path.exists(saveDir):
            os.makedirs(saveDir)
        backFile = pathDict[key][0]
        substrateFile = pathDict[key][1]
        backgroundDF = dat2DF(backFile , ",")

        chemDF = dat2DF(substrateFile , ",")
        chemDF = createPNGDF(chemDF , "SMILES" , saveDir + "/png")
        createCSV(chemDF , saveDir + "/", "chemDF")
        fig = interactiveFigGenerator(chemDF , backgroundDF , partitionVal , "Background Substrates", str(key) + " substrates")
        chemDict[key] = fig
        dfDict[key] = chemDF
        saveStrList.append([key , "Background Substrates"])
    app = create_dash_app(chemDict , dfDict , saveStrList )
    app.run(debug=True)
