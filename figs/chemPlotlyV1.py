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
from rdkit import Chem
def convertCanonical(str):
    mol = Chem.MolFromSmiles(str)
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    return canonical
def dat2List(fileDir, delimiter=None):
    with open(fileDir, 'r') as file:
        entries = list(file.readline().strip().split(delimiter))
    return entries
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
def plotly_template(): #Credit to Dylan Walsh
    template = go.layout.Template()
    template.layout.font = dict(family="Arial", size=18, color="black")
    template.layout.plot_bgcolor = "white"
    template.layout.width, template.layout.height = 1200, 600
    template.layout.xaxis.tickprefix = "<b>"
    template.layout.xaxis.ticksuffix = "<b>"
    template.layout.xaxis.showline = True
    template.layout.xaxis.linewidth = 5
    template.layout.xaxis.linecolor = "black"
    template.layout.xaxis.ticks = "outside"
    template.layout.xaxis.tickwidth = 4
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

def interactiveFigGenerator(mainDF , backgroundDF , partition , dataStr1 , dataStr2):
    symbols = np.where(mainDF["Yield"] >= partition, "circle", "circle-open")
    colors = np.where(mainDF["Yield"] >= partition, "blue", "red")
    fig = go.Figure(layout=dict(template=plotly_template()))
    fig.add_trace(go.Scatter(
        x=mainDF["PCA1"], 
        y=mainDF["PCA2"], 
        mode="markers", 
        name= str(dataStr2),
        marker=dict(symbol=symbols, size=10, color=colors),
    ))

    fig.add_trace(go.Scatter(
        x=backgroundDF["PCA1"], 
        y=backgroundDF["PCA2"], 
        mode="markers", 
        name= str(dataStr1) ,
        marker=dict(color='grey' , size = 6), 
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
            x=0.0,  
            y=0.0,  
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
            name ="Yield: " + str(df_row['Yield'])

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
def insertIntoDataframe(df1 ,  connectingStr ,connectingStr2 ,  df2 , inputStr , indexInterest  ):
    #First create an empty column of zeros at the index of interest and call it inputStr 
    df1.insert(indexInterest, inputStr, 0)
    
    col1 = list(df1[connectingStr])

    col2 = list(df2[connectingStr2])
    insertCol = list(df2[inputStr])
    

    for i , str in enumerate (col1):
        if str in col2:
            index = col2.index(str)
            insertVal = insertCol[index]
            #df1[inputStr][i] = insertVal
            df1.loc[i, inputStr] = insertVal
    return df1
def main(smilesFile ,pcaFile , outputDir , partitionVal , chemistryStr): #csv , csv 
    inputDF = pd.read_csv(smilesFile)
    smilesList = inputDF['SMILES']
    canonicals = []
    for smiles in smilesList:
        canonical = convertCanonical(smiles)
        canonicals.append(canonical)
    inputDF["Canonicals"] = canonicals
    substrateDF = pd.read_csv(pcaFile)
    substrateSmiles = substrateDF['SMILES']
    canonicalMain = []
    for substrate in substrateSmiles:
        canons = convertCanonical(substrate)
        canonicalMain.append(canons)
    substrateDF["Canonicals"] = canonicalMain
    inputDF = insertIntoDataframe(inputDF , "Canonicals" , "Canonicals" , substrateDF , "PCA1" , 0 )
    inputDF = insertIntoDataframe(inputDF , "Canonicals" , "Canonicals" , substrateDF , "PCA2" , 1 )


    greyDF = substrateDF[~substrateDF['Canonicals'].isin(canonicals)].reset_index(drop=True)
    pngDF = createPNGDF(inputDF , "SMILES" , outputDir)

    fig = interactiveFigGenerator(pngDF , greyDF , partitionVal ,  "Background Substrates" , chemistryStr )
    app = create_dash_app(fig , pngDF  ,  "Background Substrates" , chemistryStr)
    app.run(debug=True)


if __name__ == "__main__":
    #Takes in a .csv dataframe of substrates and projects matches onto a given PCA by greying everything else out that doesnt match 
    smilesFile = str(sys.argv[1]) #csv 
    pcaFile = str(sys.argv[2])
    outputDir = str(sys.argv[3])
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
        print("created Dir")
    partitionVal = float(sys.argv[4])
    chemistryStr = str(sys.argv[5])
    #print(substrateFiles)
    main(smilesFile , pcaFile , outputDir , partitionVal , chemistryStr)

