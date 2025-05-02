from dash import Dash, html, dcc, callback, Output, Input, callback_context
import plotly.graph_objects as go
import plotly 
import pandas as pd
import os 
import sys
import glob
import re
import chemdraw
import base64
import numpy as np
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from figs.chemPlotlyV2 import createPNGDF , png64 , plotly_template
from DFTWorkflow.featureMaping import savePNG , createCSV
from DFTWorkflow.pitchingATent import compressData , locateNans , eliminateNans , convertCanonical
#Danny Ruiz de Castilla 04.27.
def plotScatter( backgroundDF , xAxisStr , yAxisStr):
    fig = go.Figure(layout=dict(template=plotly_template()))
    fig.add_trace(go.Scatter(
        x=backgroundDF[xAxisStr], 
        y=backgroundDF[yAxisStr], 
        mode="markers", 
        name= str("All Substrates") ,
        marker=dict(color='grey' , size = 10), 
        opacity=0.4
    ))
    fig.update_traces(hoverinfo="none", hovertemplate=None)

    fig.update_layout(
    xaxis=dict(title=xAxisStr, scaleanchor="y"),  # Keeps x and y scales equal
    yaxis=dict(title=yAxisStr),
    plot_bgcolor='rgba(255,255,255,0.1)',  # Light background transparency
    title=dict(
        text=(str(yAxisStr) + " over " + str(xAxisStr)),  # Fixed concatenation
        font=dict(size=18, color="black"),  
        x=0.5,  # Center the title
        y=0.95  
    )
)
    return fig 
def dashScatter(figDict, pngDF):
    app = Dash(__name__)
    
    app.layout = html.Div([
        html.Div([dcc.Graph(id=fig_id, figure=fig, clear_on_unhover=True)for fig_id, fig in figDict.items()]),
        dcc.Tooltip(id="graph-tooltip"),])
    
    @callback(Output("graph-tooltip", "show"),Output("graph-tooltip", "bbox"),Output("graph-tooltip", "children"),*[Input(fig_id, "hoverData") for fig_id in figDict.keys()])
    def display_tooltip(*hover_data):
        ctx = callback_context
        if not ctx.triggered:
            return False, {}, []
            
        triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
        
        for i, fig_id in enumerate(figDict.keys()):
            if fig_id == triggered_id and hover_data[i] is not None:
                hover = hover_data[i]
                pt = hover["points"][0]
                bbox = pt["bbox"]
                trace_index = pt["curveNumber"]
                xAxis = str(fig_id.split("_")[0])
                yAxis = str(fig_id.split("_")[1])
                if trace_index == 0:
                    # Get the x,y coordinates of the point
                    x_val = pt["x"]
                    y_val = pt["y"]

                    mask = (pngDF[xAxis] == x_val) & (pngDF[yAxis] == y_val)
                    if any(mask):
                        df_row = pngDF[mask].iloc[0]  # Take the first match
                        img_file = df_row["pngPath"]

                        if os.path.exists(img_file):
                            try:
                                img_src = png64(img_file)
                                children = [
                                    html.Div([
                                        html.Img(src=img_src, style={"width": "100%"}),
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
def getFeaturePairs(featureList , motifStrs): 
    xAxisList = []
    yAxisList = []
    for feature in featureList:
        if any(substr in feature for substr in motifStrs[0]):
            xAxisList.append(feature)
        elif any(substr in feature for substr in motifStrs[1]):
            yAxisList.append(feature)
    
    masterFeatureList = [[x, y] for x in xAxisList for y in yAxisList]
    return masterFeatureList 
def removeNonLetters(text):
    return re.sub(r'[^\w\u0370-\u03FF\u1F00-\u1FFF]', '', text)
def standardCols(df):
    df.columns = [removeNonLetters(str(col)) for col in df.columns]
    return df
def main(substrateData , outputDir , elimFile , motifList ):
    if not os.path.exists(outputDir): 
        os.makedirs(outputDir)
    if os.path.exists(elimFile):
        with open(elimFile, 'r') as file:
            content = file.read()
            eliminatedPhrases = [item.strip() for item in content.split(',') if item.strip()]
    else: 
        eliminatedPhrases = ["SMILES" , "Compound_Name", "Yield", "ChemistryType"  ]
    initdataSets = glob.glob(substrateData + "/*.csv")
    initdataSets = sorted(initdataSets)
    Xdataframe , smileList  , yieldList_= compressData(initdataSets , "Yield" , eliminatedPhrases)
    standardXdf = standardCols(Xdataframe)
    Xdataframe["SMILES"] = smileList
    nanDict = locateNans(standardXdf)
    if len(nanDict) != 0:
        standardXdf = eliminateNans(standardXdf, nanDict)
    smileList = standardXdf["SMILES"].copy()
    canonicalSMILES = []
    for smile in smileList:
        canonical = convertCanonical(smile)
        canonicalSMILES.append(canonical)

    featureDF = standardXdf.drop("SMILES", axis=1)
    
    featureLabels = list(featureDF.columns)

    featurePairs = getFeaturePairs(featureLabels , motifList)
    if not os.path.exists(outputDir + "/figs/" ): 
        os.makedirs(outputDir + "/figs/" )
    pngDF = createPNGDF(standardXdf , "SMILES" , outputDir + "/figs/"  )
    figDict = {}
    print(featurePairs)
    for pair in featurePairs:
        
        xStr = str(pair[0])
        yStr = str(pair[1])
        fig = plotScatter(pngDF , xStr , yStr)

        figDict[xStr + "_" + yStr] = fig
    
    app = dashScatter(figDict, pngDF)
    app.run(debug=True)
if __name__ == "__main__":
    datasetDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    elimFile = str(sys.argv[3])

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    axisMotifs = [ ["distance" , "Buried" , "angle" , "dihedral" , "%Vbur" , "PCA1"], ["μ" , "ω" , "Dipole" , "NBO" , "polar" , "HOMO" , "NMR" , "η" , "PCA2"]] #X axis is steric features, y axis is electronic features , each has its specific strings 
    
    main(datasetDir , outputDir , elimFile , axisMotifs)
