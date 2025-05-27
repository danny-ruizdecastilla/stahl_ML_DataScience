import pandas as pd
import os 
import sys
import glob
import json
import chemdraw
import base64
import numpy as np
import html
import plotly
import plotly.io as pio
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from figs.featurePlotter import getFeaturePairs 
from figs.chemPlotlyV1 import insertIntoDataframe
from figs.chemPlotlyV2 import plotly_template , interactiveFigGenerator,createPNGDF,png64
from figs.stericvselectroPCA import pcafeatureSplitter
from DFTWorkflow.pitchingATent import compressData , locateNans , eliminateNans , convertCanonical , featureFiltering
from DFTWorkflow.featureMaping import savePNG , createCSV

def safeStringHTML(unsafeString):
    display_name = html.escape(unsafeString)
    value_name = html.escape(unsafeString, quote=True)
    return f'<option value="{value_name}">{display_name}</option>'
def main(substrateData,chemistry ,  figDir , axisMotifs, eliminatedPhrases, outputDir ):
    initdataSets = glob.glob(substrateData + "/*.csv")
    initdataSets = sorted(initdataSets)
    Xdataframe , smileList  , yieldList_= compressData(initdataSets , "Yield" , eliminatedPhrases)
    nanDict = locateNans(Xdataframe)
    if len(nanDict) != 0:
        Xdataframe["SMILES"] = smileList
        Xdataframe = eliminateNans(Xdataframe , nanDict)
    smileList = Xdataframe["SMILES"].copy()
    canonicalSMILES = []
    for smile in smileList:
        canonical = convertCanonical(smile)
        canonicalSMILES.append(canonical)
    Xdataframe = Xdataframe.drop("SMILES", axis=1)
    featureLabels = list(Xdataframe.columns)
    
    X , featureLabels  = featureFiltering(outputDir, Xdataframe ,featureLabels , chemistry)

    axisDF , axisMotifs = pcafeatureSplitter(X , axisMotifs , 1 , outputDir)
    axisDF["canonicalSMILES"] = canonicalSMILES
    axisDF["SMILES"] = smileList
    createCSV(axisDF , outputDir , "master " + str(chemistry)+ " Dataframe")
    masterDF = createPNGDF(axisDF ,"SMILES" , figDir)
    base64Col = []
    for img in list(masterDF["pngPath"]):

        base64 = png64(img)
        base64Col.append(base64)
    masterDF["base64"] = base64Col
    masterDF = masterDF.drop(columns=['pngPath'])
    print(masterDF.columns.tolist())
    jsonMAST = masterDF.to_json(orient="records" , force_ascii=False )
    jsonMAST = jsonMAST.replace('\\/', '/')
    with open(outputDir + "/" + chemistry + "MAST.json", "w" , encoding='utf-8') as f:
        f.write(jsonMAST)

    htmlAxis = "".join([safeStringHTML(axis) for axis in list(masterDF.columns)])
    columnMaps = {}
    for col in masterDF.columns:
        columnMaps[col] = col

    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Scatter Plot with Hover Images</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <script>
            const jsonData = {jsonMAST};
            const columnMapping = {json.dumps(columnMaps, ensure_ascii=False)};

            function plotData() {{
                const xCol = document.getElementById("xAxis").value;
                const yCol = document.getElementById("yAxis").value;

                const xValues = jsonData.map(p => p[xCol]);
                const yValues = jsonData.map(p => p[yCol]);

                const hoverText = jsonData.map(p =>
                    p.base64
                        ? `<img src="${{p.base64}}" style="width:150px;height:auto;">`
                        : "No image available"
                );

                const trace = {{
                    x: xValues,
                    y: yValues,
                    mode: 'markers',
                    type: 'scatter',
                    customdata: hoverText,
                    hovertemplate: '%{{customdata}}<extra></extra>',
                    marker: {{
                        size: 10,
                        color: 'rgba(55, 128, 191, 0.7)',
                        line: {{
                            color: 'rgba(55, 128, 191, 1.0)',
                            width: 1
                        }}
                    }}
                }};

                const layout = {{
                    title: 'Literature derived Epoxidation Substrate Space',
                    xaxis: {{
                        title: xCol,
                        showgrid: true,
                        zeroline: false
                    }},
                    yaxis: {{
                        title: yCol,
                        showgrid: true,
                        zeroline: false
                    }},
                    hoverlabel: {{
                        align: 'left',
                        bgcolor: 'white',
                        bordercolor: '#ccc',
                        font: {{ color: 'black', size: 12 }},
                        namelength: -1
                    }},
                    hovermode: 'closest'
                }};

                const config = {{
                    responsive: true,
                    displayModeBar: true,
                    toImageButtonOptions: {{
                        format: 'png',
                        filename: 'scatter_plot',
                        height: 600,
                        width: 800,
                        scale: 1
                    }}
                }};

                Plotly.newPlot('plot', [trace], layout, config);
            }}

            window.onload = function() {{
                document.getElementById("xAxis").selectedIndex = 0;
                document.getElementById("yAxis").selectedIndex = 0;
                plotData();
                document.getElementById("xAxis").addEventListener("change", plotData);
                document.getElementById("yAxis").addEventListener("change", plotData);
            }};
        </script>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 20px;
                background-color: #f5f5f5;
            }}
            .controls {{
                background-color: white;
                padding: 20px;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                margin-bottom: 20px;
            }}
            select {{
                padding: 8px;
                margin: 0 10px;
                border: 1px solid #ccc;
                border-radius: 4px;
                font-size: 14px;
            }}
            label {{
                font-weight: bold;
                margin-right: 5px;
            }}
            #plot {{
                background-color: white;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
        </style>
    </head>
    <body>
        <h1>Interactive Scatter Plot with Hover Images</h1>

        <div class="controls">
            <label for="xAxis">X-axis:</label>
            <select id="xAxis">
                {htmlAxis}
            </select>

            <label for="yAxis">Y-axis:</label>
            <select id="yAxis">
                {htmlAxis}
            </select>
        </div>

        <div id="plot" style="width: 100%; height: 600px;"></div>
    </body>
    </html>
    """

    # 6. Write to file
    with open("scatter_plot_with_hover_images.html", "w", encoding="utf-8") as f:
        f.write(html)


if __name__ == "__main__":
    chemistriesDir = str(sys.argv[1])
    figDir = str(sys.argv[2])
    if not os.path.exists(figDir):
        os.makedirs(figDir)
    axisMotifs = {"sterics" : ["distance" , "Buried" , "angle" , "dihedral" , "Vbur"] , "electronics" : ["fuk" , "μ" , "ω" , "Dipole" , "NBO" , "polar" , "HOMO" , "NMR" , "η" ]}
    elimFile = str(sys.argv[3])
    chemistry = str(sys.argv[4])
    outputDir = str(sys.argv[5])
    if not os.path.exists(outputDir): 
        os.makedirs(outputDir)
    if os.path.exists(elimFile):
        with open(elimFile, 'r') as file:
            content = file.read()
            eliminatedPhrases = [item.strip() for item in content.split(',') if item.strip()]
    else: 
        eliminatedPhrases = ["SMILES" , "Compound_Name", "Yield", "ChemistryType"  ]
    
    main( chemistriesDir,chemistry ,  figDir , axisMotifs, eliminatedPhrases , outputDir )

