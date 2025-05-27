import pandas as pd
import os 
import sys
import glob
import json
import numpy as np
import html
import plotly
import plotly.io as pio
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from figs.chemPlotlyV2 import createPNGDF,png64
from figs.stericvselectroPCA import pcafeatureSplitter
from DFTWorkflow.pitchingATent import compressData , locateNans , eliminateNans , convertCanonical , featureFiltering
from DFTWorkflow.featureMaping import  createCSV

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

    htmlAxis = "".join([safeStringHTML(axis) for axis in masterDF.select_dtypes(include='number').columns])
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

                if (!jsonData || !Array.isArray(jsonData)) {{
                    console.error('jsonData is not a valid array');
                    return;
                }}

                const xValues = jsonData.map(p => p[xCol]);
                const yValues = jsonData.map(p => p[yCol]);
                
                const trace = {{
                    x: xValues,
                    y: yValues,
                    mode: 'markers',
                    type: 'scatter',
                    text: jsonData.map(p => p.SMILES || 'No SMILES'),
                    customdata: jsonData.map((p, i) => ({{
                        id: p.SMILES || `point_${{i}}`, 
                        image: p.base64 || null
                    }})),
                    hovertemplate:
                        '<b>%{{text}}</b><br>' +
                        `${{xCol}}: %{{x}}<br>` +
                        `${{yCol}}: %{{y}}<br>` +
                        '<extra></extra>',

                    marker: {{
                        size: 12,
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
                    hovermode: 'closest',
                    showlegend: false
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

                Plotly.newPlot('plotContainer', [trace], layout, config).then(() => {{
                    const plotElement = document.getElementById('plotContainer');
                    
                    plotElement.on('plotly_hover', function(data) {{
                        const pointData = data.points[0];
                        const customData = pointData.customdata;
                        const imageContainer = document.getElementById('imageContainer');
                        
                        if (customData && customData.image) {{
                            const imageSrc = customData.image.startsWith('data:') ? 
                                customData.image : 
                                `data:image/png;base64,${{customData.image}}`;
                            
                            imageContainer.innerHTML = `
                                <h3>Hover Image</h3>
                                <p><strong>${{pointData.text}}</strong></p>
                                <img id="hoverImage" 
                                    src="${{imageSrc}}" 
                                    alt="Point Image" 
                                    style="max-width: 100%; height: auto; border: 1px solid #ddd;"
                                    onerror="this.style.display='none'; this.nextElementSibling.style.display='block';">
                                <div style="display:none; color: #999; font-style: italic;">Image failed to load</div>
                            `;
                        }} else {{
                            imageContainer.innerHTML = `
                                <h3>Hover Image</h3>
                                <div class="image-placeholder">No image data for this point</div>
                            `;
                        }}
                    }});
                    
                    plotElement.on('plotly_unhover', function(data) {{
                        const imageContainer = document.getElementById('imageContainer');
                        imageContainer.innerHTML = `
                            <h3>Hover Image</h3>
                            <div class="image-placeholder">Hover over a point to see its image</div>
                        `;
                    }});
                }});
            }}

            window.onload = function() {{
                const xAxisSelect = document.getElementById("xAxis");
                const yAxisSelect = document.getElementById("yAxis");
                
                if (xAxisSelect && yAxisSelect) {{
                    xAxisSelect.selectedIndex = 0;
                    yAxisSelect.selectedIndex = 0;
                    plotData();
                    xAxisSelect.addEventListener("change", plotData);
                    yAxisSelect.addEventListener("change", plotData);
                }}
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
            .plot-container {{
                display: flex;
                gap: 20px;
            }}
            #plotContainer {{
                background-color: white;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
                flex: 1;
                min-height: 500px;
            }}
            #imageContainer {{
                width: 300px;
                background-color: white;
                padding: 20px;
                border-radius: 5px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }}
            #hoverImage {{
                max-width: 100%;
                height: auto;
                border: 1px solid #ddd;
                border-radius: 4px;
            }}
            .image-placeholder {{
                color: #666;
                font-style: italic;
                text-align: center;
                padding: 50px 0;
            }}
        </style>
    </head>
    <body>
        <h1>Literature Based Alkene Epoxidation Substrate Space</h1>

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

        <div class="plot-container">
            <div id="plotContainer"></div>
            <div id="imageContainer">
                <h3>Hover Image</h3>
                <div class="image-placeholder">Hover over a point to see its image</div>
            </div>
        </div>
    </body>
    </html>
    """
    with open(outputDir + "/" + chemistry + "interactiveMASTER.html", "w", encoding="utf-8") as f:
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

