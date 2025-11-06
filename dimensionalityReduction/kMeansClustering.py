import pandas as pd
from sklearn.cluster import KMeans
import plotly
import os
import plotly.graph_objects as go
import sys
import numpy as np
import plotly.express as px
import json
from pathlib import Path
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from dimensionalityReduction.reactivityFeatures import boxGen
from reaxysProcessing.reaxysKMeans_Downselect import needleAlg
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from figs.chemPlotlyV1 import plotly_template
from figs.chemPlotlyV2 import createPNGDF,png64
#Danny Ruiz de Castilla | 11/06/2025
def htmlGeneratorCluster(masterDF, outputDir, pngDir, xCol, yCol, clusterStr):
    masterDF = createPNGDF(masterDF, "SMILES", pngDir)
    base64Col = []
    for img in list(masterDF["pngPath"]):
        base64 = png64(img)
        base64Col.append(base64)
    masterDF["base64"] = base64Col
    masterDF = masterDF.drop(columns=['pngPath'])
    
    jsonMAST = masterDF.to_json(orient="records", force_ascii=False)
    jsonMAST = jsonMAST.replace('\\/', '/')
    plotTitle = input(f"Name your clustered scatterplot: ")

    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Scatter Plot with Hover Images</title>
        <script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script>
        <script>
            const jsonData = {jsonMAST};
            function plotData() {{

                if (!jsonData || !Array.isArray(jsonData)) {{
                    console.error('jsonData is not a valid array');
                    return;
                }}

                const xValues = jsonData.map(p => p['{xCol}']);
                const yValues = jsonData.map(p => p['{yCol}']);
                const clusterValues = jsonData.map(p => p.Cluster);

                const trace = {{
                    x: xValues,
                    y: yValues,
                    mode: 'markers',
                    type: 'scatter',
                    text: jsonData.map(p => `Cluster ${{p.Cluster}}`),
                    customdata: jsonData.map(p => ({{
                        label: `Cluster ${{p.Cluster}}`,
                        image: p.base64,
                        cluster: p.Cluster
                    }})),
                    hovertemplate:
                        '<b>%{{text}}</b><br>' +
                        '{xCol}: %{{x}}<br>' +
                        '{yCol}: %{{y}}<br>' +
                        '<extra></extra>',
                    marker: {{
                        size: 12,
                        color: clusterValues,
                        colorscale: 'Viridis',
                        showscale: true,
                        line: {{
                            color: 'rgba(0,0,0,0.6)',
                            width: 0.5
                        }}
                    }}
                }};

                const layout = {{
                    title: '{plotTitle}',
                    xaxis: {{
                        title: '{xCol}',
                        showgrid: true,
                        zeroline: false
                    }},
                    yaxis: {{
                        title: '{yCol}',
                        showgrid: true,
                        zeroline: false
                    }},
                    hovermode: 'closest',
                    showlegend: false
                }};

                const config = {{
                    responsive: true,
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

        </script>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 20px;
                background-color: #f5f5f5;
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
            .image-placeholder {{
                color: #666;
                font-style: italic;
                text-align: center;
                padding: 50px 0;
            }}
        </style>
    
    </head>
    <body onload="plotData()">
        <h1>{plotTitle}</h1>

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
    outputDir = Path(outputDir)
    outputDir.mkdir(parents=True, exist_ok=True)
    with open(outputDir / (clusterStr + "interactiveMASTER.html"), "w", encoding="utf-8") as f:
        f.write(html)
def kMeansCluster(dfMAST , dim1Str , dim2Str , kMin, kMax , dimStr , outputDir):
    inertiaList = []
    clusterRange = np.arange(kMin, kMax, 1)
    reducedX = dfMAST[[dim1Str,dim2Str]].values
    for k in clusterRange:
        initInertia = []
        count = 0
        while count < 4:
            kmeans = KMeans(n_clusters=k, random_state=count)
            kmeans.fit(reducedX)
            initInertia.append(kmeans.inertia_)
            count +=1 
        inertia = np.mean(initInertia)
        inertiaList.append(inertia)
    specialK = needleAlg(list(clusterRange) , inertiaList)
    kmeansMAST = KMeans(n_clusters=specialK, random_state=42)
    
    elbowFig = go.Figure(layout=dict(template=plotly_template()))

    elbowFig.add_trace(go.Scatter(
        x=clusterRange, y=inertiaList, mode='lines+markers',
        marker=dict(size=8), line=dict(width=2),
        name='Inertia'
    ))

    elbowFig.update_layout(
        title="KMeans Elbow Method",
        xaxis_title="Number of Clusters (k)",
        yaxis_title="Inertia (Within-Cluster Sum of Squares)",
        width=600,
        height=400
    )


    dfMAST['Cluster'] = kmeansMAST.fit_predict(reducedX)

    colsKeep = [dim1Str , dim2Str , "Cluster" , "SMILES"]
    dfNew = dfMAST[colsKeep]

    htmlGeneratorCluster(dfNew, outputDir, outputDir + "/png/", dim1Str, dim2Str, dimStr + "_kMeans")

    outputDir = Path(outputDir)
    outputDir.mkdir(parents=True, exist_ok=True)

    # --- Save elbow plot ---
    elbowPath = outputDir / f"{dimStr}_elbow.html"
    elbowFig.write_html(str(elbowPath))
    print(f"✅ Elbow plot saved to: {elbowPath}")

    return dfMAST

if __name__ == "__main__":
    masterDF = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    kMin = int(sys.argv[3])
    kMax = int(sys.argv[4])
    dimStr = str(sys.argv[5])
    dfName = str(Path(masterDF).name)
    dfMAST = pd.read_csv(masterDF)
    dfCols = list(dfMAST.columns)
    colsBox = boxGen(dfCols)
    clustAxis = listInputs(f"Here are all available columns for the dataframe {colsBox}\n please enter the matching indexes corresponding to the cluster axis's:\n")
    colList = []
    for idx in clustAxis:
        column = dfCols[int(idx)]
        colList.append(column)
    axis1 = colList[0]
    axis2 = colList[1]
    dfFinal = kMeansCluster(dfMAST , axis1 , axis2 , kMin, kMax , dimStr , outputDir)
    outputDir = Path(outputDir)
    outputDir.mkdir(parents=True, exist_ok=True)
    dfFinal.to_csv(outputDir  / f"clustered_{dfName}.csv")