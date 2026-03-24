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
from figs.plotSubstrates import safeStringHTML
#Danny Ruiz de Castilla | 11/06/2025
def htmlGeneratorCluster(masterDF, outputDir, pngDir, xCol, yCol, clusterStr):
    masterDF = createPNGDF(masterDF, "SMILES", pngDir)
    base64Col = []
    for img in list(masterDF["pngPath"]):
        base64 = png64(img)
        base64Col.append(base64)
    masterDF["base64"] = base64Col
    masterDF = masterDF.drop(columns=['pngPath'])
    clusterList = [ f"Cluster {clstr}" for clstr in set(masterDF["Cluster"] )]
    htmlCluster = "".join(safeStringHTML(cluster) for cluster in clusterList)
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
            let selectedIndices = new Set();
            let currentClusteredData = [];
            const jsonData = {jsonMAST};
            function printCluster() {{
                const selectedClusterStr = document.getElementById("selectedCluster").value;
                const clusterNumber = parseInt(selectedClusterStr.match(/\d+/)[0], 10);
                currentClusteredData = jsonData.filter(p => p["Cluster"] === clusterNumber)
                selectedIndices.clear();
                showPopup(currentClusteredData);
            }}
            function showPopup(clusteredData){{
                const grid = document.getElementById("structureGrid");
                grid.innerHTML = "";

                clusteredData.forEach((point,i) => {{
                    const container = document.createElement("div");
                    container.style.border = "2px solid transparent";
                    container.style.padding = "5px";
                    container.style.cursor = "pointer";
                    const img = document.createElement("img")
                    img.src = point["base64"];
                    img.style.width = "100px";
                    img.style.height = "100px";
                    img.style.objectFit = "contain";
                    container.onclick = () => toggleSelection(i, container);

                    container.appendChild(img);
                    grid.appendChild(container);

                }});

                document.getElementById("popupOverlay").style.display = "flex"

            }}
            function toggleSelection(idx , element){{
            if (selectedIndices.has(idx)) {{
                selectedIndices.delete(idx);
                element.style.border = "2px solid transparent";
                element.style.backgroundColor = "white";

            }} else {{
                selectedIndices.add(idx);
                element.style.border = "2px solid blue";
                element.style.backgroundColor = "#e6f0ff";

                }}

            }}
            function resetSelection(){{
                selectedIndices.clear();
                const grid = document.getElementById("structureGrid").children;
                for (let el of grid) {{
                    el.style.border = "2px solid transparent";
                    el.style.backgroundColor = "white";

                }}
            }}
            function downloadSelected(clusterStr){{
                console.log("currentClusterData:", currentClusteredData);
                console.log("selectedIndices:", selectedIndices);
                const safeLabel = clusterStr.replace(/\s+/g, "_")
                if (selectedIndices.size===0){{
                    alert("No structures selected");
                }}
                const selectedSMILES = Array.from(selectedIndices).map(i => 
                    currentClusteredData[i]["SMILES"]
                );

                const smilesStr = selectedSMILES.join("\\n");
                const blob = new Blob([smilesStr] , {{type: "text/plain"}});
                const url = URL.createObjectURL(blob);

                const a = document.createElement("a");
                a.href = url;
                a.download = `selected_smiles_${{safeLabel}}.txt`;
                a.click();
                URL.revokeObjectURL(url);

            }}
            function closePopup() {{
                document.getElementById("popupOverlay").style.display = "none";
            }}
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
                        title: {{text: '{xCol}', font: {{family: 'Arial', size: 16, weight: 'bold'}}}},
                        tickfont: {{family: 'Arial', size: 14, weight: 'bold'}},
                        linecolor: 'black',
                        linewidth: 2,
                        mirror: true,
                        showgrid: true,
                        zeroline: false
                    }},
                    yaxis: {{
                        title: {{text: '{yCol}', font: {{family: 'Arial', size: 16, weight: 'bold'}}}},
                        tickfont: {{family: 'Arial', size: 14, weight: 'bold'}},
                        linecolor: 'black',
                        linewidth: 2,
                        mirror: true,
                        showgrid: true,
                        zeroline: false
                    }},
                    hovermode: 'closest',
                    showlegend: false,
                }};

                const config = {{
                    responsive: true,
                    displayModeBar: true,
                    toImageButtonOptions: {{
                        format: 'png',
                        filename: '{plotTitle}_scatterplot',
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

        </script>
        <style>
            body {{
                font-family: Arial, sans-serif;
                margin: 20px;
                background-color: #f5f5f5;
            }}
            .cluster-selection {{
                display: flex;
                align-items: center;
                margin-bottom: 15px;
                flex-wrap: wrap;
                gap: 15px;
                position: absolute;
                top: 600px;
                left: 20px;
            }}
            .generateStructures {{
                display: flex;
                align-items: center;
                flex-wrap: wrap;
            }}
            .popup-overlay {{
                position: fixed;
                top: 0;
                left: 0;
                width: 100%;
                height: 100%;
                background: rgba(0,0,0,0.6);
                display: flex;
                justify-content: center;
                align-items: center;
                z-index: 1000;
            }}

            .popup-content {{
                background: white;
                padding: 20px;
                border-radius: 10px;
                max-width: 80%;
                max-height: 80%;
                overflow-y: auto;
            }}
            .close-btn {{
                float: right;
                cursor: pointer;
                font-size: 24px;
            }}

            #structureGrid {{
                display: grid;
                grid-template-columns: repeat(auto-fill, 120px);
                gap: 10px;
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
    
        <div class="cluster-selection">
            <label for="selectedCluster">Selected Cluster:</label>
            <select id="selectedCluster">
                {htmlCluster}
            </select>
            
            <div class="generateStructures">
                <button onclick="printCluster()">Generate Substrate List</button>
            </div>

        </div>
        <div id="popupOverlay" class="popup-overlay" style="display:none;">
            <div class="popup-content">
                <span class="close-btn" onclick="closePopup()">&times;</span>
                <h3>Cluster Structures</h3>

                <button onclick="downloadSelected(document.getElementById('selectedCluster').value)">Download Selected SMILES</button>
                <button onclick="resetSelection()">Reset Selection</button>

                <div id="structureGrid"></div>
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

    htmlGeneratorCluster(dfNew, outputDir, outputDir / "png" , dim1Str, dim2Str, dimStr + "_kMeans")

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
    outputDir = Path(outputDir)
    dfFinal = kMeansCluster(dfMAST , axis1 , axis2 , kMin, kMax , dimStr , outputDir)
    outputDir.mkdir(parents=True, exist_ok=True)
    dfFinal.to_csv(outputDir  / f"clustered_{dfName}.csv")