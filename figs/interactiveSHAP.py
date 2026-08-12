#Danny Ruiz de Castilla 07.30.2026
#Prints out an interactive SHAP plot that can use hover to display chemical structures on a seperate screen 
import pandas as pd
import numpy as np
import re
from pathlib import Path
import shap
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import StandardScaler
from sklearn.neural_network import MLPRegressor
import json
import shutil
import html 
import sys 
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from figs.chemPlotlyV2 import createPNGDF,png64
def generateInteractiveSHAP(shapDF, dfMAST ,sortedCols , featureDir ):
    colList = []
    colCount = len(sortedCols)
    yColRange = [[i * 5 for i in range(colCount)]]
    for col in sortedCols[::-1]:
        colName = col[0]
        colList.append(colName)
    shapJSON = shapDF.to_dict(orient="records")
    featJSON = dfMAST.to_dict(orient="records")

    html = f"""
    <!DOCTYPE html>
    <html>
    <meta charset="UTF-8">
    <head>
        <title>Interactive SHAP for Features </title>
        <script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script>
        <script src="https://unpkg.com/simple-statistics@7.8.3/dist/simple-statistics.min.js"></script>
        <script>
            const jsonFeats = {featJSON};
            const jsonSHAP = {shapJSON};
            const globalMin = Math.min(...jsonFeats.flatMap(p => {colList}.map(c => p[c])));
            const globalMax = Math.max(...jsonFeats.flatMap(p => {colList}.map(c => p[c])));
            const figTrace = [];
            function plotData(){{
                for (let i = 0; i < {colList}.length; i++){{
                    const cols = {colList}[i];
                    const yVals = {yColRange}[i];
                    const xData = jsonSHAP.map(p => p[cols]);
                    const featVals = jsonFeats.map(p => p[cols]);
                    const yData = new Array(xData.length).fill(cols);

                    const trace = {{
                        x: xData,
                        y: yData,
                        name: cols,
                        mode: 'markers',
                        type: 'scatter',
                        customdata: jsonSHAP.map((p, i) => ({{
                            id: p.ID || `point_${{i}}`, 
                            image: p.base64 || null
                        }})),
                        hovertemplate:
                            `<b>${{cols}}</b><br>` +
                            `SHAP value: %{{x:.3f}}<br>` +
                            `Feature value: %{{marker.color:.3f}}<extra></extra>`,
                        marker: {{
                            size: 8,
                            color: featVals,
                            colorscale: 'RdBu',
                            cmin: globalMin,
                            cmax: globalMax,
                            showscale: false,
                        }}
                    }};
                    figTrace.push(trace);
                }}
                const featureList = document.getElementById("featureList");
                const traceVisible = new Array(figTrace.length).fill(true);
                figTrace.forEach((trace, idx)=>{{
                    const button = document.createElement("button");
                    button.className = "featureButton";
                    button.textContent = trace.name;
                    button.dataset.trace = idx;
                    featureList.appendChild(button);
                }});
                figTrace.push({{
                    x: [null],
                    y: [null],
                    mode: "markers",
                    hoverinfo: "skip" , 
                    showlegend: false,
                    marker: {{
                        size: 0,
                        color: [globalMin],
                        colorscale: "RdBu",
                        cmin: globalMin,
                        cmax: globalMax,
                        showscale:true,
                        colorbar: {{title: "Feature Value"}}
                    }}
                }});

            const layout = {{
                title: 'Interactive SHAP Plot',
                xaxis: {{
                    title: {{text: 'SHAP Value', font: {{family: 'Arial', size: 16, weight: 'bold'}}}},
                    tickfont: {{family: 'Arial', size: 14, weight: 'bold'}},
                    linecolor: 'black',
                    linewidth: 2,
                    mirror: true,
                    showgrid: true,
                    zeroline: false
                }},
                yaxis:{{
                    type:"category",
                    side:"left",
                    showticklabels:false,
                    autorange:"reversed",
                    automargin:true,
                    linecolor:"black",
                    linewidth:2,
                    mirror:true
                }},
                hovermode: 'closest',
                showlegend: false,
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
            Plotly.newPlot('plotContainer', figTrace, layout, config).then(() => {{
                const featureList = document.getElementById("featureList");
                const visible = new Array(figTrace.length).fill(true);
                featureList.addEventListener("click" , e=>{{
                    if(!e.target.matches(".featureButton"))
                        return;
                    const idx = Number(e.target.dataset.trace);
                    visible[idx] = !visible[idx];
                    Plotly.restyle(
                        "plotContainer",
                        {{
                            visible: visible[idx] ? true : "legendonly"
                        }},
                        [idx]
                    );
                    e.target.classList.toggle("inactive");
                }});
                const plotElement = document.getElementById('plotContainer');
                const imageContainer = document.getElementById('imageContainer');
                
                plotElement.on('plotly_hover', function(data) {{
                    const pointData = data.points[0];
                    const customData = pointData.customdata;
                    
                    if (customData && customData.image) {{
                        const imageSrc = customData.image.startsWith('data:') ? 
                            customData.image : 
                            `data:image/png;base64,${{customData.image}}`;
                        
                        imageContainer.innerHTML = `
                            <h3>Hover Image</h3>
                            <p><strong>${{pointData.customdata.id}}</strong></p>
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
                    imageContainer.innerHTML = `
                        <h3>Hover Image</h3>
                        <div class="image-placeholder">Hover over a point to see its image</div>
                    `;
                }});
            }});
            }}
        document.addEventListener("DOMContentLoaded", plotData);
        </script>
        <style>
            body{{
                font-family: Arial, sans-serif;
                margin: 20px;
                background-color: #f5f5f5;
            }}
            .plot-container{{
                display:flex;
                gap:20px;
                align-items:flex-start;
            }}

            #featureList{{
                width:180px;
                display:flex;
                flex-direction:column;
                gap:3px;
            }}

            .featureButton{{

                background:white;

                border:1px solid #ccc;

                border-radius:3px;

                cursor:pointer;

                text-align:right;

                padding:4px 8px;

                font-size:13px;

                transition:.15s;
            }}

            .featureButton:hover{{
                background:#eee;
            }}

            .featureButton.inactive{{

                color:#999;

                text-decoration:line-through;

                background:#f7f7f7;
            }}

            #plotContainer{{

                flex:1;

                min-height:700px;

                background:white;
            }}

            #imageContainer{{

                width:300px;
            }}
        </style>
    </head>
    <body>
        <h1>Interactive SHAP Plot</h1>

        <div class="plot-container">

            <div id="featureList"></div>

            <div id="plotContainer"></div>

            <div id="imageContainer">
                <h3>Hover Image</h3>
                <div class="image-placeholder">
                    Hover over a point to see its image
                </div>
            </div>

        </div>
    </body>
    </html>
    """
    with open(featureDir / "shapInteractive.html", "w", encoding="utf-8") as f:
        f.write(html)
def main(df , yCol  , idStr , SMILESStr , outputDir):
    while True:
        modelInt = input(f"Please Select the index for the model you want to fit to:\n\n[1]  XGBOOST\n\n[2]   RandomForest\n\n[3]    Linear Regression\n\n[4]    MLP\n\n").strip()
        if modelInt == "1":
            model = GradientBoostingRegressor(n_estimators=300 ,max_depth=4 , learning_rate=0.05 )
            break
        elif modelInt == "2" : 
            model = RandomForestRegressor(n_estimators=300,max_depth=4,random_state=42,n_jobs=-1)
            break
        elif modelInt == "3":
            model = LinearRegression()
            break
        elif modelInt == "4":
            model = MLPRegressor(activation = 'relu' , alpha =  0.01, hidden_layer_sizes = (100, 50), learning_rate_init =  0.005)
            break
        else:
            print("Invalid input, enter 1 or 2 only")

    yVals = df[yCol]
    idVals = df[idStr]    
    figDir = Path(outputDir) / "png" 
    dfMAST = createPNGDF(df ,"SMILES" , str(figDir))
    base64Col = []
    for img in list(dfMAST["pngPath"]):

        base64 = png64(img)
        base64Col.append(base64)
    dfMAST = dfMAST.drop(columns = [yCol , idStr , SMILESStr , "pngPath"])
    scaler = StandardScaler()

    dfMAST = pd.DataFrame(
        scaler.fit_transform(dfMAST),
        columns=dfMAST.columns,
        index=dfMAST.index
    )
    if modelInt == "3":
        model.fit(dfMAST , yVals)
        modelExplainer = shap.LinearExplainer(model , dfMAST)
    elif modelInt == "4":
        model.fit(dfMAST , yVals )
        modelExplainer = shap.KernelExplainer(model, dfMAST)
    else:
        model.fit(dfMAST , yVals)
        modelExplainer = shap.TreeExplainer(model)
    shapValues = modelExplainer.shap_values(dfMAST)
    shapDF = pd.DataFrame(shapValues,index=dfMAST.index,columns=dfMAST.columns)
    shutil.rmtree(figDir)

    featureColsHash = {}
    for idx , feature in enumerate(shapDF.columns):
        shapVal = np.mean([np.abs(val) for val in list(shapDF[feature])])
        featureColsHash[feature] = shapVal
    sortedCols = sorted(featureColsHash.items(),key=lambda kv: kv[1])
    print(shapDF.head())
    print(dfMAST.head())
    shapDF["ID"] = idVals
    shapDF["base64"] = base64Col
    featureDir = Path(outputDir) / "features" 
    featureDir.mkdir(parents=True, exist_ok=True)

    generateInteractiveSHAP(shapDF, dfMAST ,sortedCols , featureDir)

if __name__ == "__main__":
    featureDFStr = str(sys.argv[1])
    yStr = str(sys.argv[2])
    idStr = str(sys.argv[3])
    smilesStr = str(sys.argv[4])
    outputDir = str(sys.argv[5])
    featureDF = pd.read_csv(featureDFStr)
    main(featureDF , yStr , idStr , smilesStr , outputDir )