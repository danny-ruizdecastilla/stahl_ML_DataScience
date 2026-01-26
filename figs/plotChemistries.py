import pandas as pd
import os 
import sys
import glob
import json
import numpy as np
import html
import plotly
import shutil
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
from figs.chemPlotlyV2 import createPNGDF,png64
from figs.chemPlotlyV1 import convertCanonical
from figs.plotSubstrates import safeStringHTML , plotSubstratesMain
#from DFTWorkflow.pitchingATent import compressData , locateNans , eliminateNans , convertCanonical , featureFiltering
#from DFTWorkflow.featureMaping import  createCSV
def checkSubstratePath(substrateDF):
    if not os.path.exists(substrateDF):
        print(f"\nWarning: Could not locate master substrate dataset at:\n   → {substrateDF}")
        print("It looks like it doesn't exist yet.")
        print("Would you like to build it from scratch?")
        while True:
            userInput = input("➡️  Enter 1 to build from scratch, or 2 to abort: ").strip()
            if userInput == "1":
                print("Proceeding to build the dataset from scratch...\n")
                return 1
            elif userInput == "2":
                print("Aborting the process\n")
                return 2
            else:
                userInput = input("⚠️  Invalid input. Please enter 1 (yes) or 2 (no). ").strip()

    else:
        print(f"Found existing master dataset at: {substrateDF}")
        return 0
def createAxisMotifs(axisNum):
    axisDict = {}
    for num in range(1,axisNum+1):
        motifList = listInputs(f"Enter motif names for axis {num}, Ex: [distance,Buried,angle,dihedral,Vbur]: ")
        naming = input(f"Do you want to name this axis? | Enter 1 for yes, or 2 to accept {num} as the name for axis  {num}:")
        while True:
            if naming == "1":
                name = input(f"Enter the name for axis {num}:")
                break
            elif naming == "2":
                name = str(num)
                break
            else:
                naming = input("⚠️  Invalid input. Please enter 1 (yes) or 2 (no).").strip()
        axisDict[name] = motifList
    return axisDict
def htmlGenerator2(jsonDict, axisList, chemStr, outputDir, partitionStr):
    plotChemStr = input("Name your new functional scatterplot: ")
    htmlAxis = "".join([safeStringHTML(axis) for axis in axisList])
    
    all_values = []
    for group in jsonDict.values():
        all_values.extend(p[partitionStr] for p in group if partitionStr in p and isinstance(p[partitionStr], (int, float)))

    if all_values:
        sliderMin = min(all_values)
        sliderMax = max(all_values)
        sliderStart = round((sliderMin + sliderMax) / 2, 2)
    else:
        sliderMin = 0
        sliderMax = 100
        sliderStart = 50
    
    html = f"""
    <!DOCTYPE html>
    <html>
    <meta charset="UTF-8">
    <head>
        <title>Scatter Plot with Hover Images</title>
        <script src="https://cdn.plot.ly/plotly-3.0.1.min.js"></script>
        <script src="https://unpkg.com/simple-statistics@7.8.3/dist/simple-statistics.min.js"></script>
        <script>
            const jsonDict = {json.dumps(jsonDict)};
            const groupedData = jsonDict;

            function plotData() {{
                const selectedGroup = document.getElementById("groupDropDown").value;
                const threshold = parseFloat(document.getElementById("yieldSlider").value);
                const xCol = document.getElementById("xAxis").value;
                const yCol = document.getElementById("yAxis").value;

                // Update threshold display
                document.getElementById("thresholdValue").textContent = threshold;

                let traces = [];
                let allXData = [];
                let allYData = [];

                for (const groupKey in groupedData) {{
                    const jsonData = groupedData[groupKey];
                    
                    if (groupKey === selectedGroup) continue;
                    const xData = jsonData.map(p => p[xCol]);
                    const yData = jsonData.map(p => p[yCol]);
                    allXData.push(...xData);
                    allYData.push(...yData);
                    
                    const trace = {{
                        x: xData,
                        y: yData,
                        mode: 'markers',
                        type: 'scatter',
                        name: `${{groupKey}} (background)`,
                        text: jsonData.map(p => `"{partitionStr}": ${{p["{partitionStr}"]}}`),
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
                            size: 8,
                            color: 'rgba(188, 188, 188, 0.4)'
                        }}
                    }};
                    traces.push(trace);
                }}
                const selectedData = groupedData[selectedGroup] || [];
                const above = selectedData.filter(p => p["{partitionStr}"] > threshold);
                const below = selectedData.filter(p => p["{partitionStr}"] <= threshold);
                if (below.length > 0) {{
                    const xData = below.map(p => p[xCol]);
                    const yData = below.map(p => p[yCol]);
                    allXData.push(...xData);
                    allYData.push(...yData);
                    
                    const traceBelow = {{
                        x: xData,
                        y: yData,
                        mode: 'markers',
                        type: 'scatter',
                        name: `${{selectedGroup}} ("{partitionStr}" ≤ ${{threshold}})`,
                        text: below.map(p => `"{partitionStr}": ${{p["{partitionStr}"]}}`),
                        customdata: below.map((p, i) => ({{
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
                            color: 'rgba(255, 255, 255, 0.5)',
                            line: {{
                                color: 'rgba(153, 000, 000, 1.0)',
                                width: 1
                            }}
                        }},
                    }};
                    traces.push(traceBelow);
                }}
                if (above.length > 0) {{
                    const xData = above.map(p => p[xCol]);
                    const yData = above.map(p => p[yCol]);
                    allXData.push(...xData);
                    allYData.push(...yData);
                    
                    const traceAbove = {{
                        x: xData,
                        y: yData,
                        mode: 'markers',
                        type: 'scatter',
                        name: `${{selectedGroup}} ("{partitionStr}" > ${{threshold}})`,
                        text: above.map(p => `"{partitionStr}": ${{p["{partitionStr}"]}}`),
                        customdata: above.map((p, i) => ({{
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
                            color: 'rgba(001, 031, 091, 0.8)',
                        }}
                    }};
                    traces.push(traceAbove);
                    
                }}

                // Linear regression calculation
                let lineTrace = null;
                let equationAnnotation = null;
                
                if (allXData.length > 1 && document.getElementById("toggleFit").checked) {{
                    try {{
                        const dataPoints = allXData.map((x, i) => [x, allYData[i]]).filter(([x, y]) => 
                            !isNaN(x) && !isNaN(y) && isFinite(x) && isFinite(y)
                        );
                        
                        if (dataPoints.length > 1) {{
                            const linearRegression = ss.linearRegression(dataPoints);
                            const linearFn = ss.linearRegressionLine(linearRegression);
                            const r2 = ss.rSquared(dataPoints, linearFn);
                            
                            const xFit = [Math.min(...allXData), Math.max(...allXData)];
                            const yFit = xFit.map(x => linearFn(x));
                            const slope = linearRegression.m.toFixed(3);
                            const intercept = linearRegression.b.toFixed(3);
                            const yTrue = dataPoints.map(([x, y]) => y);
                            const yPred = dataPoints.map(([x, y]) => linearFn(x));
                            const rmse = Math.sqrt(ss.mean(yTrue.map((y, i) => Math.pow(y - yPred[i], 2))));
                            const rmseText = rmse.toFixed(4);
                            const r2Text = r2.toFixed(4);
                            const equation = `y = ${{slope}}x + ${{intercept}}<br> RMSE = ${{rmseText}}`;

                            lineTrace = {{
                                x: xFit,
                                y: yFit,
                                mode: 'lines',
                                type: 'scatter',
                                name: 'Fit Line',
                                line: {{ dash: 'dot', width: 2, color: 'red' }},
                                showlegend: false
                            }};
                            
                            equationAnnotation = {{
                                x: 0.05,
                                y: 0.95,
                                xref: 'paper',
                                yref: 'paper',
                                text: equation,
                                showarrow: false,
                                align: 'left',
                                font: {{ size: 14 }},
                                bordercolor: 'blue',
                                borderwidth: 1,
                                borderpad: 4,
                                bgcolor: 'white',
                                opacity: 0.9
                            }};
                        }}
                    }} catch (error) {{
                        console.warn('Linear regression failed:', error);
                    }}
                }}
                
                if (lineTrace) {{
                    traces.push(lineTrace);
                }}

                const layout = {{
                    title: 'Literature derived Epoxidation Substrate Space',
                    xaxis: {{
                        title: {{text: xCol, font: {{family: 'Arial', size: 16, weight: 'bold'}}}},
                        tickfont: {{family: 'Arial', size: 14, weight: 'bold'}},
                        linecolor: 'black',
                        linewidth: 2,
                        mirror: true,
                        showgrid: true,
                        zeroline: false
                    }},
                    yaxis: {{
                        title: {{text: yCol, font: {{family: 'Arial', size: 16, weight: 'bold'}}}},
                        tickfont: {{family: 'Arial', size: 14, weight: 'bold'}},
                        linecolor: 'black',
                        linewidth: 2,
                        mirror: true,
                        showgrid: true,
                        zeroline: false
                    }},
                    hovermode: 'closest',
                    showlegend: true,
                    annotations: equationAnnotation ? [equationAnnotation] : []
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

                Plotly.newPlot('plotContainer', traces, layout, config).then(() => {{
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
                        imageContainer.innerHTML = `
                            <h3>Hover Image</h3>
                            <div class="image-placeholder">Hover over a point to see its image</div>
                        `;
                    }});
                }});
            }}

            window.onload = function() {{
                // Populate group dropdown
                const dropdown = document.getElementById("groupDropDown");
                for (const groupName in groupedData) {{
                    const option = document.createElement("option");
                    option.value = groupName;
                    option.textContent = groupName;
                    dropdown.appendChild(option);
                }}

                // Set initial selections and plot
                const xAxisSelect = document.getElementById("xAxis");
                const yAxisSelect = document.getElementById("yAxis");
                const yieldSlider = document.getElementById("yieldSlider");
                const toggleFit = document.getElementById("toggleFit");
                
                if (xAxisSelect && yAxisSelect) {{
                    xAxisSelect.selectedIndex = 0;
                    yAxisSelect.selectedIndex = 1;
                    plotData();
                    
                    // Add event listeners
                    xAxisSelect.addEventListener("change", plotData);
                    yAxisSelect.addEventListener("change", plotData);
                    dropdown.addEventListener("change", plotData);
                    yieldSlider.addEventListener("input", plotData);
                    toggleFit.addEventListener("change", plotData);
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
            .control-row {{
                display: flex;
                align-items: center;
                margin-bottom: 15px;
                flex-wrap: wrap;
                gap: 15px;
            }}
            select {{
                padding: 8px;
                margin: 0 10px;
                border: 1px solid #ccc;
                border-radius: 4px;
                font-size: 14px;
                min-width: 120px;
            }}
            label {{
                font-weight: bold;
                margin-right: 5px;
            }}
            .slider-container {{
                display: flex;
                align-items: center;
                gap: 10px;
            }}
            .linRegChk {{
                display: flex;
                align-items: center;
                gap: 10px;
            }}
            input[type="range"] {{
                width: 200px;
            }}
            #thresholdValue {{
                font-weight: bold;
                color: #007acc;
                min-width: 40px;
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
        <h1>{plotChemStr}</h1>

        <div class="controls">
            <div class="control-row">
                <label for="groupDropDown">Dataset:</label>
                <select id="groupDropDown">
                    <!-- Options will be populated by JavaScript -->
                </select>

                <label for="xAxis">X-axis:</label>
                <select id="xAxis">
                    {htmlAxis}
                </select>

                <label for="yAxis">Y-axis:</label>
                <select id="yAxis">
                    {htmlAxis}
                </select>
            </div>
            
            <div class="control-row">
                <div class="slider-container">
                    <label for="yieldSlider">"{partitionStr}" Threshold:</label>
                    <input type="range" id="yieldSlider" min="{sliderMin}" max="{sliderMax}" value="{sliderStart}" step="1">
                    <span id="thresholdValue">50</span>%
                </div>
            </div>
            <div class="control-row">
                <div class="linRegChk">
                    <label>
                        <input type="checkbox" id="toggleFit" checked/>
                        Show Fit line
                    </label>
                </div>
            </div>
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
    
    with open(outputDir + "/" + chemStr + "interactiveMASTER.html", "w", encoding="utf-8") as f:
        f.write(html)
def getBackgrounds(df , smileList , smileStr):
    return df[~df[smileStr].isin(smileList)].copy()
def main(chemistryFiles, masterDF, chemistriesDict , chemistry, outputDir , partitionStr):
    cols = masterDF.columns
    if "base64" in cols:
        pass
    else:
        figDir = input("Enter the figure Directory where PNG's of the substrates will be stored (temporarily): ")
        masterDF = createPNGDF(masterDF ,"SMILES" , figDir)
        base64Col = []
        for img in list(masterDF["pngPath"]):

            base64 = png64(img)
            base64Col.append(base64)
        masterDF["base64"] = base64Col
        masterDF = masterDF.drop(columns=['pngPath'])
        if "Unnamed: 0" in masterDF.columns:
            masterDF = masterDF.drop(columns=['Unnamed: 0'])
        shutil.rmtree(figDir)
    jsonDict = {}
    allSmiles = []
    for file in chemistryFiles:
        chemName = file.split("/")[-1].split(".")[0]
        columnList = list(masterDF.columns)
        if not partitionStr in columnList:
            columnList.append(partitionStr)
        df = pd.DataFrame(columns= columnList)
        chemDF = pd.read_csv(file)  
        for _, row in chemDF.iterrows():
            smiles = row["SMILES"]
            canonical = convertCanonical(smiles)
            matches = masterDF[masterDF["canonicalSMILES"] == canonical]
            if not matches.empty:
                allSmiles.append(canonical)
                rowMatch = matches.head(1).copy()
                rowMatch[partitionStr] = float(row[partitionStr])
                df = pd.concat([df, rowMatch], ignore_index=True)
            else:
                continue
        jsonChem = df.to_json(orient="records" , force_ascii=False )
        jsonChem = jsonChem.replace('\\/', '/')
        dfName = str(chemistriesDict[chemName])
        jsonDict[dfName] = json.loads(jsonChem)
    allSmiles = list(set(allSmiles))
    masterDF = getBackgrounds(masterDF , allSmiles , "canonicalSMILES")
    jsonRemain = masterDF.to_json(orient='records' , force_ascii= False)
    jsonRemain = jsonRemain.replace('\\/', '/')
    jsonDict['Background Molecules'] = json.loads(jsonRemain)
    axisList = masterDF.select_dtypes(include='number').columns
    htmlGenerator2(jsonDict , axisList , chemistry , outputDir, partitionStr)


if __name__ == "__main__":
    chemistriesDir = str(sys.argv[1])
    substrateFile = str(sys.argv[2])
    chemistry = str(sys.argv[3])
    outputDir = str(sys.argv[4])
    partitionStr = str(sys.argv[5])
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    initDataSets = glob.glob(chemistriesDir + "/*.csv")
    initdataSets = sorted(initDataSets)
    chemistryNames = [name.split("/")[-1].split(".")[0] for name in initdataSets]
    chemistriesDict = {}
    for name in chemistryNames:
        nameList = listInputs(f"For Chemistry: {name} Please type the name you want this chemistry to be represented by: ")
        name_ = nameList[0]
        chemistriesDict[name] = name_
    #print(chemistriesDict)
    dataOption = checkSubstratePath(substrateFile)
    if dataOption == 2:
        sys.exit()
    elif dataOption == 1:
        datasetsDir = input("Enter the directory where the substrate Data resides: ")
        eliminatedFile = listInputs(f"Enter the dataframe eliminated phrases for {chemistry} chemistry: ")
        elimFile = eliminatedFile[0]
        if os.path.exists(elimFile):
            with open(elimFile, 'r') as file:
                content = file.read()
                eliminatedPhrases = [item.strip() for item in content.split(',') if item.strip()]
        else: 
            eliminatedPhrases = ["SMILES" , "Compound_Name", "Yield", "ChemistryType",  "Unnamed" , "ID" , "Canonicals" , "canonicalSMILES"]
        figDir = input("Enter the figure Directory: ")
        if not os.path.exists(figDir):
            os.makedirs(figDir)
        axisMotifs = createAxisMotifs(2)        
        masterDF = plotSubstratesMain(datasetsDir,chemistry , figDir  , axisMotifs, eliminatedPhrases , outputDir)
    elif dataOption== 0:
        masterDF = pd.read_csv(substrateFile , encoding='utf-8')
    if not "canonicalSMILES" in list(masterDF.columns):
        smilesList = masterDF["SMILES"]
        canonicals = []
        for smiles in smilesList:
            canonical = convertCanonical(smiles)
            canonicals.append(canonical)
        masterDF["canonicalSMILES"] = canonicals
    main(initDataSets , masterDF , chemistriesDict, chemistry, outputDir , partitionStr)