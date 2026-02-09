#3dscatterplot
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
def htmlGenerator3D(jsonDict, axisList, chemStr, outputDir, partitionStr):
    from figs.plotSubstrates import safeStringHTML
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
                const zCol = document.getElementById("zAxis").value;

                // Update threshold display
                document.getElementById("thresholdValue").textContent = threshold;

                let traces = [];
                let allXData = [];
                let allYData = [];
                let allZData = [];

                for (const groupKey in groupedData) {{
                    const jsonData = groupedData[groupKey];
                    
                    if (groupKey === selectedGroup) continue;
                    const xData = jsonData.map(p => p[xCol]);
                    const yData = jsonData.map(p => p[yCol]);
                    const zData = jsonData.map(p => p[zCol]);
                    allXData.push(...xData);
                    allYData.push(...yData);
                    allZData.push(...zData);
                    
                    const trace = {{
                        x: xData,
                        y: yData,
                        z: zData,
                        mode: 'markers',
                        type: 'scatter3d',
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
                            `${{zCol}}: %{{z}}<br>` +
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
                    const zData = below.map(p => p[zCol]);
                    allXData.push(...xData);
                    allYData.push(...yData);
                    allZData.push(...zData);
                    
                    const traceBelow = {{
                        x: xData,
                        y: yData,
                        z: zData,
                        mode: 'markers',
                        type: 'scatter3d',
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
                            `${{zCol}}: %{{z}}<br>` +
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
                    const zData = above.map(p => p[zCol]);
                    allXData.push(...xData);
                    allYData.push(...yData);
                    allZData.push(...zData);
                    
                    const traceAbove = {{
                        x: xData,
                        y: yData,
                        z: zData,
                        mode: 'markers',
                        type: 'scatter3d',
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
                            `${{zCol}}: %{{z}}<br>` +
                            '<extra></extra>',
                        marker: {{
                            size: 12,
                            color: 'rgba(001, 031, 091, 0.8)',
                        }}
                    }};
                    traces.push(traceAbove);
                    
                }}

                const layout = {{
                title: 'Literature derived Epoxidation Substrate Space',
                scene: {{
                    xaxis: {{
                    title: {{ text: xCol, font: {{ family: 'Arial', size: 16 }} }},
                    tickfont: {{ family: 'Arial', size: 14 }},
                    showgrid: true,
                    zeroline: false
                    }},
                    yaxis: {{
                    title: {{ text: yCol }},
                    showgrid: true
                    }},
                    zaxis: {{
                    title: {{ text: zCol }},
                    showgrid: true
                    }}
                }},
                hovermode: 'closest',
                showlegend: true
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
                const zAxisSelect = document.getElementById("zAxis");
                const yieldSlider = document.getElementById("yieldSlider");
                
                if (xAxisSelect && yAxisSelect && zAxisSelect) {{
                    xAxisSelect.selectedIndex = 0;
                    yAxisSelect.selectedIndex = 1;
                    zAxisSelect.selectedIndex = 2;
                    plotData();
                    
                    // Add event listeners
                    xAxisSelect.addEventListener("change", plotData);
                    yAxisSelect.addEventListener("change", plotData);
                    zAxisSelect.addEventListener("change", plotData);
                    dropdown.addEventListener("change", plotData);
                    yieldSlider.addEventListener("input", plotData);
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

                <label for="zAxis">Z-axis:</label>
                <select id="zAxis">
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
    
    with open(outputDir + "/" + chemStr + "3dinteractiveMASTER.html", "w", encoding="utf-8") as f:
        f.write(html)