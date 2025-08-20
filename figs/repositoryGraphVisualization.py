#Run this code periodically to do a deep analysis on all the repo dependancies to create a graph of the file structure
import sys
import json
import os
from pathlib import Path
from PIL import ImageColor
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs

class TreeNode:
    def __init__(self, name ):
        self.name = name
        self.children = []
    def add_child(self, name):
        self.children.append(name)
def gatherScripts(repoDir, scriptStrs):
    allScripts = []
    repoDir = Path(repoDir)
    root = TreeNode(repoDir.name)
    directoryDict = {}  # Path → list of files
    visited = set()
    nodeMap = {repoDir: root}

    def build_tree(currentPath):
        if currentPath in visited or currentPath.name.startswith('.'):
            return
        visited.add(currentPath)

        currentNode = nodeMap[currentPath]
        scriptFiles = []

        for p in currentPath.iterdir():
            if p.name.startswith('.'):
                continue 
            if p.is_dir():
                childNode = TreeNode(p.name)
                currentNode.add_child(childNode)
                nodeMap[p] = childNode
                build_tree(p)  # Recurse
            elif any(str(p).endswith(ext) for ext in scriptStrs):
                scriptFiles.append(p)
                script = Path(p)
                scriptName = script.name
                allScripts.append(scriptName)
        directoryDict[currentPath] = scriptFiles

    build_tree(repoDir)
    return root, directoryDict , allScripts
def str_to_rgb(colorStr):
    try:
        rgb = ImageColor.getrgb(colorStr)
        return rgb
    except:
        return False
def brighten(rgb, factor):
    return tuple(min(int(c * factor), 255) for c in rgb)
def colorNodes(root, nodesColors=None, notFirst=False, prevRGB=None):
    if nodesColors is None:
        nodesColors = {}
    
    firstNodes = [rnode.name for rnode in root.children]

    for node in firstNodes:
        if notFirst:
            # Brighten from parent's color
            newRGB = brighten(prevRGB, 1.2)
            nodesColors[node] = newRGB
        else:
            colorStr = input(f"Please enter the color name for this subdirectory {node}: ").strip()
            while True:
                rgb = str_to_rgb(colorStr)
                if not rgb:
                    colorStr = input(f"Please enter a new color name for this subdirectory {node}: ").strip()
                else:
                    break
            nodesColors[node] = rgb
    for rnode in root.children:
        colorNodes(
            rnode,
            nodesColors=nodesColors,
            notFirst=True,
            prevRGB=nodesColors[rnode.name]
        )
    return nodesColors
def readScript(file , scriptNames):
    eligibleLines = []
    with open(file, 'r') as f:
        for line in f:
            if " import " in line:
                eligibleLines.append(line)
    newLines = [line for line in eligibleLines if any(scriptName in line for scriptName in scriptNames)]
    pathFile = Path(file)
    dependencies = {}
    for line in newLines:
        import_ = line.split(" import ")[-1].split(",")

        if "." in line.split(" import ")[0]:
            node = line.split(" import ")[0].split(".")[-1].strip()
        else:
            node = line.split(" import ")[0].split("from")[-1].strip()
        allEdges = [edge.strip() for edge in import_]
        if len(allEdges) != 0:
            dependencies[node] = [edge.strip() for edge in import_]
    return dependencies   
from pathlib import Path

def networkGenerator(graphJSON, outputDir): 
    plotTitle = input("Name this new network graph: ")
    html = f"""
    <!doctype html>
    <html lang="en">
    <head>
        <meta charset="utf-8" />
        <meta name="viewport" content="width=device-width, initial-scale=1" />
        <title>Basic Network Graph (D3)</title>
        <script src="https://cdn.jsdelivr.net/npm/d3@7"></script>
        <style>
            :root {{
                --bg: #0f172a;
                --card: #11827e;
                --text: #ad0c0c;
                --muted: #94a3b8;
                --accent: #38bdf8;
            }}
            html, body {{
                height: 100%;
                margin: 0;
                background: radial-gradient(1200px 600px at 70% -10%, #060e19 0%, var(--bg) 50%);
                color: var(--text);
                font-family: ui-sans-serif, system-ui, Segoe UI, Roboto, Arial, Helvetica;
            }}
            .app {{
                max-width: 1100px;
                margin: 24px auto;
                padding: 10px;
                background: var(--card);
                border-radius: 12px;
                box-shadow: 0 20px 60px rgba(0,0,0,.50);
                border: 1px solid rgba(148,163,184,0.18);
            }}
            .header {{
                display: flex; align-items: center; justify-content: space-between; gap: 12px;
            }}
            .title {{
                margin: 0; font-size: 20px; letter-spacing: .2px; font-weight: 650;
            }}
            .sub {{ color: var(--muted); font-size: 12px; }}
            .btn {{
                appearance: none; 
                border: 1px solid rgba(148,163,184,.5); 
                background: transparent; 
                color: var(--text);
                padding: 8px 12px; 
                border-radius: 12px; 
                cursor: pointer; 
                font-weight: 600; 
                font-size: 12px; 
                letter-spacing: .2px;
            }}
            .btn:hover {{ border-color: rgba(148,163,184,.45); }}

            #graph {{
                width: 100%; 
                height: 620px; 
                border-radius: 12px; 
                overflow: hidden; 
                position: relative;
                background: linear-gradient(180deg, rgba(148,163,184,.5), rgba(148,163,184,.5));
            }}
            svg {{ width: 100%; height: 100%; display: block; }}

            .tooltip {{
                position: absolute; 
                pointer-events: none; 
                background: rgba(15,23,42,.92); 
                color: var(--text);
                border: 1px solid rgba(148,163,184,.25); 
                padding: 6px 8px; 
                border-radius: 10px; 
                font-size: 12px; 
                opacity: 0; 
                transform: translate(-50%,-140%);
                white-space: nowrap; 
                backdrop-filter: blur(6px);
            }}
        </style>
    </head>
    <body>
        <div class="app">
            <div class="header">
                <div>
                    <h1 class="title">{plotTitle}</h1>
                    <div class="sub"> Drag nodes • Scroll/Tap to zoom • Hover for details</div>
                </div>
                <div>
                    <button class="btn" id="shuffle">Shuffle layout</button>
                    <button class="btn" id="reset">Reset zoom</button>
                </div>
            </div>
            <div id="graph"></div>
        </div>
        <script>
            const graphData = {graphJSON};
            console.log(graphData);

            const container = document.getElementById('graph');
            const svg = d3.select(container).append('svg');
            const g = svg.append('g');
            const tooltip = d3.select(container).append('div').attr('class', 'tooltip');

            svg.append("defs").selectAll("marker")
            .data(["arrow"])
            .enter().append("marker")
                .attr("id", d => d)
                .attr("viewBox", "0 -5 10 10")
                .attr("refX", 20)
                .attr("refY", 0)
                .attr("markerWidth", 6)
                .attr("markerHeight", 6)
                .attr("orient", "auto")
                .append("path")
                    .attr("d", "M0,-5L10,0L0,5")
                    .attr("fill", "#999");

            const defs = svg.append('defs');
            const glow = defs.append('filter').attr('id','glow');
            glow.append('feGaussianBlur').attr('stdDeviation', 3.2).attr('result', 'coloredBlur');
            const feMerge = glow.append('feMerge');
            feMerge.append('feMergeNode').attr('in','coloredBlur');
            feMerge.append('feMergeNode').attr('in','sourceGraphic');
            
            const simulation = d3.forceSimulation(graphData.nodes)
                .force('link', d3.forceLink(graphData.edges).id(d => d.id).distance(100).strength(0.2))
                .force('charge', d3.forceManyBody().strength(-280))
                .force('collide', d3.forceCollide().radius(22))
                .force('center', d3.forceCenter())
                .force('x', d3.forceX().strength(0.05))
                .force('y', d3.forceY().strength(0.05));

            const link = g.selectAll('.link')
                .data(graphData.edges)
                .enter().append('line')
                .attr('class', 'link')
                .style('stroke', '#999')
                .style('stroke-opacity', 0.6)
                .style('stroke-width', d => d.weight || 1.5)
                .attr("marker-end", "url(#arrow)");
        
            const linkLabels = g.selectAll('.link-label')
                .data(graphData.edges)
                .enter().append('text')
                .attr('class', 'link-label')
                .attr('text-anchor', 'middle')
                .attr('dy', -5)
                .style('font-size', '12px')
                .style('fill', '#555')
                .style('opacity', 0) 
                .text(d => d.label);

            link.on('mousemove', function(event, d) {{
                d3.select(this).style('stroke-width', 3);
                tooltip.style('opacity', 1)
                    .html(`<strong>${{d.label}}</strong>`)
                    .style('left', (event.pageX + 10) + 'px')
                    .style('top', (event.pageY - 20) + 'px');
            }})
            .on('mouseout', function() {{
                d3.select(this).style('stroke-width', 1.5);
                tooltip.style('opacity', 0);
            }});

            const node = g.selectAll('.node')
                .data(graphData.nodes)
                .enter().append('circle')
                .attr('class', 'node')
                .attr('r', 10)
                .attr('fill', d => d.color)
                .style('filter','url(#glow)')
                .call(d3.drag()
                    .on('start', dragstarted)
                    .on('drag', dragged)
                    .on('end', dragended))
                .on('mousemove', (event, d) => {{
                    tooltip.style('opacity', 1)
                        .html(`<strong>${{d.id}}</strong><br/><span style="color:#94a3b8">${{d.group}}</span>`)
                        .style('left', (event.pageX + 10) + 'px')
                        .style('top', (event.pageY - 20) + 'px');
                }})
                .on('mouseout', () => tooltip.style('opacity', 0));

            simulation.on("tick", () => {{
                link
                    .attr("x1", d => d.source.x)
                    .attr("y1", d => d.source.y)
                    .attr("x2", d => d.target.x)
                    .attr("y2", d => d.target.y);

                node
                    .attr("cx", d => d.x)
                    .attr("cy", d => d.y);

                linkLabels
                    .attr("x", d => (d.source.x + d.target.x) / 2)
                    .attr("y", d => (d.source.y + d.target.y) / 2);
            }});

            const zoom = d3.zoom().scaleExtent([0.2, 4]).on('zoom', (event) => {{
                g.attr('transform', event.transform);
            }});
            svg.call(zoom);

            function resize() {{
                const {{ width, height }} = container.getBoundingClientRect();
                svg.attr('viewBox', `${{-width/2}} ${{-height/2}} ${{width}} ${{height}}`);
                simulation.force('center', d3.forceCenter(0, 0)).alpha(0.2).restart();
            }}
            window.addEventListener('resize', resize);
            resize();

            document.getElementById('shuffle').addEventListener('click', () => {{
                graphData.nodes.forEach(n => {{ n.vx = (Math.random()-.5)*3; n.vy = (Math.random()-.5)*3; }});
                simulation.alpha(0.6).restart();
            }});
            document.getElementById('reset').addEventListener('click', () => {{
                svg.transition().duration(350).call(zoom.transform, d3.zoomIdentity);
            }});

            function dragstarted(event, d) {{
                if (!event.active) simulation.alphaTarget(0.3).restart();
                d.fx = d.x; d.fy = d.y;
            }}
            function dragged(event, d) {{
                d.fx = event.x; d.fy = event.y;
            }}
            function dragended(event, d) {{
                if (!event.active) simulation.alphaTarget(0);
                d.fx = null; d.fy = null;
            }}
        </script>
    </body>
    </html>
    """
    
    networkHTMLPath = Path(outputDir) / f"{plotTitle}.html"
    with open(networkHTMLPath, "w", encoding="utf-8") as f:
        f.write(html)

def main(repoDir , scriptStrs , outputDir):
    trunkMAST = Path(repoDir).name
    tree , scriptsHash , scriptNames = gatherScripts(repoDir , scriptStrs) #gather code scripts in a hash table with directories 
    print(scriptsHash)
    edgeMAST = []
    treeColors = colorNodes(tree) #gives each level of the tree structure a color 
    print(treeColors)
    nodeMAST = []
    for key, val in scriptsHash.items():
        if len(val) != 0:
            mainFolder = str(key.name)
            if mainFolder in list(treeColors.keys()):
                color = treeColors[mainFolder]
                for file in val:
                    print(file)
                    fileName = file.name
                    fileStr = str(file)
                    firstDir = str(Path(fileStr.split(trunkMAST)[1]).parts[1])
                    node = { "id": fileName , "color": f"rgba({color[0]},{color[1]},{color[2]},0.9)", "group": firstDir}
                    nodeMAST.append(node)
                    dependenciesHash = readScript(fileStr, scriptNames)
                    for source , shared in dependenciesHash.values():
                        weight = len(shared)
                        allScripts = "\n".join(script for script in shared) + "\n"
                        edge = {"source" : source , "target": fileName , "label": allScripts , "weight" : 2*weight}
                        edgeMAST.append(edge)
    graphData = {"nodes": nodeMAST , "edges": edgeMAST}
    print(graphData)
    graphJSON = json.dumps(graphData)
    networkGenerator(graphJSON , outputDir)

if __name__ == "__main__":
    repoDir = str(sys.argv[1])
    outputDir = str(sys.argv[2])
    scriptStrs = listInputs(f"Please enter all the file types for the scripts in this repository: {repoDir} | Ex: .py,.cpp,.html,.js")
    main(repoDir , scriptStrs , outputDir)