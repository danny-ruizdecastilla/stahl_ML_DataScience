#Run this code periodically to do a deep analysis on all the repo dependancies to create a graph of the file structure
import sys
import glob
import os
from pathlib import Path
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs

class TreeNode:
    def __init__(self, value):
        self.value = value
        self.children = []

    def add_child(self, node):
        self.children.append(node)
    
def gatherScripts(repoDir , scriptStrs):
    currentDir = Path(repoDir)
    connections = [] #overlooks the .git path
    depthDict = {}
    branching = {}
    directoryDict = {}
    while True:
        
        immediateSubdirs = []
        immediateFiles = []
        for p in currentDir.iterdir():
            if p.is_dir() and not p.name.startswith('.') and p not in connections:
                immediateSubdirs.append(p)
            elif any(str(p).endswith(ext) for ext in scriptStrs):
                immediateFiles.append(p)
        directoryDict[currentDir] = immediateFiles
        connections.append(currentDir)
        headFolder = currentDir.name
        depthDict[headFolder] = len(immediateSubdirs)
        branching[headFolder] = len(immediateSubdirs)

        if len(immediateSubdirs) > 0:
            nextDir = immediateSubdirs[0]
            connections.append(nextDir)
            currentDir = nextDir
        else:
            # If all branches are exhausted, break
            if all(value == 0 for value in depthDict.values()):
                break
            else:
                # Get folders with unvisited children
                remainingFolders = [k for k, v in branching.items() if v != 0]

                if not remainingFolders:
                    break  # failsafe in case something gets stuck

                # Choose one (you can change strategy here)
                nextFolderName = min(remainingFolders)  # alphabetically smallest

                # Find path to it among connections
                for visited in connections:
                    if visited.name == nextFolderName:
                        currentDir = visited
                        break
    return directoryDict
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
        
        for edge in import_:
            dependencies[edge.strip()] = [node , str(pathFile.stem)]
    return dependencies     
        


def main(repoDir , scriptStrs):
    scripts = gatherScripts(repoDir , scriptStrs)
    allFiles = []
    for key, val in scripts.items():
        if len(val) != 0:
            for file in val:
                print(file)
                allFiles.append(file)


if __name__ == "__main__":
    repoDir = str(sys.argv[1])
    scriptStrs = listInputs(f"Please enter all the file types for the scripts in this repository: {repoDir} | Ex: .py,.cpp,.html,.js")
    main(repoDir , scriptStrs)