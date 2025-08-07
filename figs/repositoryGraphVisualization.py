#Run this code periodically to do a deep analysis on all the repo dependancies to create a graph of the file structure
import sys
import glob
import os
from pathlib import Path
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

        directoryDict[currentPath] = scriptFiles

    build_tree(repoDir)
    return root, directoryDict
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
    tree , scripts = gatherScripts(repoDir , scriptStrs)
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