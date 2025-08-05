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
    connections = []
    depthDict = {}
    branching = {}

    while True:
        
        # Get unvisited subdirectories
        immediateSubdirs = [p for p in currentDir.iterdir() if p.is_dir()]
        immediateSubdirs = [d for d in immediateSubdirs if d not in connections]
        connections.append(currentDir)
        headFolder = currentDir.name
        depthDict[headFolder] = len(immediateSubdirs)
        branching[headFolder] = len(immediateSubdirs)

        # If there are subdirectories, go deeper
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


def main(repoDir , scriptStrs):
    shellTree  , folderDict = gatherScripts(repoDir , scriptStrs)

if __name__ == "__main__":
    repoDir = str(sys.argv[1])
    scriptStrs = listInputs(f"Please enter all the file types for the scripts in this repository: {repoDir} | Ex: .py,.cpp,.html,.js")
    main(repoDir , scriptStrs)