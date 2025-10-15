import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
from rdkit.Chem.PandasTools import LoadSDF
import ast
import networkx as nx
from networkx import Graph
import matplotlib.pyplot as plt
from itertools import combinations
from collections import deque 
import random
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__) , "../.."))
sys.path.append(parentDir)

#Take a step
#Check if you've entered a cyclic
#If so, skip to edges, add 1 step
#If not then add 1 step
#Check termination
#When terminated return final submolec 
def isCyclic(mol, g , atomID , currentPath):
    cyclit = False
    