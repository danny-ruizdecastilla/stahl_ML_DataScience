#Danny Ruiz de Castilla 09.22.25
#takes in a smiles of a molecule and returns indexes of a specific sub molec
import sys
import glob
import numpy as np
import os
from rdkit import Chem
import pandas as pd
from rdkit.Chem.PandasTools import LoadSDF
import matplotlib.pyplot as plt
import random
import re
import networkx as nx
from networkx import Graph
from itertools import combinations
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)

def texturedAdjency(smiles):


    return adjencyMatrix