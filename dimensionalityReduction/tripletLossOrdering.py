import pandas as pd
import os 
import sys
import glob
import re
import numpy as np
from pathlib import Path

parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from dimensionalityReduction.reactivityFeatures import boxGen
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
#Danny Ruiz de Castilla 01.23.26