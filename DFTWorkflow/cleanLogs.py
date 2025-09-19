#clean up log files before extracting information
#Danny Ruiz de Castilla 09.19.2025
import sys
import glob
import numpy as np
import os
import pandas as pd
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(parentDir)

def basicTerm(logFile, errorPhrase):
    with open(logFile, 'r') as f:
        for line in f:
            if errorPhrase in line:
                return True
    return False