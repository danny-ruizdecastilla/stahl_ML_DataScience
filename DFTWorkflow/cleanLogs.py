#clean up log files before extracting information
#Danny Ruiz de Castilla 09.19.2025
import sys
import glob
import numpy as np
import os
import pandas as pd
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.append(parentDir)

def basicTerm(logFile, errorPhrase , termPhrase):
    termCount = 0
    with open(logFile, 'r') as f:
        for line in f:
            if errorPhrase in line:
                return True
            elif termPhrase in line:
                termCount +=1
    if termCount == 0:
        return True
    else:
        return False