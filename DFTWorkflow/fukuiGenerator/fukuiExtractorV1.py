import os 
import sys
import glob
import shutil
import re
import pandas as pd 
#Danny Ruiz de Castilla
#Generates Fukui maps for each molecule 
def fukuiMap(substrate:dict):
    substrateDF = pd.DataFrame(columns = substrate.keys())
    for atom , presences in enumerate(substrate):
        #presences = [neutral , cation , anion]
        f_0 = 0.5*(float(presences[2]) - float(presences[1]) )
        f_plus = float(presences[2]) - float(presences[0])
        f_minus = float(presences[0]) - float(presences[1])

