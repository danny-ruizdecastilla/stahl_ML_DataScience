from dash import Dash, html, dcc, callback, Output, Input, callback_context
import plotly.graph_objects as go
import plotly 
import pandas as pd
import os 
import sys
import glob
import re
import chemdraw
import base64
import numpy as np
parentDir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(parentDir)
from figs.chemPlotlyV2 import createPNGDF , png64 , plotly_template
from DFTWorkflow.featureMaping import savePNG , createCSV
from DFTWorkflow.pitchingATent import compressData , locateNans , eliminateNans , convertCanonical
#Danny Ruiz de Castilla 05.02.25
