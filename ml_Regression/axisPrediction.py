#Predicts the rank positioning on a 1-d axis based on a feature matrix
import torch
import torch.nn as nn
import os 
import sys
import glob
import pandas as pd 
import numpy as np
import chemdraw
import re
from pathlib import Path
