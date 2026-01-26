import pandas as pd
import os 
import sys
import glob
import re
import numpy as np
from pathlib import Path
import torch
import torch.nn as nn
import torch.optim as optim
import torchvision
import torchvision.transforms as transforms
import torch.nn.functional as F
parentDir = Path(__file__).resolve().parents[1]
sys.path.append(str(parentDir))
from dimensionalityReduction.reactivityFeatures import boxGen
from reaxysProcessing.reaxysSubstrateExtractorV2 import listInputs
#Danny Ruiz de Castilla 01.23.26
class Net(nn.Module): #initialize neural net object 
    def __init__(self):
        super(Net, self).__init__()
        self.conv1 = nn.Conv2d(3, 6, 5)
        self.pool = nn.MaxPool2d(2, 2)
        self.conv2 = nn.Conv2d(6, 16, 5)
        self.fc1 = nn.Linear(16 * 5 * 5, 120)
        self.fc2 = nn.Linear(120, 84)
        self.fc3 = nn.Linear(84, 10)

    def forward(self, x):
        x = self.pool(F.relu(self.conv1(x)))
        x = self.pool(F.relu(self.conv2(x)))
        x = x.view(-1, 16 * 5 * 5)
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = self.fc3(x)
        return x
def main():
    triplet_loss = nn.TripletMarginLoss(margin=1.0, p=2, eps=1e-7)