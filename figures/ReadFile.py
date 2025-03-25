

import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import matplotlib.patheffects as PathEffects
import sys
from numpy import random

# sys.path.append('/hdd/glaciome/models/glaciome1D')
# sys.path.insert(0, '')
sys.path.append('/Users/psummers8/Documents/glaciome1D')
from glaciome1D import constants, glaciome


matplotlib.rc('lines',linewidth=1) 

font = {'family' : 'serif',
        'weight' : 'normal',
        'size'   : 8}

matplotlib.rc('font', **font)

import cmasher as cmr

cmap = cmr.get_sub_cmap('viridis', 0, 0.95)

import glob
import pickle

constant = constants()




files = glob.glob('./*.pickle')
with open(files[0], 'rb') as file:
    data = pickle.load(file)
    file.close()
print(files[0])
print('Glacier thickness %.3f m' %data.Ht)
print('Glacier velocity  %.3f m/a' %data.Ut)
print('Glacier velocity  %.3f m/a' %data.Uc)

