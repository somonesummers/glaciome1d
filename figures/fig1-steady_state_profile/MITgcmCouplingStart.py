#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 27 12:39:00 2023

@author: jason
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 19 19:15:01 2023

@author: jason
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 16 14:04:35 2023

@author: jason
"""

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

#%%
run_simulations = 'y'

if run_simulations == 'y':
        
    n_pts = 21 # number of grid points
    L = 25e3 # ice melange length
    Ut = 6e3 # glacier terminus velocity [m/a]; treated as a constant
    Uc = 6e3 # glacier calving rate [m/a]; treated as a constant
    Ht = 800 # terminus thickness
    n = 101 # number of time steps
    dt = 0.002# 1/(n_pts-1)/10 # time step [a]; needs to be quite small for this to work
    H0_manual = 180 #posit H0 start value, set to None to use default (75 m)

    # specifying fjord geometry
    X_fjord = np.linspace(0e3,400e3,101)
    Wt = 4000
    W_fjord = Wt + 0/10000*X_fjord
    B = -0.6*constant.daysYear
    
    # first run to steady state
    # data = glaciome(n_pts, dt, L, Ut, Uc, Ht, B, X_fjord, W_fjord, H0 = H0_manual)
    files = glob.glob('./tempfile.pickle')
    with open(files[0], 'rb') as file:
        data = pickle.load(file)
        file.close()

    b_mitgcm = np.load('./meltRatesB.npy')*.4
    x_mitgcm = np.load('./meltRatesX.npy')
    print("Melt Rates [m/day]:", [b_mitgcm[0],b_mitgcm[-1]])
    # plt.figure()
    # plt.plot(x_mitgcm,b_mitgcm)
    # plt.show()
    # plt.close()
    data.X_externalGrid = x_mitgcm
    data.B_externalGrid = b_mitgcm*constant.daysYear
    data.dt = 0.0005
    data.steadystate(method='lm')
    data.save('MITgcmRun_0_5.pickle')

