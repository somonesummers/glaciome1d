import numpy as np
import os
import sys
# from glaciome1D_dimensional import glaciome, basic_figure, plot_basic_figure, constants
sys.path.append('/Users/psummers8/Documents/glaciome1D')
from glaciome1D import glaciome, basic_figure, plot_basic_figure, constants
from scipy.integrate import trapz
import pickle
import glob

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: 
    Paul Summers
    May 2025
    Example of running glaciome1d

"""


# basic parameters needed for setting up the model; later will modify this so that 
# the fjord geometry can be passed through
constant = constants()

n_pts = 21 # number of grid points
L = 1e4 # ice melange length
Ut = 0.6e4 # glacier terminus velocity [m/a]; treated as a constant
Uc = 0.6e4 # glacier calving rate [m/a]; treated as a constant
Ht = 600 # terminus thickness
n = 101 # number of time steps
dt = 0.01# 1/(n_pts-1)/10 # time step [a]; needs to be quite small for this to work

# specifying fjord geometry
X_fjord = np.linspace(-200e3,200e3,101)
Wt = 4000
W_fjord = Wt + 0/10000*X_fjord

B_const = -0.6*constant.daysYear


## Run to steady state
if(False):
    data = glaciome(n_pts, dt, L, Ut, Uc, Ht, B_const, X_fjord, W_fjord)

    data.steadystate(method='hybr')
    data.save('steadystate.pickle')

files = sorted(glob.glob('../steadystate.pickle'))

for j in np.arange(0,len(files)):
    file = open(files[j], 'rb')
    data = pickle.load(file)
    file.close()

#Reset to time = 0
data.t = 0
dt = 10.0/365.0 # 10 days
data.dt = dt

for i in range(37): 
    data.X_externalGrid = X_fjord
    data.B_externalGrid = np.ones_like(X_fjord) * (-0.6 - .2 * np.sin(i*2*np.pi/(36.5))) * constant.daysYear
    data.param.muS = 0.25 - .05 * np.sin(i*2*np.pi/(36.5))
    data.prognostic(method='hybr') # lm or hybr
    if(i % 1 == 0):
        print(f"Step {i:03d} with H0:{data.H0:7.2f} m and L:{data.L:9.2f} m")
        data.save(f'bVarMelt{i:05d}.pickle')

print("Done!")







