import numpy as np
from matplotlib import pyplot as plt
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
Wt = 5600
W_fjord = Wt + 0/10000*X_fjord

B_const = -0.3*constant.daysYear


## Run to steady state
if(False):
    # data = glaciome(n_pts, dt, L, Ut, Uc, Ht, B_const, X_fjord, W_fjord)
    files = sorted(glob.glob('steadystate_muS30_w560_sgd2000.pickle'))
    for j in np.arange(0,len(files)):
        file = open(files[j], 'rb')
        print(file)
        data = pickle.load(file)
        file.close()
    data.X_externalGrid = X_fjord
    data.B_externalGrid = np.ones_like(X_fjord) * -.6 * constant.daysYear
    print(data)
    data.steadystate(method='lm') # lm or hybr
    data.save('steadystate_B060.pickle')

if(True):
    #Reset to time = 0
    dt = 10/365.0 # 10 days
    yearsToSimulate = 10
    endTime = int(yearsToSimulate/dt)
    time = np.arange(endTime)*dt
    iList = np.arange(0,endTime,1)
    BSin = np.sin(iList*2*np.pi/(36.5)) - np.sqrt(3)/2 # ~60 days of year
    BSin[BSin < 0 ] = 0
    BSin = BSin/np.max(BSin)
    Bview = (-0.6 - .0 * BSin)
    muSin = np.sin(iList*2*np.pi/(36.5)) - 1/2 # ~120 days of year
    muSin[muSin < 0] = 0
    muSin = muSin/np.max(muSin)
    muSview = 0.3 - .00 * muSin

    plt.plot(time,Bview,color='red')
    ax_2 =plt.gca().twinx()
    ax_2.plot(time,muSview,color='aqua')
    plt.xlabel('Time [year]')
    # plt.show()
    plt.close()

    for alpha in [1e-6]:
        files = sorted(glob.glob('steadystate_muS30_w560_sgd2000.pickle'))
        for j in np.arange(0,len(files)):
            file = open(files[j], 'rb')
            data = pickle.load(file)
            file.close()
        data.t = 0
        data.dt = dt

        U0 = 7300 #unbuttressed calving [m/yr]
        # alpha =10.0e-5 #buttressing coefficient (5,10,15e-5 so far have been good) [m^2 yr^-1 N ^-1]
        beta = 50.0e-3 #calving rate coefficient [50m y^-1 meter^-1]
        print(f"alpha {alpha:.2e}, beta {beta:.2e}, F0 is {U0/alpha:3.2e} N/m")
        print(f"{iList[-1]} steps to run")
        for i in iList: 
            data.X_externalGrid = X_fjord
            data.B_externalGrid = np.ones_like(X_fjord) * Bview[i] * constant.daysYear
            data.param.muS = muSview[i]
            F = data.H0*data.pressure(data.H0)
            data.Uc = U0 - alpha * F - beta * data.X[0]
            data.Ht = 600 - data.X[0] * 50e-3 ## meter increase in thickness for 1km of retreat 
            data.prognostic(method='hybr') # lm or hybr
            if(i % 10 == 0):
                print(f"Step {i:03d} with H0:{data.H0:7.2f} m and L:{data.L:9.2f} m")
            if(i % 1 == 0):
                data.save(f'linearBed_a{int(alpha*1e6):03d}_b{int(beta*1e3):03d}_{i:05d}.pickle')

print("Done!")







