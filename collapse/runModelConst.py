import numpy as np
from matplotlib import pyplot as plt
import os
import sys
# from glaciome1D_dimensional import glaciome, basic_figure, plot_basic_figure, constants
sys.path.append('/Users/psummers8/Documents/glaciome1D')
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
from glaciome1D import glaciome, basic_figure, plot_basic_figure, constants
from scipy.integrate import trapz
import pickle
import glob

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: 
    Paul Summers
    August 2025
    Example of running glaciome1d, now with umelt shape

"""

def meltRate(strength,n):
    x = np.linspace(0,1,n)
    baseLine = .6 -.4*x
    arc = .5 + 50*(x-.46)**6

    min = .4
    max = .6587
    if(strength < min):
        return baseLine*strength/min
    elif(strength < max):
        r = (strength - min)/(max-min)
        return baseLine*(1-r) + arc*(r)
    else:
        return arc*strength/max

# basic parameters needed for setting up the model; later will modify this so that 
# the fjord geometry can be passed through
constant = constants()

n_pts = 21 # number of grid points
L = 15e3 # ice melange length
Ut = 0.6e4 # glacier terminus velocity [m/a]; treated as a constant
Uc = 0.6e4 # glacier calving rate [m/a]; treated as a constant
Ht = 600 # terminus thickness
n = 101 # number of time steps
dt = 0.01# 1/(n_pts-1)/10 # time step [a]; needs to be quite small for this to work

# specifying fjord geometry
X_fjord = np.linspace(-200e3,200e3,101)
Wt = 5600
W_fjord = Wt + 0/10000*X_fjord
B_const = -0.4 * constant.daysYear


## Run to steady state
if(False):
    data = glaciome(n_pts, dt, L, Ut, Uc, Ht, B_const, X_fjord, W_fjord)
    # files = sorted(glob.glob('steadystate.pickle'))
    # for j in np.arange(0,len(files)):
    #     file = open(files[j], 'rb')
    #     data = pickle.load(file)
    # # file.close()
    # data.X_externalGrid = data.X
    # data.B_externalGrid = -1*meltRate(B_const,n_pts) * constant.daysYear
    data.steadystate()
    data.save('steadystate.pickle')

files = sorted(glob.glob('chkpt035uMelt.pickle'))

for j in np.arange(0,len(files)):
    file = open(files[j], 'rb')
    data = pickle.load(file)
    file.close()

data.W_fjord = 5600 + 0*data.X_fjord

#Reset to time = 0
data.t = 0
dt = 10/365.0 # [years]
data.dt = dt

yearsToSimulate = 5
endTime = int(yearsToSimulate/dt)
time = np.arange(endTime)*dt
iList = np.arange(0,endTime,1) 
# print(iList)
Bview = np.linspace(.35,.55,endTime)
# print(Bview)
plt.plot(time,Bview,color='red')
plt.xlabel('Time [years]')
plt.ylabel('Average Melt Rate [m/day]')
plt.show()
plt.close()

for i in iList: # 5 years 
    data.X_externalGrid = data.X
    data.B_externalGrid = -1*meltRate(Bview[i],n_pts) * constant.daysYear
    data.prognostic(method='lm') # lm or hybr
    if(i % 1 == 0):
        print(f"Step {i:03d} with H0:{data.H0:7.2f} m and L:{data.L:9.2f} m")
        data.save(f'uMelt{i:05d}.pickle')

print("Done!")







