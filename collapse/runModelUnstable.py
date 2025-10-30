import numpy as np
from matplotlib import pyplot as plt
import os
import sys
import time
# from glaciome1D_dimensional import glaciome, basic_figure, plot_basic_figure, constants
sys.path.append('/Users/psummers8/Documents/glaciome1D')
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
from glaciome1D import glaciome, basic_figure, plot_basic_figure, constants
import pickle
import glob
from scipy.integrate import simpson
sys.path.append('.')
# from localVars import *
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
    arc = .5 + 50*(x-.5)**6 #x-.46 was 
    minV = np.mean(baseLine) #0.4 is calibrated value
    maxV = np.mean(arc) - .015 #extra bit seems to help glaciome get right value for U shape
    if(strength < minV):
        # print(baseLine*strength/minV)
        return baseLine*strength/minV
    elif(strength < maxV):
        r = (strength - minV)/(maxV-minV)
        # print(np.mean(baseLine*(1-r) + arc*(r)))
        return baseLine*(1-r) + arc*(r)
    else:
        # print('max')
        # print(np.mean(arc*strength/maxV))
        return arc*strength/maxV

def meltFlatRate(strength,n):
    x = np.linspace(1,1,n)
    return x*strength

def meltLinearRate(strength,n):
    x = np.linspace(0,1,n)
    baseLine = .6 -.4*x 
    return baseLine*strength/np.mean(baseLine)

def meltLessLinearRate(strength,n):
    x = np.linspace(0,1,n)
    baseLine = .5 -.2*x 
    return baseLine*strength/np.mean(baseLine)

def meltBackLinearRate(strength,n):
    x = np.linspace(0,1,n)
    baseLine = .3 +.2*x 
    return baseLine*strength/np.mean(baseLine)

def meltURate(strength,n):
    x = np.linspace(0,1,n)
    arc = .5 + 50*(x-.5)**6
    maxV = np.mean(arc) - .015 #extra bit seems to help glaciome get right value for U shape
    return arc*strength/maxV


# basic parameters needed for setting up the model; later will modify this so that 
# the fjord geometry can be passed through
constant = constants()


## You define experiment here, just these 
def meltFun(strength,n):
    return meltLessLinearRate(strength,n)

targetLength = 5000
targetMelt = .40
calvingRate = 6000


# Make some nice plots

# # print(Bview)
# plt.subplot(211)
# plt.plot(timePlot,Bview,color='red')
# plt.xlabel('Time [years]')
# plt.ylabel('Average Melt Rate [m/day]')
# plt.subplot(212)
# plt.plot(meltFun(Bview[0],50),color='red')
# plt.xlabel('Distance []')
# plt.ylabel('Average Melt Rate [m/day]')
# plt.savefig('meltForcing.png',format='png',dpi=200)
# # plt.show()
# plt.close()


# files = sorted(glob.glob('../chkpt035LessLinear.pickle'))
files = sorted(glob.glob('../PACE/F6000/*.pickle'))
# files = sorted(glob.glob('../RecoverLinear_405/recover00100.pickle'))

for j in np.arange(0,len(files)):
    file = open(files[j], 'rb')
    data = pickle.load(file)
    file.close()
    if(abs(data.L - targetLength) < 20):
        break

print(data)

#Reset to time = 0, set dt
data.t = 0

lastX = data.L

alpha = 0e-5 #buttressing coefficient (25e-5 so far have been good) [m^2 yr^-1 N ^-1]
U0 = data.Uc + alpha*data.force() #initialize stable
print(f"\talpha {alpha:.2e}, U0 is {U0:3.2e} m/yr")

V = 0
V_old = 0

print("Done!")







