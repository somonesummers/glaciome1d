import numpy as np
from matplotlib import pyplot as plt
import os
import sys
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

n_pts = 21 # number of grid points
L = 15e3 # ice melange length
Ut = 0.6e4 # glacier terminus velocity [m/a]; treated as a constant
Uc = 0.6e4 # glacier calving rate [m/a]; treated as a constant
Ht = 600 # terminus thickness
# n = 101 # number of time steps
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

# files = sorted(glob.glob('../chkpt035LessLinear.pickle'))
files = sorted(glob.glob('chkpt035LL6100.pickle'))
# files = sorted(glob.glob('../RecoverLinear_405/recover00100.pickle'))

for j in np.arange(0,len(files)):
    file = open(files[j], 'rb')
    data = pickle.load(file)
    file.close()

# data.W_fjord = 5600 + 0*data.X_fjord

#Reset to time = 0
data.t = 0
dt = 10/365.25 # [years]
data.dt = dt

yearsToSimulate = 30 #saved as years, but you input days/365.25
endTime = int(yearsToSimulate/dt)
time = np.arange(endTime)*dt
iList = np.arange(0,endTime,1) 
if(data.t > 0):
    iList[:] = iList[:] + int(data.t/dt) + 1 
print(iList)


minMelt = .35
maxMelt = .65
def meltFun(strength,n):
    return meltLessLinearRate(strength,n)

Bview = np.linspace(minMelt,maxMelt,endTime)

# print(Bview)
plt.subplot(211)
plt.plot(time,Bview,color='red')
plt.xlabel('Time [years]')
plt.ylabel('Average Melt Rate [m/day]')
plt.subplot(212)
plt.plot(meltFun(Bview[0],50),color='red')
plt.xlabel('Distance []')
plt.ylabel('Average Melt Rate [m/day]')
plt.savefig('meltForcing.png',format='png',dpi=200)
plt.show()
plt.close()

lastX = data.X[-1]

data.Uc = 6100
alpha = 0e-5 #buttressing coefficient (25e-5 so far have been good) [m^2 yr^-1 N ^-1]
U0 = data.Uc + alpha*data.force() #initialize stable
print(f"\talpha {alpha:.2e}, U0 is {U0:3.2e} m/yr")

vol = 0
volOld = 0
for i in iList: 
    # F = data.force()
    data.Uc = U0 - alpha * data.force() 
    # data.Ht = 600 - data.X[0] * beta ## increase in thickness with retreat
    meltGridScale = (2*data.X[-1] - lastX)/data.X[-1]
    data.X_externalGrid = np.linspace(0,1,50)*data.X[-1] * meltGridScale
    vol = simpson(data.H, x=data.X_) #m^2
    # print(f'scale grid by {meltGridScale:3.2f}, guess Length {data.X_externalGrid[-1]:7.0f} m')
    # print(data.X)
    # print(-1*meltRate(Bview[i],n_pts) * constant.daysYear)

    data.B_externalGrid = -1*meltFun(Bview[i],50) * constant.daysYear
    lastX = data.X[-1]
    data.prognostic(method='lm') # lm or hybr
    if(i % 1 == 0):
        print(f"t {data.t*constant.daysYear:5.1f} day with H0:{data.H0:7.2f} m, L:{data.L:7.0f} m, Uc:{data.Uc:5.0f} m/yr, Uf {data.U[-1]/constant.daysYear:5.1f} m/day, ∆Vol: {vol-volOld:5.2g} m^2, melt {np.mean(data.B/constant.daysYear):0.5f} m/day")
        data.save(f'sweep{i:05d}.pickle')
    volOld = vol
    

print("Done!")







