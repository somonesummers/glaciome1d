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
import warnings
import argparse


# from localVars import *
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: 
    Paul Summers
    August 2025
    Example of running glaciome1d, now with umelt shape

"""
parser = argparse.ArgumentParser(description='Run model to find unstable fix points')
parser.add_argument('-s','--shape', nargs=1, default=[0],type=int,
                    help='shape of melt profile [defaut = 0] 0 flat, 1 linear, 2 less linear, 3 U')
parser.add_argument('-v','--verbose', action='count', default=0,
                    help='how verbose to be')
parser.add_argument('-f','--fileSave', action='count', default=0,
                    help='how many details to save')
args = parser.parse_args()

print(f'input args: {args}')

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


## You define experiment in args input
if(args.shape[0]==0):
    def meltFun(strength,n):
        return meltFlatRate(strength,n)
elif(args.shape[0]==1):
    def meltFun(strength,n):
        return meltLinearRate(strength,n)
elif(args.shape[0]==2):
    def meltFun(strength,n):
        return meltLessLinearRate(strength,n)
elif(args.shape[0]==3):
    def meltFun(strength,n):
        return meltURate(strength,n)


verboseLevel = args.verbose

n_up = 5 #how long to rise to confirm behavior 
n_down = 10 #how long to shrink to confirm
m = 1 #how many extra collapses are needed to confirm unstable fix point
lookBack = 3 #how much longer than last collapse do we start check
libraryDirectory = 'zTemp'

#Load library
files = sorted(glob.glob(f'{libraryDirectory}*.pickle'))
lengthDict = np.zeros(len(files))
meltDict = np.zeros(len(files))
for i in np.arange(0,len(files)):
    file = open(files[i], 'rb')
    data = pickle.load(file)
    lengthDict[i] = data.L
    meltDict[i] = np.mean(data.B)*-1/365.25
    file.close()
initMelt = np.mean(meltDict) - .006 #initaly take a reasonable step away

meltToTry = initMelt - np.linspace(0,.2,100)

print(f'melt to try {meltToTry[:]}')
collapseHigh = np.zeros(len(meltToTry)) #last stable length
collapseLow = np.zeros(len(meltToTry)) #first unstable length
libraryIndex = 0
for j in range(len(meltToTry)):
    targetMelt = meltToTry[j]
    collapseCount = 0
    growCount = 0
    libraryIndex = np.max([libraryIndex - lookBack,0]) - 1
    while collapseCount < (m+1):
        libraryIndex += 1 #look at next longer melange from library
        if(libraryIndex >= len(lengthDict)):
            warnings.warn("Warning: Exceeded length library range")
            print("** WARN ** Exceeded length library range ** WARN **")
            break
        file = open(files[libraryIndex], 'rb')
        data = pickle.load(file)
        if(verboseLevel > 0):
            print(f'load index {libraryIndex:05d}, H0{data.H0:7.2f} m, L{data.L:7.0f} m')
        #Reset to time = 0, set dt
        data.t = 0
        lastX = data.L
        lastH = data.H0
        V = simpson(data.H, x=data.X_) #m^2
        V_old = 0
        dLdt = 100
        dHdt = 100
        i = 0
        negCounter = 0
        posCounter = 0
        if(verboseLevel > 1):
            print(f"\tt {data.t*constant.daysYear:6.1f} d, index {i:05d}, H0{data.H0:7.2f} m, L{data.L:7.0f} m, dLdt{dLdt:10.2f} m/y, dHdt{dHdt:7.2f} m/y, Uc{data.Uc:7.0f} m/y, Uf{data.U[-1]/constant.daysYear:6.1f} m/d, ∆Vol{V-V_old:9.2g} m^2, melt {np.mean(data.B/constant.daysYear):6.4f} m/d ({targetMelt*-1})")
        if(args.fileSave > 0):
            data.save(f'output_{data.L:05.0f}_{i:05d}.pickle')
        i += 1
        while(negCounter < n_down and posCounter < n_up):
        # while(np.abs(dLdt) > 10 and data.L > 3000):
            data.dt = (data.dx*data.L)/np.max(data.U) * 1.0 #target CFL
            meltGridScale = (2*data.L - lastX)/data.L
            data.X_externalGrid = np.linspace(0,1,50)*data.L * meltGridScale
            data.B_externalGrid = -1*meltFun(targetMelt,50) * constant.daysYear
            lastX = data.L
            lastH = data.H0
            t_old = data.t
            V_old = V
            data.prognostic(method='lm') # lm or hybr
            V = simpson(data.H, x=data.X_) #m^2
            dLdt = (data.L-lastX)/(data.t-t_old)
            dHdt = (data.H0-lastH)/(data.t-t_old)
            if(verboseLevel > 1):
                print(f"\tt {data.t*constant.daysYear:6.1f} d, index {i:05d}, H0{data.H0:7.2f} m, L{data.L:7.0f} m, dLdt{dLdt:10.2f} m/y, dHdt{dHdt:7.2f} m/y, Uc{data.Uc:7.0f} m/y, Uf{data.U[-1]/constant.daysYear:6.1f} m/d, ∆Vol{V-V_old:9.2g} m^2, melt {np.mean(data.B/constant.daysYear):6.4f} m/d ({targetMelt*-1})")
            if(args.fileSave > 0):
                data.save(f'output_{data.L:05.0f}_{i:05d}.pickle')
            if(dHdt < 0 and dLdt < 0):
                negCounter += 1
            elif(dHdt > 0 and dLdt > 0):
                posCounter += 1
            i += 1
        if(verboseLevel > 0):
            print(f"  L {data.L:5.0f}, melt {targetMelt:7.4f}: neg/pos count {negCounter}/{posCounter}")
        if(negCounter == n_down):
            collapseCount += 1
        else:
            growCount += 1
    if(growCount < 1):
        warnings.warn("Warning: Failed to find positive side of unstable fix point")
        print("** WARN ** Failed to find positive side of unstable fix point ** WARN **")
    collapseHigh[j] = lengthDict[libraryIndex - (m+2)]
    collapseLow[j] = lengthDict[libraryIndex - (m+1)]
    print(f'unstable fix point: {meltToTry[j]:7.4f}, {collapseLow[j]:7.4f}, {collapseHigh[j]:7.4f}')

unstableNodes = np.zeros([len(meltToTry),3])
for j in range(len(meltToTry)):
    print(f'{meltToTry[j]:7.4f}: {collapseLow[j]:7.4f}, {collapseHigh[j]:7.4f}')
    unstableNodes[j,0] = meltToTry[j]
    unstableNodes[j,1] = collapseLow[j]
    unstableNodes[j,2] = collapseHigh[j]
np.save('unstableNodes',unstableNodes)








