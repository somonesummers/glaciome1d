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
    return meltLinearRate(strength,n)


n = 10 #how long to shrink to confirm behavior 
m = 1 #how many extra collapses are needed to confirm unstable fix point

targetLength = 5000
meltToTry = [.425, .40, .38, .35, .33, .31] #L6000 collapse at  0.4283, L 10905
# H0:  41.17 m, L:  10905 m, Uc: 6000 m/yr, Uf 213.4 m/day, ∆Vol:        0 m^2, melt -0.4250 m/day (-0.425)
lengthsToTry = [10000,9000,8000,7000,6500,6000,5520,5100,4600,4300,3900]
print(meltToTry)
print(lengthsToTry)
for targetMelt in meltToTry:
    collapseCount = 0
    for targetLength in lengthsToTry:
        # files = sorted(glob.glob('../PACE/F6000/*.pickle'))
        files = sorted(glob.glob('BaseFiles/*.pickle'))
        # files = sorted(glob.glob('../RecoverLinear_405/recover00100.pickle'))
        found = False
        minMiss = [1000,0]
        for j in np.arange(0,len(files)):
            file = open(files[j], 'rb')
            data = pickle.load(file)
            file.close()
            err = abs(data.L - targetLength)
            if(err < minMiss[0]):
                minMiss[0] = err
                minMiss[1] = data.L
            if(err < 20):
                found = True
                break

        if( not found):
            raise Exception(f"no matching length found for L = {targetLength} m, closest is {minMiss[1]:.0f}")

        # print(data)

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
        # print(f"\tt {data.t*constant.daysYear:6.1f} d, index {i:05d}, H0{data.H0:7.2f} m, L{data.L:7.0f} m, dLdt{dLdt:10.2f} m/y, dHdt{dHdt:7.2f} m/y, Uc{data.Uc:7.0f} m/y, Uf{data.U[-1]/constant.daysYear:6.1f} m/d, ∆Vol{V-V_old:9.2g} m^2, melt {np.mean(data.B/constant.daysYear):6.4f} m/d ({targetMelt*-1})")
        data.save(f'output_{targetLength}_{i:05d}.pickle')
        i += 1
        while(negCounter < n and posCounter < n):
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
            # print(f"\tt {data.t*constant.daysYear:6.1f} d, index {i:05d}, H0{data.H0:7.2f} m, L{data.L:7.0f} m, dLdt{dLdt:10.2f} m/y, dHdt{dHdt:7.2f} m/y, Uc{data.Uc:7.0f} m/y, Uf{data.U[-1]/constant.daysYear:6.1f} m/d, ∆Vol{V-V_old:9.2g} m^2, melt {np.mean(data.B/constant.daysYear):6.4f} m/d ({targetMelt*-1})")
            data.save(f'output_{targetLength}_{i:05d}.pickle')
            if(dHdt < 0 and dLdt < 0):
                negCounter += 1
            elif(dHdt > 0 and dLdt > 0):
                posCounter += 1
            i += 1
        if(negCounter == n):
            collapseCount += 1
        print(f"  L {targetLength:7.0f}, melt {targetMelt:7.4f}: neg/pos count {negCounter}/{posCounter}")
        if(collapseCount > m):
            break #3 collapses for a given melt rate is good enough







