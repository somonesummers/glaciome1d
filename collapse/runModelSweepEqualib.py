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

minMelt = .35
maxMelt = .7
calvingRate = 6000
yearsToSimulate = 10 #saved as years, but you input days/365.25
dt = 10/365.25 # [years]

# Make some nice plots
endTime = int(yearsToSimulate/dt)
timePlot = np.arange(endTime)*dt
iList = np.arange(0,endTime,1) 
Bview = np.linspace(minMelt,maxMelt,endTime)
# print(Bview)
plt.subplot(211)
plt.plot(timePlot,Bview,color='red')
plt.xlabel('Time [years]')
plt.ylabel('Average Melt Rate [m/day]')
plt.subplot(212)
plt.plot(meltFun(Bview[0],50),color='red')
plt.xlabel('Distance []')
plt.ylabel('Average Melt Rate [m/day]')
plt.savefig('meltForcing.png',format='png',dpi=200)
# plt.show()
plt.close()


## Run to steady state
if(False):
    #initial mélange values
    n_pts = 21 # number of grid points
    L = 15e3 # ice melange length
    Ut = calvingRate # glacier terminus velocity [m/a]; 
    Uc = calvingRate # glacier calving rate [m/a]; 
    Ht = 600 # terminus thickness
    # n = 101 # number of time steps
    # specifying fjord geometry
    X_fjord = np.linspace(-200e3,200e3,101)
    Wt = 5600
    W_fjord = Wt + 0/10000*X_fjord
    B_const = -1*minMelt * constant.daysYear # will replace with grid later
#    data = glaciome(n_pts, dt, L, Ut, Uc, Ht, B_const, X_fjord, W_fjord)
    files = sorted(glob.glob('steadystate.pickle'))
    for j in np.arange(0,len(files)):
          file = open(files[j], 'rb')
          data = pickle.load(file)
          file.close()
    data.Uc = calvingRate
    data.X_externalGrid = data.X
    data.B_externalGrid = meltFun(B_const,n_pts)

    #prep loop variables
    L_old = data.X[-1]
    V_old = 0
    V = 1e7
    flag = 0
    t_step_old = time.time()
    k = 0 
    k_step = 5
    t_old = 0
    while(flag < 3):
        data.prognostic(method='lm')
        
        # calculate rate of volume change
        X_ = np.concatenate(([data.X[0]], data.X_, [data.X[-1]]))
        H = np.concatenate(([data.H0], data.H, [data.HL]))
        W = np.concatenate(([data.W0], data.W, [data.WL]))
        V = simpson(H*W, x=X_)*1e-9
    
        dVdt = (V-V_old)/data.dt
        if((1.5*data.H[-1]-0.5*data.H[-2]) < 24):
            print('\t WARNING LOW H_L: ' + "{:.2f}".format(1.5*data.H[-1]-0.5*data.H[-2]) + ' m')
        print('CFL: ' + "{:.4f}".format(np.max(data.U)*data.dt/(data.dx*data.L)))
        if (k%k_step) == 0:
            #Adjust timesteps
            data.dt = (data.dx*data.L)/np.max(data.U) * 3.0 #set this to what you want in CLF
            
            #Print
            t_step = time.time()
            print('Step: ' + str(int(k)) )
            print('\tAdjust dt to %.4g days' % (data.dt * constant.daysYear))
            print('\tTime per step: ' + "{:.2f}".format((t_step-t_step_old)/k_step) + ' s')
            print('\tSimulation time: ' + "{:.3f}".format(data.t) + ' yr')
            print('\tLength: ' + "{:.2f}".format(data.L) + ' m')
            print('\tdL/dt: ' + "{:.2f}".format((data.L-L_old)/(data.t-t_old)) + ' m/yr') # over last time step
            print('\tVolume: ' + "{:.4f}".format(V) + ' km^3')
            print('\tdV/dt: ' + "{:.2f}".format(dVdt) + ' km^3/yr' )
            print('\tH_L: ' + "{:.2f}".format(1.5*data.H[-1]-0.5*data.H[-2]) + ' m') 
            print('\tH_0: ' + "{:.2f}".format(data.H0) + ' m') 
            print('\tCFL: ' + "{:.4f}".format(np.max(data.U)*data.dt/(data.dx*data.L)))
            if( np.abs(data.L-L_old)/(data.t-t_old) > 10e3 or np.abs(dVdt) > 5):
                print('\t\t****** WARNING ******')
                print('\t\tVery far from equalib')
            tmpFileName = 'tempfile_%05i.pickle' % k
            print('Saving intermediate: %s' %tmpFileName)
            data.save(tmpFileName)        
            #check for convergence
            dLdt = (data.L-L_old)/(data.t-t_old)
            if np.abs(dLdt) < 10:
                flag += 1
            
        #update vars and meltrate
        meltGridScale = (2*data.L - L_old)/data.L
        data.X_externalGrid = np.linspace(0,1,50)*data.L * meltGridScale
        data.B_externalGrid = meltFun(B_const,50) 
        L_old = data.L
        t_step_old = t_step
        V_old = V
        t_old = data.t
        k += 1
    print(f'Converged! {data}')
    data.save('steadystate.pickle')
else: #direct load
     files = sorted(glob.glob('steadystate.pickle'))

     for j in np.arange(0,len(files)):
         file = open(files[j], 'rb')
         data = pickle.load(file)
         file.close()


#Reset to time = 0, set dt
data.t = 0
data.dt = dt

lastX = data.L

alpha = 0e-5 #buttressing coefficient (25e-5 so far have been good) [m^2 yr^-1 N ^-1]
U0 = data.Uc + alpha*data.force() #initialize stable
print(f"\talpha {alpha:.2e}, U0 is {U0:3.2e} m/yr")

V = 0
V_old = 0
dLdt = 100
for i in iList: 
    j = 0
    while(np.abs(dLdt) > 10):
        V = simpson(data.H, x=data.X_) #m^2
        data.dt = (data.dx*data.L)/np.max(data.U) * 3.0 #target CFL
        meltGridScale = (2*data.L - lastX)/data.L
        data.X_externalGrid = np.linspace(0,1,50)*data.L * meltGridScale
        data.B_externalGrid = -1*meltFun(Bview[i],50) * constant.daysYear
        lastX = data.L
        t_old = data.t
        data.prognostic(method='lm') # lm or hybr
        data.save(f'zTemp{j:05d}.pickle')
        j += 1
        V_old = V 
        if(data.L < 1000 or np.min(data.H) < 24.5):
            break
        dLdt = (data.L-lastX)/(data.t-t_old)
        print(f'dt {data.dt*365.25:4.2f} days, dLdt {dLdt:4.2f}, L:{data.L:7.0f} m, L:{data.L:7.0f} m, Uf {data.U[-1]/constant.daysYear:5.1f} m/day')
    if(data.L < 1000 or np.min(data.H) < 24.5):
        break
    #At this point, we're at equalib on stable node. Delete temp files, save point, reset for next loop
    os.system(f'rm -f zTemp*') 
    print(f"t {data.t*constant.daysYear:5.1f} day, index {i:05d}, with H0:{data.H0:7.2f} m, L:{data.L:7.0f} m, Uc:{data.Uc:5.0f} m/yr, Uf {data.U[-1]/constant.daysYear:5.1f} m/day, ∆Vol: {V-V_old:8.2g} m^2, melt {np.mean(data.B/constant.daysYear):6.4f} m/day ({Bview[i]*-1})")
    data.save(f'sweep{i:05d}.pickle')
    dLdt = 100 #kick it back into while loop
print("Done!")







