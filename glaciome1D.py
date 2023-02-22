#### WARNING ####
# Because of how the global variables are currently imported, the interpreter
# needs to be re-started if any changes are made to the them.

import numpy as np
from matplotlib import pyplot as plt
import matplotlib

from scipy.optimize import root
from scipy.optimize import minimize
from scipy.optimize import fsolve

from general_utilities import second_invariant, width

import config

import importlib

# here you should specify the rheology that should be imported
rheology = 'granular_fluidity'
model = importlib.import_module(rheology)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: jason
"""


#%

dt = config.secsDay*30 # time step [s]
n = 21 # number of time steps

B = -0/config.secsYear # mass balance rate [m s^-1]

x0 = 0 # left boundary of the melange [m]
L = 10000 # initial ice melange length [m]
dx = 0.1  # grid spacing [m]
x = np.arange(0,1+dx,dx) # longitudinal grid
X = x*L # coordinates of unstretched grid


W = width((X[:-1]+X[1:])/2)
#W = 4000*np.ones(len(x)-1) # fjord width [m]; treated as constant for now
# W needs to move with the grid...

H = np.ones(len(x)-1)*config.d # initial ice melange thickness [m]


Ut = 10000/config.secsYear # width-averaged glacier terminus velocity [m/s]
U = Ut*(1-x) # initial guess for the averaged velocity [m/s]; the model
# determine velocity profile that is consistent with initial thickness; unlike 
# subsequent steps this does not involve an implicit time step

#U = fsolve(model.spinup, U, (x,X,L,Ut,H,W,dx,dt))

#%%


fig_width = 12
fig_height = 6.5
plt.figure(figsize=(fig_width,fig_height))

ax_width = 3/fig_width
ax_height = 2/fig_height
left = 1/fig_width
bot = 0.5/fig_height
ygap = 0.75/fig_height
xgap= 1/fig_width

ax1 = plt.axes([left, bot+ax_height+2.25*ygap, ax_width, ax_height])
ax1.set_xlabel('Longitudinal coordinate [m]')
ax1.set_ylabel('Speed [m/d]')
ax1.set_ylim([0,40])
ax1.set_xlim([0,10000])

ax2 = plt.axes([left+ax_width+xgap, bot+ax_height+2.25*ygap, ax_width, ax_height])
ax2.set_xlabel('Longitudinal coordinate [m]')
ax2.set_ylabel('Thickness [m]')
ax2.set_ylim([-75, 100])
ax2.set_xlim([0,10000])

ax3 = plt.axes([left, bot+1.25*ygap, ax_width, ax_height])
ax3.set_xlabel('Longitudinal coordinate [m]')
ax3.set_ylabel('$\mu$')
ax3.set_ylim([0, 0.3])
ax3.set_xlim([0,10000])

ax4 = plt.axes([left+ax_width+xgap, bot+1.25*ygap, ax_width, ax_height])
ax4.set_xlabel('Longitudinal coordinate [m]')
ax4.set_ylabel('$\mu_w$')
ax4.set_ylim([0, 1])
ax4.set_xlim([0,10000])

ax5 = plt.axes([left+2*(ax_width+xgap), bot+1.25*ygap, 0.75*ax_width, 2*ax_height+ygap])
ax5.set_xlabel('Transverse coordinate [m]')
ax5.set_ylabel(r'Speed at $\chi=0.2$ [m/d]')
ax5.set_xlim([0,4000])
ax5.set_ylim([0,40])

ax_cbar = plt.axes([left, bot, 2*(ax_width+xgap)+0.75*ax_width, ax_height/15])

cbar_ticks = np.linspace(0, (n-1)*dt/config.secsDay, 11, endpoint=True)
cmap = matplotlib.cm.viridis
bounds = cbar_ticks
norm = matplotlib.colors.Normalize(vmin=0, vmax=(n-1)*dt/config.secsDay)
cb = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap=cmap, norm=norm,
                                orientation='horizontal')#,extend='min')
cb.set_label("Time [d]")

#%%


color_id = np.linspace(0,1,n) 


# plot initial time step
#ax1.plot(X,U*config.secsDay,color=plt.cm.viridis(color_id[0]))
#ax2.plot((X[:-1]+X[1:])/2,H,color=plt.cm.viridis(color_id[0]))

#mu, muW = model.get_mu(x,U,H,W,X[-1]-X[0],dx)

#ax3.plot((X[:-1]+X[1:])/2,mu,color=plt.cm.viridis(color_id[0]))
#ax4.plot(X[1:-1],muW,color=plt.cm.viridis(color_id[0]))

# concatenate U and H since the implicit time step requires that they are
# iteratively solved simultaneously
UH = np.append(U,H)
#UHL = np.concatenate((U,H,[L]))

for k in np.arange(1,n):
    print('Time: ' + "{:.0f}".format(k*dt/config.secsDay) + ' days')     

    UH = fsolve(model.convergence, UH, (x,X,Ut,H,W,dx,dt,U,H,B))

    
    # Note: I think all that needs to be saved are UH, W, and the initial fjord length
    
    # the following is for plotting purposes
    U = UH[:len(x)]
    H = UH[len(x):]
    #L = UHL[-1]
    
    xt = X[0] + U[0]*dt
    xL = X[-1] + U[-1]*dt
    X = np.linspace(xt,xL,len(x))-xt
    X_ = (X[:-1]+X[1:])/2
    
    mu, muW = model.get_mu(x,U,H,W,X[-1]-X[0],dx)
    
    W = width((X[:-1]+X[1:])/2)
    
    ind = int(len(W)/5) # index of midpoint in fjord
    y, u_transverse, _ = model.transverse(W[ind],muW[ind],H[ind])
    
    
    ax1.plot(X,U*config.secsDay,color=plt.cm.viridis(color_id[k]))
    
    
    ax2.plot(np.append(X_,X_[::-1]),np.append(-config.rho/config.rho_w*H,(1-config.rho/config.rho_w)*H[::-1]),color=plt.cm.viridis(color_id[k]))
    
    ax3.plot(X_,mu,color=plt.cm.viridis(color_id[k]))
    ax4.plot(X[1:-1],muW,color=plt.cm.viridis(color_id[k]))
    ax5.plot(np.append(y,y+y[-1]),np.append(u_transverse,u_transverse[-1::-1])*config.secsDay,color=plt.cm.viridis(color_id[k]))

ax3.plot(np.array([0,1e4]),np.array([config.muS,config.muS]),'k:')
ax4.plot(np.array([0,1e4]),np.array([config.muS,config.muS]),'k:')
