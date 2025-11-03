#!/usr/bin/env python
# coding: utf-8
import sys
import numpy as np
from matplotlib import pyplot as plt
sys.path.append('/Users/psummers8/Documents/glaciome1D')
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
from glaciome1D import constants, glaciome
import scipy

import glob
import pickle

folders = ['W4000_L6000','W4800_L6000','W6800_L6000','zOLD/L3000_slow','zOLD/L4500_slow','zOLD/L6000_slow']
fileEnding = 'L'
# folders = ['W4000_L6000','W4800_L6000','L6000_slow','W6800_L6000','']
# fileEnding = 'widths'

files = folders.copy()

for i in range(len(folders)):
	fold = folders[i]
	fileTmp = sorted(glob.glob(f'PACE/{fold}/swee*.pickle'))
	if(len(fileTmp) >0):
		files[i] = fileTmp[-1]

nameList = folders


toIterate = np.arange(len(files))
H0Time = np.zeros(np.shape(toIterate))
UfTime = np.zeros(np.shape(toIterate))
lengthTime = np.zeros(np.shape(toIterate))
bTime =  np.zeros(np.shape(toIterate))
UcTime =  np.zeros(np.shape(toIterate))
wTime = np.zeros(np.shape(toIterate))
for j in toIterate:
    file = files[j]
    with open(files[j], 'rb') as fileName:
        data = pickle.load(fileName)
        fileName.close()
    H0Time[j]=data.H0
    lengthTime[j]=data.X[-1] - data.X[0]
    UfTime[j]=data.U[-1]
    bTime[j] = -1 * np.mean(data.B) /365.25
    UcTime[j] = data.Uc
    wTime[j] = np.mean(data.W_fjord)


fig, axes = plt.subplots(2, 2, figsize=(12, 6), layout="constrained")
ax1 = axes[0,0]
ax2 = axes[0,1]
ax3 = axes[1,0]
ax4 = axes[1,1]

sc = ax1.scatter(UcTime,lengthTime/1e3,bTime*500,wTime/1e3,marker='o',cmap='viridis',alpha=.8)
ax1.set_axisbelow(True)
# ax1.set_facecolor('xkcd:light grey')
ax1.grid(color='gray',alpha=.2)
cbar = plt.colorbar(sc)
cbar.set_label('Fjord Width [km]')
ax1.set_ylabel('Minimum Length [km]')
ax1.set_xlabel('Cavling Rate [m/y]')
ax1.scatter([],[],[.3*500],label='0.3 m/d', color='gray')
ax1.scatter([],[],[.5*500],label='0.5 m/d', color='gray')
ax1.legend( title="Avg Melt Rate")


sc = ax2.scatter(wTime/1e3,lengthTime/1e3,bTime*500,UcTime,marker='o',cmap='cividis',alpha=.8)
# ax1.legend()
cbar = plt.colorbar(sc)
ax2.set_axisbelow(True)
ax2.grid(color='gray',alpha=.2)
cbar.set_label('Calving Speed [m/y]')
ax2.set_ylabel('Minimum Length [km]')
ax2.set_xlabel('Fjord Width [km]')
ax2.scatter([],[],[.3*500],label='0.3 m/d', color='gray')
ax2.scatter([],[],[.5*500],label='0.5 m/d', color='gray')
ax2.legend( title="Avg Avg Melt Rate")

sc = ax3.scatter(bTime,lengthTime/1e3,(UcTime)/20,wTime/1e3,marker='o',cmap='viridis',alpha=.8)
# ax1.set_alpha(.25)
ax3.set_axisbelow(True)
# ax1.set_facecolor('xkcd:light grey')
ax3.grid(color='gray',alpha=.2)
cbar = plt.colorbar(sc)
cbar.set_label('Fjord Width [km]')
ax3.set_ylabel('L min [km]')
ax3.set_xlabel('Maximum Avg Melt Rate [m/d]')
ax3.scatter([],[],[3000/20],label='3000 m/y', color='gray')
ax3.scatter([],[],[4500/20],label='4500 m/y', color='gray')
ax3.scatter([],[],[6000/20],label='6000 m/y', color='gray')
ax3.legend( title="Calving Speed ")


sc = ax4.scatter(bTime,lengthTime/1e3,wTime/20,UcTime,marker='o',cmap='cividis',alpha=.8)
# ax1.set_alpha(.25)
ax4.set_axisbelow(True)
# ax1.set_facecolor('xkcd:light grey')
ax4.grid(color='gray',alpha=.2)
cbar = plt.colorbar(sc)
cbar.set_label('Calving Speed [m/y]')
ax4.set_ylabel('L min [km]')
ax4.set_xlabel('Maximum Avg Melt Rate [m/d]')
ax4.scatter([],[],[4000/20],label='4 km', color='gray')
ax4.scatter([],[],[6800/20],label='6.8 km', color='gray')
ax4.legend( title="Fjord Width")


plt.savefig(f"figs/Fig3_{fileEnding}.png", format='png', dpi=400)
plt.show()
