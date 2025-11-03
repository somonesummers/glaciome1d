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


folders = ['L3000_slow','L4500_slow','L6000_slow','','']
fileEnding = 'widths'
UcList = np.array([3000,4500,6000])


files = sorted(glob.glob(f'PACE/{folders[0]}/swee*.pickle'))
ind1 = len(files)
files += sorted(glob.glob(f'PACE/{folders[1]}/swee*.pickle'))
ind2 = len(files)
files += sorted(glob.glob(f'PACE/{folders[2]}/swee*.pickle'))
ind3 = len(files)
files += sorted(glob.glob(f'PACE/{folders[3]}/sw*.pickle'))
ind4 = len(files)
# files += sorted(glob.glob(f'PACE/{folders[4]}/sw*.pickle'))
ind5 = len(files)
print(ind1,ind2,ind3,ind4,ind5)
colorList = ['xkcd:blue','xkcd:green','xkcd:red','xkcd:indigo','xkcd:lavender']
nameList = folders

fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout="constrained")
ax1 = axes[0]
ax2 = axes[1]

# files += sorted(glob.glob('sweepTemp/*.pickle'))
# print(f'index 1,2: {ind1},{ind2}')
toIterate = np.arange(len(files))
H0Time = np.zeros(np.shape(toIterate))
UfTime = np.zeros(np.shape(toIterate))
lengthTime = np.zeros(np.shape(toIterate))
bTime =  np.zeros(np.shape(toIterate))
for j in toIterate:
    file = files[j]
    with open(files[j], 'rb') as fileName:
        data = pickle.load(fileName)
        fileName.close()
    H0Time[j]=data.H0
    lengthTime[j]=data.X[-1] - data.X[0]
    UfTime[j]=data.U[-1]
    bTime[j] = -1 * np.mean(data.B) /365.25

Ht = 600 #m
Hf = 25 #m


uf = np.linspace(50, 400, 100)
l = np.linspace(1, 25000, 100)
Uf, L = np.meshgrid(uf,l)


# plt.figure(3,figsize=(8, 5))
# # plt.plot(l/1e3,fluxPlug/Hf/365.25,color='gray')
# plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]/365.25,[2],alpha=1,color=colorList[0])
# fun1 = scipy.interpolate.interp1d(lengthTime[:ind1],UfTime[:ind1]/365.25,bounds_error=False,fill_value=np.nan)
# plt.plot(l/1e3,fun1(l),color=colorList[0],linestyle='-',alpha=.2,label=nameList[0])
loadData1 = np.load(f'PACE/{folders[0]}/unstableNodes.npy')

# plt.scatter(lengthTime[ind1:ind2]/1e3,UfTime[ind1:ind2]/365.25,[2],alpha=1,color=colorList[1])
# fun2 = scipy.interpolate.interp1d(lengthTime[ind1:ind2],UfTime[ind1:ind2]/365.25,bounds_error=False,fill_value=np.nan)
# plt.plot(l/1e3,fun2(l),color=colorList[1],linestyle='-',alpha=.2,label=nameList[1])
loadData2 = np.load(f'PACE/{folders[1]}/unstableNodes.npy')

if(ind3>ind2):
	a = ind2
	b = ind3
	# plt.scatter(lengthTime[a:b]/1e3,UfTime[a:b]/365.25,[2],alpha=1,color=colorList[2])
	# fun3 = scipy.interpolate.interp1d(lengthTime[a:b],UfTime[a:b]/365.25,bounds_error=False,fill_value=np.nan)
	# plt.plot(l/1e3,fun3(l),color=colorList[2],linestyle='-',alpha=.2,label=nameList[2])
	loadData3 = np.load(f'PACE/{folders[2]}/unstableNodes.npy')

if(ind4>ind3):
	a = ind3
	b = ind4
	# plt.scatter(lengthTime[a:b]/1e3,UfTime[a:b]/365.25,[2],alpha=1,color=colorList[3])
	# fun4 = scipy.interpolate.interp1d(lengthTime[a:b],UfTime[a:b]/365.25,bounds_error=False,fill_value=np.nan)
	# plt.plot(l/1e3,fun4(l),color=colorList[3],linestyle='-',alpha=.2,label=nameList[3])
	loadData4 = np.load(f'PACE/{folders[3]}/unstableNodes.npy')

if(ind5>ind4):
	a = ind4
	b = ind5
	# plt.scatter(lengthTime[a:b]/1e3,UfTime[a:b]/365.25,[2],alpha=1,color=colorList[4])
	# fun5 = scipy.interpolate.interp1d(lengthTime[a:b],UfTime[a:b]/365.25,bounds_error=False,fill_value=np.nan)
	# plt.plot(l/1e3,fun5(l),color=colorList[4],linestyle='-',alpha=.2,label=nameList[4])
	loadData5 = np.load(f'PACE/{folders[4]}/unstableNodes.npy')




for b in np.asarray([bTime[-1]+.1,bTime[-1],bTime[-1]-.1])*365.25:
	ax1.plot(l/1e3,UcList[2]*Ht - b*l,color='gray',alpha=.2,linewidth=1)#label=f'Flux in - Melt (b = {b/365.25:.2f})')

ax1.scatter(lengthTime[ind2:]/1e3,UfTime[ind2:]*Hf,[25],alpha=.8,color=colorList[2])

lengthCutOff = 3e3
lengthFromUnstable, unq_ind = np.unique((loadData3[loadData3[:,1] > lengthCutOff,1]+loadData3[loadData3[:,1] > lengthCutOff,2])/2,return_index=True)

## Parabolic Fit for fun
# l_toFit = np.append(lengthFromUnstable,lengthTime[ind2:])
# f_toFit = np.append((UcList[2]*Ht - loadData3[unq_ind,0]*365.25*lengthFromUnstable),UfTime[ind2:]*Hf)
# g = np.ones([len(l_toFit),3])
# g[:,1] = l_toFit/1e3
# g[:,2] = (l_toFit/1e3)**2
# m = np.linalg.pinv(g)@f_toFit
# ax1.plot(l/1e3,m[0] + m[1]*l/1e3 + m[2]*(l/1e3)**2,color='gray',alpha=.5,linewidth=1,linestyle='--')

ax1.scatter(lengthFromUnstable/1e3,(UcList[2]*Ht - loadData3[unq_ind,0]*365.25*lengthFromUnstable),[25],alpha=.8,edgecolors=colorList[2],color='white')
ax1.set_ylim([0,4e6])
ax1.set_title('Flux Diagram')
ax1.set_xlabel('Length [km]')
ax1.set_ylabel('Volume Flux [m^2/year]')
# ax1.legend()

lengthFromUnstable, unq_ind = np.unique((loadData1[loadData1[:,1] > lengthCutOff,1]+loadData1[loadData1[:,1] > lengthCutOff,2])/2,return_index=True)
ax2.plot(np.append(loadData1[unq_ind,0],bTime[ind1-1]),np.append((loadData1[unq_ind,1]+loadData1[unq_ind,2])/2.0,
	lengthTime[ind1-1]),color=colorList[0],linestyle='--')
ax2.plot(bTime[:ind1],lengthTime[:ind1],color=colorList[0],linestyle='-',label=nameList[0])

lengthFromUnstable, unq_ind = np.unique((loadData2[loadData2[:,1] > lengthCutOff,1]+loadData2[loadData2[:,1] > lengthCutOff,2])/2,return_index=True)
ax2.plot(np.append(loadData2[unq_ind,0],bTime[ind2-1]),np.append((loadData2[unq_ind,1]+loadData2[unq_ind,2])/2.0,
	lengthTime[ind2-1]),color=colorList[1],linestyle='--')
ax2.plot(bTime[ind1:ind2],lengthTime[ind1:ind2],color=colorList[1],linestyle='-',label=nameList[1])

if(ind3>ind2):
	ax2.plot(np.append(loadData3[unq_ind,0],bTime[ind3-1]),np.append((loadData3[unq_ind,1]+loadData2[unq_ind,2])/2.0,
	lengthTime[ind3-1]),color=colorList[2],linestyle='--')
ax2.plot(bTime[ind2:ind3],lengthTime[ind2:ind3],color=colorList[2],linestyle='-',label=nameList[2])
if(ind4>ind3):
	ax2.plot(np.insert(loadData4[:,0],0,bTime[ind4-1]),np.insert((loadData4[:,1]+loadData4[:,2])/2.0,0,
	lengthTime[ind4-1]),color=colorList[3],linestyle='--')
ax2.plot(bTime[ind3:ind4],lengthTime[ind3:ind4],color=colorList[3],linestyle='-',label=nameList[3])
if(ind5>ind4):
	ax2.plot(np.insert(loadData5[:,0],0,bTime[ind5-1]),np.insert((loadData5[:,1]+loadData5[:,2])/2.0,0,
	lengthTime[ind5-1]),color=colorList[4],linestyle='--')
ax2.plot(bTime[ind4:ind5],lengthTime[ind4:ind5],color=colorList[4],linestyle='-',label=nameList[4])

ax2.legend()
ax2.set_title('Stability Diagram')
ax2.set_ylabel('Length [km]')
ax2.set_xlabel('Melt Rate [m/day]')
plt.savefig(f"figs/fig2{fileEnding}.png", format='png', dpi=400)
plt.show()
