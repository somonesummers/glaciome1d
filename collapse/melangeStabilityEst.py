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


folders = ['W4000_L6000','W4800_L6000','L6000_slow','W6800_L6000','']
fileEnding = 'widths'
UcList = np.array([3000,4500,6000,4500,4500])


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
colorList = ['xkcd:red','xkcd:green','xkcd:blue','xkcd:indigo','xkcd:lavender']
nameList = folders



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

Ut = 4500 #m/yr
Ht = 600 #m
Hf = 25 #m
bMelt  = 0.25*365.25 #m/yr (quoted in m/day)

uf = np.linspace(50, 400, 100)
l = np.linspace(1, 25000, 100)
Uf, L = np.meshgrid(uf,l)

output = Ht*Ut - Uf*365.25*Hf - bMelt*L
fluxPlug = Ht*Ut/(l/230)*Hf

pwrShift = 50e3
fitLineL = 25

plt.figure(3,figsize=(8, 5))
# plt.plot(l/1e3,fluxPlug/Hf/365.25,color='gray')
plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]/365.25,[2],alpha=1,color=colorList[0])
fun1 = scipy.interpolate.interp1d(lengthTime[:ind1],UfTime[:ind1]/365.25,bounds_error=False,fill_value=np.nan)
plt.plot(l/1e3,fun1(l),color=colorList[0],linestyle='-',alpha=.2,label=nameList[0])
loadData1 = np.load(f'PACE/{folders[0]}/unstableNodes.npy')

plt.scatter(lengthTime[ind1:ind2]/1e3,UfTime[ind1:ind2]/365.25,[2],alpha=1,color=colorList[1])
fun2 = scipy.interpolate.interp1d(lengthTime[ind1:ind2],UfTime[ind1:ind2]/365.25,bounds_error=False,fill_value=np.nan)
plt.plot(l/1e3,fun2(l),color=colorList[1],linestyle='-',alpha=.2,label=nameList[1])
loadData2 = np.load(f'PACE/{folders[1]}/unstableNodes.npy')

if(ind3>ind2):
	a = ind2
	b = ind3
	plt.scatter(lengthTime[a:b]/1e3,UfTime[a:b]/365.25,[2],alpha=1,color=colorList[2])
	fun3 = scipy.interpolate.interp1d(lengthTime[a:b],UfTime[a:b]/365.25,bounds_error=False,fill_value=np.nan)
	plt.plot(l/1e3,fun3(l),color=colorList[2],linestyle='-',alpha=.2,label=nameList[2])
	loadData3 = np.load(f'PACE/{folders[2]}/unstableNodes.npy')

if(ind4>ind3):
	a = ind3
	b = ind4
	plt.scatter(lengthTime[a:b]/1e3,UfTime[a:b]/365.25,[2],alpha=1,color=colorList[3])
	fun4 = scipy.interpolate.interp1d(lengthTime[a:b],UfTime[a:b]/365.25,bounds_error=False,fill_value=np.nan)
	plt.plot(l/1e3,fun4(l),color=colorList[3],linestyle='-',alpha=.2,label=nameList[3])
	loadData4 = np.load(f'PACE/{folders[3]}/unstableNodes.npy')

if(ind5>ind4):
	a = ind4
	b = ind5
	plt.scatter(lengthTime[a:b]/1e3,UfTime[a:b]/365.25,[2],alpha=1,color=colorList[4])
	fun5 = scipy.interpolate.interp1d(lengthTime[a:b],UfTime[a:b]/365.25,bounds_error=False,fill_value=np.nan)
	plt.plot(l/1e3,fun5(l),color=colorList[4],linestyle='-',alpha=.2,label=nameList[4])
	loadData5 = np.load(f'PACE/{folders[4]}/unstableNodes.npy')

for i in range(len(UcList)):
	plt.plot(0*UcList[i],Ht*UcList[i]/Hf/365.25,marker='*',color=colorList[i],linewidth=0)
# print(f' fit is ax^p a: {m1[0]}, p: {m1[1]}')
# print(f' fit is ax^p a: {m2[0]}, p: {m2[1]}')
# print(f' fit is ax^p a: {m3[0]}, p: {m3[1]}')
# print(f' fit is ax^p a: {m4[0]}, p: {m4[1]}')
plt.title('Fitting')
plt.xlabel('Length [km]')
plt.ylabel('End Speed [m/day]')
plt.ylim([50,np.max([300,np.max(UcList)*Ht/Hf/365.25+20])])
plt.legend()
plt.savefig(f"figs/Fitting_{fileEnding}.png", format='png', dpi=400)
# plt.show()

plt.figure(1,figsize=(8, 5))
# plt.subplot(121)
cp = plt.contourf(L/1e3,Uf,output,np.linspace(-4e6,4e6,127),cmap='RdBu')
cc = plt.contour(L/1e3,Uf,output,[0],linestyles='--',colors='black')
plt.plot(l/1e3,fun1(l),color=colorList[0],linestyle='--',label=nameList[0])
plt.plot(l/1e3,fun2(l),color=colorList[1],linestyle='--',label=nameList[1])
plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]/365.25,[2],alpha=.25,color=colorList[0])
plt.scatter(lengthTime[ind1:ind2]/1e3,UfTime[ind1:ind2]/365.25,[2],alpha=.25,color=colorList[1])
if(ind3>ind2):
	plt.plot(l/1e3,fun3(l),color=colorList[2],linestyle='--',label=nameList[2])
	plt.scatter(lengthTime[ind2:ind3]/1e3,UfTime[ind2:ind3]/365.25,[2],alpha=.25,color=colorList[2])
if(ind4>ind3):
	plt.plot(l/1e3,fun4(l),color=colorList[3],linestyle='--',label=nameList[3])
	plt.scatter(lengthTime[ind3:ind4]/1e3,UfTime[ind3:ind4]/365.25,[2],alpha=.25,color=colorList[3])
if(ind5>ind4):
	plt.plot(l/1e3,fun5(l),color=colorList[4],linestyle='--',label=nameList[4])
	plt.scatter(lengthTime[ind4:ind5]/1e3,UfTime[ind4:ind5]/365.25,[2],alpha=.25,color=colorList[4])
plt.clabel(cc, inline=3, fontsize=8)
cbar = plt.colorbar(cp)
cbar.set_label('$\\partial V / \\partial t$')
plt.xlabel('Length [km]')
plt.ylabel('End Speed [m/day]')
plt.title(f'Change in Volume [m^2/yr] (melt = {bMelt/365.25})')
plt.ylim([np.min(uf),np.max(uf)])
plt.xlim([np.min(l)*1e-3,np.max(l)*1e-3])
plt.legend()
plt.tight_layout()
plt.savefig(f"figs/dVdt_{fileEnding}.png", format='png', dpi=400)

# plt.show()


plt.figure(2,figsize=(8, 5))
# plt.plot(l/1e3,fluxPlug,color='gray',label='Previous estimate of discharge')
for i in range(len(UcList)):
	if(i == 0 or UcList[i] != UcList[0]):
		for b in np.linspace(.2,.6,15)*365.25:
			plt.plot(l/1e3,UcList[i]*Ht - b*l,color=colorList[i],alpha=.2,linewidth=.5)#label=f'Flux in - Melt (b = {b/365.25:.2f})')


plt.plot(l/1e3,fun1(l)*365.25*Hf,color=colorList[0],linestyle='--',label=nameList[0])
plt.plot(l/1e3,fun2(l)*365.25*Hf,color=colorList[1],linestyle='--',label=nameList[1])
plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]*Hf,[2],alpha=.25,color=colorList[0])
plt.scatter(lengthTime[ind1:ind2]/1e3,UfTime[ind1:ind2]*Hf,[2],alpha=.25,color=colorList[1])
if(ind3>ind2):
	plt.plot(l/1e3,fun3(l)*365.25*Hf,color=colorList[2],linestyle='--',label=nameList[2])
	plt.scatter(lengthTime[ind2:ind3]/1e3,UfTime[ind2:ind3]*Hf,[2],alpha=.25,color=colorList[2])
if(ind4>ind3):
	plt.plot(l/1e3,fun4(l)*365.25*Hf,color=colorList[3],linestyle='--',label=nameList[3])
	plt.scatter(lengthTime[ind3:ind4]/1e3,UfTime[ind3:ind4]*Hf,[2],alpha=.25,color=colorList[3])
if(ind5>ind4):
	plt.plot(l/1e3,fun5(l)*365.25*Hf,color=colorList[4],linestyle='--',label=nameList[4])
	plt.scatter(lengthTime[ind4:ind5]/1e3,UfTime[ind4:ind5]*Hf,[2],alpha=.25,color=colorList[4])
# plt.plot(l/1e3,fluxFast,color='green',linestyle='--')
# plt.plot(l/1e3,fluxSlow,color='green',linestyle='--')

plt.ylim([0,1e7])
plt.title('Flux Diagram')
plt.xlabel('Length [km]')
plt.ylabel('Volume Flux [m^2/year]')

plt.legend()
plt.savefig(f"figs/FluxDiagram_{fileEnding}.png", format='png', dpi=400)


plt.figure(4,figsize=(8, 5))
## See above comment explaining why eq1y is plotted on the x-axis. This is intentional. 



plt.plot(np.insert(loadData1[:,0],0,bTime[ind1-1]),np.insert((loadData1[:,1]+loadData1[:,2])/2.0,0,
	lengthTime[ind1-1]),color=colorList[0],linestyle='--')
plt.plot(bTime[:ind1],lengthTime[:ind1],color=colorList[0],linestyle='-',label=nameList[0])

plt.plot(np.insert(loadData2[:,0],0,bTime[ind2-1]),np.insert((loadData2[:,1]+loadData2[:,2])/2.0,0,
	lengthTime[ind2-1]),color=colorList[1],linestyle='--')
plt.plot(bTime[ind1:ind2],lengthTime[ind1:ind2],color=colorList[1],linestyle='-',label=nameList[1])

if(ind3>ind2):
	plt.plot(np.insert(loadData3[:,0],0,bTime[ind3-1]),np.insert((loadData3[:,1]+loadData2[:,2])/2.0,0,
	lengthTime[ind3-1]),color=colorList[2],linestyle='--')
plt.plot(bTime[ind2:ind3],lengthTime[ind2:ind3],color=colorList[2],linestyle='-',label=nameList[2])
if(ind4>ind3):
	plt.plot(np.insert(loadData4[:,0],0,bTime[ind4-1]),np.insert((loadData4[:,1]+loadData4[:,2])/2.0,0,
	lengthTime[ind4-1]),color=colorList[3],linestyle='--')
plt.plot(bTime[ind3:ind4],lengthTime[ind3:ind4],color=colorList[3],linestyle='-',label=nameList[3])
if(ind5>ind4):
	plt.plot(np.insert(loadData5[:,0],0,bTime[ind5-1]),np.insert((loadData5[:,1]+loadData5[:,2])/2.0,0,
	lengthTime[ind5-1]),color=colorList[4],linestyle='--')
plt.plot(bTime[ind4:ind5],lengthTime[ind4:ind5],color=colorList[4],linestyle='-',label=nameList[4])

plt.legend()
plt.title('Stability Diagram')
plt.ylabel('Length [km]')
plt.xlabel('Melt Rate [m/day]')
plt.savefig(f"figs/Stability_{fileEnding}.png", format='png', dpi=400)
plt.show()
