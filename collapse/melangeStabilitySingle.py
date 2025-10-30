#!/usr/bin/env python
# coding: utf-8
import sys
import numpy as np
from matplotlib import pyplot as plt
sys.path.append('/Users/psummers8/Documents/glaciome1D')
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
from glaciome1D import constants, glaciome
import scipy
import matplotlib as mpl
import glob
import pickle
import cmocean

folders = ['LL6000_slow']
fileEnding = 'LL6000'
filePrefix = 'figs/Single'
UcList = np.array([6000])


files = sorted(glob.glob(f'PACE/{folders[0]}/swee*.pickle'))
ind1 = len(files)
# files += sorted(glob.glob(f'PACE/{folders[1]}/swee*.pickle'))
ind2 = len(files)
# files += sorted(glob.glob(f'PACE/{folders[2]}/swee*.pickle'))
ind3 = len(files)
# files += sorted(glob.glob(f'PACE/{folders[3]}/sw*.pickle'))
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

#cleaning
index = len(lengthTime)-1
overGrab = 1
if(np.min(lengthTime) < 1000 or np.min(H0Time) < 24):
	print([np.argmin(abs(lengthTime-1000)), np.argmin(abs(H0Time-24))])
	index = np.min([np.argmin(abs(lengthTime-1000)), np.argmin(abs(H0Time-24))])
	toIterateTmp = toIterate[:index + overGrab]
	H0TimeTmp = H0Time[:index + overGrab]
	UfTimeTmp = UfTime[:index + overGrab]
	bTimeTmp = bTime[:index + overGrab]
	lengthTimeTmp = lengthTime[:index + overGrab]
elif(False):
	ind1 += 1
	toIterateTmp = np.append(toIterate, toIterate[-1]+1)
	H0TimeTmp = np.append(H0Time, 2*H0Time[-1] -   H0Time[-2]) #linear
	UfTimeTmp = np.append(UfTime, UcList[0]*600/25 + 5*365.25)
	bTimeTmp  = np.append(bTime,  2*bTime[-1]  -    bTime[-2]) #linear
	lengthTimeTmp = np.append(lengthTime, 0)
	print(f'Appending record with false post collapse at (Uf,Length): ({UfTimeTmp[-1]},{lengthTimeTmp[-1]})')
else:
	toIterateTmp = toIterate
	H0TimeTmp = H0Time
	UfTimeTmp = UfTime
	bTimeTmp  = bTime
	lengthTimeTmp = np.append(lengthTime, 0)
print(f'Collapse at index {index} of {len(files)-1}, melt {bTime[index]:.04f}, overplot {overGrab}')
toIterate = toIterateTmp
H0Time = H0TimeTmp
UfTime = UfTimeTmp
lengthTime = lengthTimeTmp
bTime =  bTimeTmp


Ut = UcList[0] #m/yr
Ht = 600 #m
Hf = 25 #m
bMelt  = bTime[0]*365.25 #m/yr 

uf = np.linspace(50, 400, 100)
l = np.linspace(1, 25000, 100)
Uf, L = np.meshgrid(uf,l)

output = Ht*Ut - Uf*365.25*Hf - bMelt*L
fluxPlug = Ht*Ut/(l/230)*Hf

plt.figure(3,figsize=(8, 5))
# plt.plot(l/1e3,fluxPlug/Hf/365.25,color='gray')
plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]/365.25,[2],alpha=1,color=colorList[0])
fun1 = scipy.interpolate.interp1d(lengthTime[:ind1],UfTime[:ind1]/365.25,fill_value='extrapolate')
plt.plot(l/1e3,fun1(l),color=colorList[0],linestyle='-',alpha=.2,label=nameList[0])

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
plt.savefig(f"{filePrefix}Fitting_{fileEnding}.png", format='png', dpi=400)
# plt.show()

plt.figure(1,figsize=(8, 5))
# plt.subplot(121)
cp = plt.contourf(L/1e3,Uf,output,np.linspace(-4e6,4e6,127),cmap='RdBu')
cc = plt.contour(L/1e3,Uf,output,[0],linestyles='--',colors='black',alpha=.5)
cc2 = plt.contour(L/1e3,Uf,np.abs(L*L/500/(output+1e3))*365.25,[10,100,1000],colors='purple',linestyles='--',alpha=.25)
plt.plot(l/1e3,fun1(l),color=colorList[0],linestyle='-',label=nameList[0])
# plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]/365.25,[2],alpha=.25,color=colorList[0])

plt.clabel(cc, inline=3, fontsize=8)
plt.clabel(cc2, inline=3, fontsize=8)
cbar = plt.colorbar(cp)
cbar.set_label('$\\partial V / \\partial t$')
plt.xlabel('Length [km]')
plt.ylabel('End Speed [m/day]')
plt.title(f'Change in Volume [m^2/yr] (melt = {bMelt/365.25:.3f})')
plt.ylim([np.min(uf),np.max(uf)])
plt.xlim([np.min(l)*1e-3,np.max(l)*1e-3])
plt.legend()
plt.tight_layout()
plt.savefig(f"{filePrefix}dVdt_{fileEnding}.png", format='png', dpi=400)

# plt.show()


plt.figure(2,figsize=(8, 5))
# plt.plot(l/1e3,fluxPlug,color='gray',label='Previous estimate of discharge')
mapName = 'RdYlBu_r'
n_lines = 11
colorMapForLine = mpl.colormaps[mapName]
meltRates = np.linspace(bTime[0]*.9,bTime[-1]*1.1,n_lines,endpoint=True)
fluxColors = colorMapForLine(np.linspace(0,1,n_lines))
for j in range(n_lines):
	b = meltRates[j]
	plt.plot(l/1e3,UcList[0]*Ht - b*l*365.25,color=fluxColors[j],alpha=1,linewidth=1)#label=f'Flux in - Melt (b = {b/365.25:.2f})')
# plt.scatter(np.zeros(n_lines), meltRates, s=30, c=meltRates, cmap='Reds')



plt.plot(l/1e3,fun1(l)*365.25*Hf,color='black',alpha=.25,linestyle='-',label=nameList[0])
plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]*Hf,[10],c=bTime[:ind1],cmap=mapName)
plt.clim([bTime[0]*.9,bTime[-1]*1.1])
cbar = plt.colorbar()
cbar.set_label('Melt Rate [m/day]')

plt.ylim([np.min(fun1(l))*.9*365.25*Hf,np.max(fun1(l))*1.1*365.25*Hf])
plt.title('Flux Diagram')
plt.xlabel('Length [km]')
plt.ylabel('Volume Flux [m^2/year]')

plt.legend()
plt.savefig(f"{filePrefix}FluxDiagram_{fileEnding}.png", format='png', dpi=400)


l = np.linspace(1, 25000, 300)

# def eq(b,length,m,calvSpd):
# 	return np.exp(m[0])*(length  + pwrShift )**(m[1])*365.25*(25) + b*365.25*length - calvSpd*600

## So I originally had the wrong thing on the x/y axis, so we plot 'eqx' on the y axis below to ensure the variable we change (melt) is on the x-axis. 
# Please forgive me, for I know I have erred, but to leave a code comment explaining this almost as good as getting it right...right?

def eq(b,length,fun,calvSpd):
	return fun(length)*365.25*(25) + b*365.25*length - calvSpd*600
eq1x = []
eq1y = []

print(UcList)
# We take the product to find a zero crossing, then we a linear root finder

for b in np.linspace(.2,1,300):
	for li in range(1,np.shape(l)[0]):
		if(eq(b,l[li],fun1,UcList[0]) * eq(b,l[li-1],fun1,UcList[0]) < 0):
			rt = l[li-1] + eq(b,l[li],fun1,UcList[0])/(eq(b,l[li-1],fun1,UcList[0]) - eq(b,l[li],fun1,UcList[0]))* (l[li]-l[li-1])
			eq1x.append(rt)
			eq1y.append(b)
		# if(eqFun(b,l[li],fun1,UcList[0]) * eqFun(b,l[li-1],fun1,UcList[0]) < 0):
		# 	rt = l[li-1] + eqFun(b,l[li],fun1,UcList[0])/(eqFun(b,l[li-1],fun1,UcList[0]) - eqFun(b,l[li],fun1,UcList[0]))* (l[li]-l[li-1])
		# 	eq1x.append(rt)
		# 	eq1y.append(b)

eq1x = np.asarray(eq1x)
eq1y = np.asarray(eq1y)

sorted_indices = np.argsort(eq1x)
eq1x = eq1x[sorted_indices]
eq1y = eq1y[sorted_indices]
maxInd1 = np.argmax(eq1y)

plt.figure(4,figsize=(8, 5))
## See above comment explaining why eq1y is plotted on the x-axis. This is intentional. 

plt.plot(eq1y[:maxInd1+1],eq1x[:maxInd1+1],color=colorList[0],linestyle='--')

plt.plot(eq1y[maxInd1:],eq1x[maxInd1:],color=colorList[0],linestyle='-',label=nameList[0])


plt.legend()
plt.title('Stability Diagram')
plt.ylabel('Length [km]')
plt.xlabel('Melt Rate [m/day]')
plt.savefig(f"{filePrefix}Stability_{fileEnding}.png", format='png', dpi=400)
plt.show()
