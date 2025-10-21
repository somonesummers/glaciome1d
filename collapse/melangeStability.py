#!/usr/bin/env python
# coding: utf-8
import sys
import numpy as np
from matplotlib import pyplot as plt
sys.path.append('/Users/psummers8/Documents/glaciome1D')
sys.path.append('/storage/home/hcoda1/2/psummers8/glaciome1d')
from glaciome1D import constants, glaciome


import glob
import pickle

files = sorted(glob.glob('sweepLin/*.pickle'))
ind1 = len(files)
files += sorted(glob.glob('sweepU/*.pickle'))
ind2 = len(files)
files += sorted(glob.glob('sweepFlat/*.pickle'))
ind3 = len(files)
files += sorted(glob.glob('sweepLessLinear/*.pickle'))
ind4 = len(files)
files += sorted(glob.glob('sweepBackLinear/*.pickle'))
ind5 = len(files)
# files += sorted(glob.glob('sweepTemp/*.pickle'))
# print(f'index 1,2: {ind1},{ind2}')
toIterate = np.arange(len(files))
H0Time = np.zeros(np.shape(toIterate))
UfTime = np.zeros(np.shape(toIterate))
lengthTime = np.zeros(np.shape(toIterate))
for j in toIterate:
    file = files[j]
    with open(files[j], 'rb') as fileName:
        data = pickle.load(fileName)
        fileName.close()
    H0Time[j]=data.H0
    lengthTime[j]=data.X[-1] - data.X[0]
    UfTime[j]=data.U[-1]

Ut = 6000 #m/yr
Ht = 600 #m
Hf = 25 #m
b  = .75*365.25 #m/yr (quoted in m/day)

uf = np.linspace(50, 350, 100)
l = np.linspace(1, 25000, 100)
Uf, L = np.meshgrid(uf,l)

output = Ht*Ut - Uf*365.25*Hf - b*L
fluxPlug = Ht*Ut/(l/230)*Hf

pwrShift = 5e3

plt.figure(3,figsize=(8, 5))
plt.plot(l/1e3,fluxPlug/Hf/365.25,color='gray')
plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]/365.25,[2],alpha=.25,color='xkcd:red')
g = np.ones([np.shape(lengthTime[:ind1])[0],2])
g[:,1] = np.log(lengthTime[:ind1] + pwrShift)
m1 = np.linalg.pinv(g) @ (np.log(UfTime[:ind1]/365.25))
plt.plot(l/1e3,np.exp(m1[0])*(l + pwrShift)**m1[1],color='red',linestyle='--',label=f'Linear Melt \n  ({np.exp(m1[0]):.3g})(L+{pwrShift:1.1g})^({m1[1]:.3g})')

plt.scatter(lengthTime[ind1:ind2]/1e3,UfTime[ind1:ind2]/365.25,[2],alpha=.25,color='xkcd:blue')
g = np.ones([np.shape(lengthTime[ind1:ind2])[0],2])
g[:,1] = np.log(lengthTime[ind1:ind2] + pwrShift)
m2 = np.linalg.pinv(g) @ (np.log(UfTime[ind1:ind2]/365.25))
plt.plot(l/1e3,np.exp(m2[0])*(l + pwrShift)**m2[1],color='blue',linestyle='--',label=f'U Melt \n  ({np.exp(m2[0]):.3g})(L+{pwrShift:1.1g})^({m2[1]:.3g})')

plt.scatter(lengthTime[ind2:ind3]/1e3,UfTime[ind2:ind3]/365.25,[2],alpha=.25,color='xkcd:cyan')
g = np.ones([np.shape(lengthTime[ind2:ind3])[0],2])
g[:,1] = np.log(lengthTime[ind2:ind3] + pwrShift)
m3 = np.linalg.pinv(g) @ (np.log(UfTime[ind2:ind3]/365.25))
plt.plot(l/1e3,np.exp(m3[0])*(l + pwrShift)**m3[1],color='xkcd:cyan',linestyle='--',label=f'Flat Melt \n  ({np.exp(m3[0]):.3g})(L+{pwrShift:1.1g})^({m3[1]:.3g})')

plt.scatter(lengthTime[ind3:ind4]/1e3,UfTime[ind3:ind4]/365.25,[2],alpha=.25,color='xkcd:indigo')
g = np.ones([np.shape(lengthTime[ind3:ind4])[0],2])
g[:,1] = np.log(lengthTime[ind3:ind4] + pwrShift)
m4 = np.linalg.pinv(g) @ (np.log(UfTime[ind3:ind4]/365.25))
plt.plot(l/1e3,np.exp(m4[0])*(l + pwrShift)**m4[1],color='xkcd:indigo',linestyle='--',label=f'Less Linear Melt \n  ({np.exp(m4[0]):.3g})(L+{pwrShift:1.1g})^({m4[1]:.3g})')

plt.scatter(lengthTime[ind4:ind5]/1e3,UfTime[ind4:ind5]/365.25,[2],alpha=.25,color='xkcd:apple green')
g = np.ones([np.shape(lengthTime[ind4:ind5])[0],2])
g[:,1] = np.log(lengthTime[ind4:ind5] + pwrShift)
m5 = np.linalg.pinv(g) @ (np.log(UfTime[ind4:ind5]/365.25))
plt.plot(l/1e3,np.exp(m5[0])*(l + pwrShift)**m5[1],color='xkcd:apple green',linestyle='--',label=f'Back Linear Melt \n  ({np.exp(m5[0]):.3g})(L+{pwrShift:1.1g})^({m5[1]:.3g})')

# print(f' fit is ax^p a: {m1[0]}, p: {m1[1]}')
# print(f' fit is ax^p a: {m2[0]}, p: {m2[1]}')
# print(f' fit is ax^p a: {m3[0]}, p: {m3[1]}')
# print(f' fit is ax^p a: {m4[0]}, p: {m4[1]}')

plt.xlabel('Length [km]')
plt.ylabel('End Speed [m/day]')
plt.ylim([75,350])
plt.legend()
plt.savefig("Fitting.png", format='png', dpi=400)


plt.figure(1,figsize=(8, 5))
# plt.subplot(121)
cp = plt.contourf(L/1e3,Uf,output,np.linspace(-4e6,4e6,127),cmap='RdBu_r')
cc = plt.contour(L/1e3,Uf,output,[0],linestyles='--',colors='black')
plt.plot(l/1e3,np.exp(m1[0])*(l + pwrShift)**(m1[1]),color='xkcd:red',linestyle='--',label='Linear Melt')
plt.plot(l/1e3,np.exp(m2[0])*(l + pwrShift)**(m2[1]),color='xkcd:blue',linestyle='--',label='U Melt')
plt.plot(l/1e3,np.exp(m3[0])*(l + pwrShift)**(m3[1]),color='xkcd:cyan',linestyle='--',label='Flat Melt')
plt.plot(l/1e3,np.exp(m4[0])*(l + pwrShift)**(m4[1]),color='xkcd:indigo',linestyle='--',label='Less Linear Melt')
plt.plot(l/1e3,np.exp(m5[0])*(l + pwrShift)**(m5[1]),color='xkcd:apple green',linestyle='--',label='Back Linear Melt')
plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]/365.25,[2],alpha=.25,color='xkcd:red')
plt.scatter(lengthTime[ind1:ind2]/1e3,UfTime[ind1:ind2]/365.25,[2],alpha=.25,color='xkcd:blue')
plt.scatter(lengthTime[ind2:ind3]/1e3,UfTime[ind2:ind3]/365.25,[2],alpha=.25,color='xkcd:cyan')
plt.scatter(lengthTime[ind3:ind4]/1e3,UfTime[ind3:ind4]/365.25,[2],alpha=.25,color='xkcd:indigo')
plt.scatter(lengthTime[ind4:ind5]/1e3,UfTime[ind4:ind5]/365.25,[2],alpha=.25,color='xkcd:apple green')
plt.clabel(cc, inline=3, fontsize=8)
cbar = plt.colorbar(cp)
cbar.set_label('$\\partial V / \\partial t$')
plt.xlabel('Length [km]')
plt.ylabel('End Speed [m/day]')
plt.title(f'Change in Volume [m^2/yr] (melt = {b/365.25})')
plt.ylim([np.min(uf),np.max(uf)])
plt.xlim([np.min(l)*1e-3,np.max(l)*1e-3])
plt.legend()
plt.tight_layout()
plt.savefig("dVdt.png", format='png', dpi=400)

# plt.show()


plt.figure(2,figsize=(8, 5))
plt.plot(l/1e3,fluxPlug,color='gray',label='Previous estimate of discharge')
for b in np.linspace(.2,2,10)*365.25:
	plt.plot(l/1e3,Ut*Ht - b*l,color='maroon',alpha=.1,linewidth=.5)#label=f'Flux in - Melt (b = {b/365.25:.2f})')

plt.plot(l/1e3,np.exp(m1[0])*(l + pwrShift)**(m1[1])*365.25*Hf,color='xkcd:red',linestyle='--',label='Linear Melt')
plt.plot(l/1e3,np.exp(m2[0])*(l + pwrShift)**(m2[1])*365.25*Hf,color='xkcd:blue',linestyle='--',label='U Melt')
plt.plot(l/1e3,np.exp(m3[0])*(l + pwrShift)**(m3[1])*365.25*Hf,color='xkcd:cyan',linestyle='--',label='Flat Melt')
plt.plot(l/1e3,np.exp(m4[0])*(l + pwrShift)**(m4[1])*365.25*Hf,color='xkcd:indigo',linestyle='--',label='Less Linear Melt')
plt.plot(l/1e3,np.exp(m5[0])*(l + pwrShift)**(m5[1])*365.25*Hf,color='xkcd:apple green',linestyle='--',label='Back Linear Melt')
# plt.plot(l/1e3,fluxFast,color='green',linestyle='--')
# plt.plot(l/1e3,fluxSlow,color='green',linestyle='--')
plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]*Hf,[2],alpha=.25,color='xkcd:red')
plt.scatter(lengthTime[ind1:ind2]/1e3,UfTime[ind1:ind2]*Hf,[2],alpha=.25,color='xkcd:blue')
plt.scatter(lengthTime[ind2:ind3]/1e3,UfTime[ind2:ind3]*Hf,[2],alpha=.25,color='xkcd:cyan')
plt.scatter(lengthTime[ind3:ind4]/1e3,UfTime[ind3:ind4]*Hf,[2],alpha=.25,color='xkcd:indigo')
plt.scatter(lengthTime[ind4:ind5]/1e3,UfTime[ind4:ind5]*Hf,[2],alpha=.25,color='xkcd:apple green')
plt.ylim([0,1e7])
plt.xlabel('Length [km]')
plt.ylabel('Volume Flux [m^2/year]')

plt.legend()
plt.savefig("FluxDiagram.png", format='png', dpi=400)


l = np.linspace(1, 25000, 250)

def eq(b,length,m):
	return np.exp(m[0])*(length  + pwrShift )**(m[1])*365.25*(25) + b*365.25*length - 6000*600

eq1x = []
eq1y = []
eq2x = []
eq2y = []
eq3x = []
eq3y = []
eq4x = []
eq4y = []
eq5x = []
eq5y = []
# We take the product to find a zero crossing, then we a linear root finder
for b in np.linspace(.2,1.6,250):
	for li in range(1,np.shape(l)[0]):
		if(eq(b,l[li],m3) * eq(b,l[li-1],m3) < 0):
			rt = l[li-1] + eq(b,l[li],m3)/(eq(b,l[li-1],m3) - eq(b,l[li],m3))* (l[li]-l[li-1])
			eq3x.append(rt)
			eq3y.append(b)
		if(eq(b,l[li],m2) * eq(b,l[li-1],m2) < 0):
			rt = l[li-1] + eq(b,l[li],m2)/(eq(b,l[li-1],m2) - eq(b,l[li],m2))* (l[li]-l[li-1])
			eq2x.append(rt)
			eq2y.append(b)
		if(eq(b,l[li],m1) * eq(b,l[li-1],m1) < 0):
			rt = l[li-1] + eq(b,l[li],m1)/(eq(b,l[li-1],m1) - eq(b,l[li],m1))* (l[li]-l[li-1])
			eq1x.append(rt)
			eq1y.append(b)
		if(eq(b,l[li],m4) * eq(b,l[li-1],m4) < 0):
			rt = l[li-1] + eq(b,l[li],m4)/(eq(b,l[li-1],m4) - eq(b,l[li],m4))* (l[li]-l[li-1])
			eq4x.append(rt)
			eq4y.append(b)
		if(eq(b,l[li],m5) * eq(b,l[li-1],m5) < 0):
			rt = l[li-1] + eq(b,l[li],m5)/(eq(b,l[li-1],m5) - eq(b,l[li],m5))* (l[li]-l[li-1])
			eq5x.append(rt)
			eq5y.append(b)
eq1x = np.asarray(eq1x)
eq1y = np.asarray(eq1y)
eq2x = np.asarray(eq2x) 
eq2y = np.asarray(eq2y) 
eq3x = np.asarray(eq3x) 
eq3y = np.asarray(eq3y) 
eq4x = np.asarray(eq4x) 
eq4y = np.asarray(eq4y) 
eq5x = np.asarray(eq5x) 
eq5y = np.asarray(eq5y) 

sorted_indices = np.argsort(eq1x)
eq1x = eq1x[sorted_indices]
eq1y = eq1y[sorted_indices]
maxInd1 = np.argmax(eq1y)

sorted_indices = np.argsort(eq2x)
eq2x = eq2x[sorted_indices]
eq2y = eq2y[sorted_indices]
maxInd2 = np.argmax(eq2y)

sorted_indices = np.argsort(eq3x)
eq3x = eq3x[sorted_indices]
eq3y = eq3y[sorted_indices]
maxInd3 = np.argmax(eq3y)

sorted_indices = np.argsort(eq4x)
eq4x = eq4x[sorted_indices]
eq4y = eq4y[sorted_indices]
maxInd4 = np.argmax(eq4y)

sorted_indices = np.argsort(eq5x)
eq5x = eq5x[sorted_indices]
eq5y = eq5y[sorted_indices]
maxInd5 = np.argmax(eq5y)

plt.figure(4,figsize=(8, 5))

# plt.plot(l,eq1(.3,l))

plt.plot([np.min(l),np.max(l)],[.422, .422],color='xkcd:red',linestyle=':')
plt.plot([np.min(l),np.max(l)],[.425, .425],color='xkcd:blue',linestyle=':')
plt.plot([np.min(l),np.max(l)],[.2, .2],color='xkcd:cyan',linestyle=':')
plt.plot([np.min(l),np.max(l)],[.480, .480],color='xkcd:indigo',linestyle=':')
plt.plot([np.min(l),np.max(l)],[.2, .2],color='xkcd:apple green',linestyle=':')

plt.plot(eq1x[:maxInd1+1],eq1y[:maxInd1+1],color='xkcd:red',linestyle='--')
plt.plot(eq2x[:maxInd2+1],eq2y[:maxInd2+1],color='xkcd:blue',linestyle='--')
plt.plot(eq3x[:maxInd3+1],eq3y[:maxInd3+1],color='xkcd:cyan',linestyle='--')
plt.plot(eq4x[:maxInd4+1],eq4y[:maxInd4+1],color='xkcd:indigo',linestyle='--')
plt.plot(eq5x[:maxInd5+1],eq5y[:maxInd5+1],color='xkcd:apple green',linestyle='--')
plt.plot(eq1x[maxInd1:],eq1y[maxInd1:],color='xkcd:red',label='Linear Melt')
plt.plot(eq2x[maxInd2:],eq2y[maxInd2:],color='xkcd:blue',label='U Melt')
plt.plot(eq3x[maxInd3:],eq3y[maxInd3:],color='xkcd:cyan',label='Flat Melt')
plt.plot(eq4x[maxInd4:],eq4y[maxInd4:],color='xkcd:indigo',label='Less Linear Melt')
plt.plot(eq5x[maxInd5:],eq5y[maxInd5:],color='xkcd:apple green',label='Back Linear Melt')
plt.legend()
plt.xlabel('Length [km]')
plt.ylabel('Melt Rate [m/day]')
plt.savefig("Stability.png", format='png', dpi=400)
plt.show()
