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

folders = []

files = sorted(glob.glob('PACE/L6000/swee*.pickle'))
ind1 = len(files)
files += sorted(glob.glob('PACE/L6000a05/swee*.pickle'))
ind2 = len(files)
files += sorted(glob.glob('PACE/L6000a25/swee*.pickle'))
ind3 = len(files)
files += sorted(glob.glob('PACE/F6000/sw*.pickle'))
ind4 = len(files)
# files += sorted(glob.glob('PACE/L6000a25/sw*.pickle'))
ind5 = len(files)
print(ind1,ind2,ind3,ind4,ind5)
colorList = ['xkcd:red','xkcd:green','xkcd:blue','xkcd:indigo','xkcd:lavendar']
nameList = ['F6000','L6000','L6000a25','F6000','']
UcList = np.array([6000,6330,7800,6000])
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

Ut = 6000 #m/yr
Ht = 600 #m
Hf = 25 #m
b  = 0.35*365.25 #m/yr (quoted in m/day)

uf = np.linspace(50, 400, 100)
l = np.linspace(1, 25000, 100)
Uf, L = np.meshgrid(uf,l)

output = Ht*Ut - Uf*365.25*Hf - b*L
fluxPlug = Ht*Ut/(l/230)*Hf

pwrShift = 50e3
fitLineL = 25

plt.figure(3,figsize=(8, 5))
# plt.plot(l/1e3,fluxPlug/Hf/365.25,color='gray')
plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]/365.25,[2],alpha=1,color=colorList[0])
g = np.ones([np.shape(lengthTime[:ind1])[0],2])
g[:,1] = np.log(lengthTime[:ind1] + pwrShift)
m1 = np.linalg.pinv(g) @ (np.log(UfTime[:ind1]/365.25))
plt.plot(l/1e3,np.exp(m1[0])*(l + pwrShift)**m1[1],color=colorList[0],linestyle='--',alpha=.5,
	label=f'{nameList[0]} \n  ({np.exp(m1[0]):.3g})(L+{pwrShift:1.1g})^({m1[1]:.3g})')
g = np.ones([np.shape(lengthTime[ind1-fitLineL:ind1])[0],2])
g[:,1] = lengthTime[ind1-fitLineL:ind1]
l_fit = np.linalg.pinv(g) @ (UfTime[ind1-fitLineL:ind1]/365.25)
plt.plot(l/1e3,l_fit[0] + l*l_fit[1],color=colorList[0],linestyle=':',alpha=.2)
fun1 = scipy.interpolate.interp1d(lengthTime[:ind1],UfTime[:ind1]/365.25,fill_value='extrapolate')
plt.plot(l/1e3,fun1(l),color='gray',linestyle='-',alpha=.2)

plt.scatter(lengthTime[ind1:ind2]/1e3,UfTime[ind1:ind2]/365.25,[2],alpha=1,color=colorList[1])
g = np.ones([np.shape(lengthTime[ind1:ind2])[0],2])
g[:,1] = np.log(lengthTime[ind1:ind2] + pwrShift)
m2 = np.linalg.pinv(g) @ (np.log(UfTime[ind1:ind2]/365.25))
plt.plot(l/1e3,np.exp(m2[0])*(l + pwrShift)**m2[1],color=colorList[1],linestyle='--',alpha=.5,
	label=f'{nameList[1]} \n  ({np.exp(m2[0]):.3g})(L+{pwrShift:1.1g})^({m2[1]:.3g})')
g = np.ones([np.shape(lengthTime[ind2-fitLineL:ind2])[0],2])
g[:,1] = lengthTime[ind2-fitLineL:ind2]
l_fit = np.linalg.pinv(g) @ (UfTime[ind2-fitLineL:ind2]/365.25)
plt.plot(l/1e3,l_fit[0] + l*l_fit[1],color=colorList[1],linestyle=':',alpha=.2)

if(ind3>ind2):
	plt.scatter(lengthTime[ind2:ind3]/1e3,UfTime[ind2:ind3]/365.25,[2],alpha=1,color=colorList[2])
	g = np.ones([np.shape(lengthTime[ind2:ind3])[0],2])
	g[:,1] = np.log(lengthTime[ind2:ind3] + pwrShift)
	m3 = np.linalg.pinv(g) @ (np.log(UfTime[ind2:ind3]/365.25))
	plt.plot(l/1e3,np.exp(m3[0])*(l + pwrShift)**m3[1],color=colorList[2],linestyle='--',alpha=.5,
		label=f'{nameList[2]} \n  ({np.exp(m3[0]):.3g})(L+{pwrShift:1.1g})^({m3[1]:.3g})')
	g = np.ones([np.shape(lengthTime[ind3-fitLineL:ind3])[0],2])
	g[:,1] = lengthTime[ind3-fitLineL:ind3]
	l_fit = np.linalg.pinv(g) @ (UfTime[ind3-fitLineL:ind3]/365.25)
	plt.plot(l/1e3,l_fit[0] + l*l_fit[1],color=colorList[2],linestyle=':',alpha=.2)

if(ind4>ind3):
	plt.scatter(lengthTime[ind3:ind4]/1e3,UfTime[ind3:ind4]/365.25,[2],alpha=1,color=colorList[3])
	g = np.ones([np.shape(lengthTime[ind3:ind4])[0],2])
	g[:,1] = np.log(lengthTime[ind3:ind4] + pwrShift)
	m4 = np.linalg.pinv(g) @ (np.log(UfTime[ind3:ind4]/365.25))
	plt.plot(l/1e3,np.exp(m4[0])*(l + pwrShift)**m4[1],color=colorList[3],linestyle='--',alpha=.5,
		label=f'{nameList[3]} \n  ({np.exp(m4[0]):.3g})(L+{pwrShift:1.1g})^({m4[1]:.3g})')
	g = np.ones([np.shape(lengthTime[ind4-fitLineL:ind4])[0],2])
	g[:,1] = lengthTime[ind4-fitLineL:ind4]
	l_fit = np.linalg.pinv(g) @ (UfTime[ind4-fitLineL:ind4]/365.25)
	plt.plot(l/1e3,l_fit[0] + l*l_fit[1],color=colorList[3],linestyle=':',alpha=.2)

if(ind5>ind4):
	plt.scatter(lengthTime[ind4:ind5]/1e3,UfTime[ind4:ind5]/365.25,[2],alpha=1,color=colorList[4])
	g = np.ones([np.shape(lengthTime[ind4:ind5])[0],2])
	g[:,1] = np.log(lengthTime[ind4:ind5] + pwrShift)
	m5 = np.linalg.pinv(g) @ (np.log(UfTime[ind4:ind5]/365.25))
	plt.plot(l/1e3,np.exp(m5[0])*(l + pwrShift)**m5[1],color=colorList[4],linestyle='--',alpha=.5,
		label=f'{nameList[4]} \n  ({np.exp(m5[0]):.3g})(L+{pwrShift:1.1g})^({m5[1]:.3g})')
	g = np.ones([np.shape(lengthTime[ind5-fitLineL:ind5])[0],2])
	g[:,1] = lengthTime[ind5-fitLineL:ind5]
	l_fit = np.linalg.pinv(g) @ (UfTime[ind5-fitLineL:ind5]/365.25)
	plt.plot(l/1e3,l_fit[0] + l*l_fit[1],color=colorList[4],linestyle=':',alpha=.2)

for i in range(len(UcList)):
	plt.plot(0*UcList[i],Ht*UcList[i]/Hf/365.25,marker='*',color=colorList[i],linewidth=0)
# print(f' fit is ax^p a: {m1[0]}, p: {m1[1]}')
# print(f' fit is ax^p a: {m2[0]}, p: {m2[1]}')
# print(f' fit is ax^p a: {m3[0]}, p: {m3[1]}')
# print(f' fit is ax^p a: {m4[0]}, p: {m4[1]}')
plt.title('Fitting')
plt.xlabel('Length [km]')
plt.ylabel('End Speed [m/day]')
plt.ylim([75,np.max([300,np.max(UcList)*Ht/Hf/365.25+20])])
plt.legend()
plt.savefig("Fitting.png", format='png', dpi=400)
# plt.show()

plt.figure(1,figsize=(8, 5))
# plt.subplot(121)
cp = plt.contourf(L/1e3,Uf,output,np.linspace(-4e6,4e6,127),cmap='RdBu_r')
cc = plt.contour(L/1e3,Uf,output,[0],linestyles='--',colors='black')
plt.plot(l/1e3,np.exp(m1[0])*(l + pwrShift)**(m1[1]),color=colorList[0],linestyle='--',label=nameList[0])
plt.plot(l/1e3,np.exp(m2[0])*(l + pwrShift)**(m2[1]),color=colorList[1],linestyle='--',label=nameList[1])
plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]/365.25,[2],alpha=.25,color=colorList[0])
plt.scatter(lengthTime[ind1:ind2]/1e3,UfTime[ind1:ind2]/365.25,[2],alpha=.25,color=colorList[1])
if(ind3>ind2):
	plt.plot(l/1e3,np.exp(m3[0])*(l + pwrShift)**(m3[1]),color=colorList[2],linestyle='--',label=nameList[2])
	plt.scatter(lengthTime[ind2:ind3]/1e3,UfTime[ind2:ind3]/365.25,[2],alpha=.25,color=colorList[2])
if(ind4>ind3):
	plt.plot(l/1e3,np.exp(m4[0])*(l + pwrShift)**(m4[1]),color=colorList[3],linestyle='--',label=nameList[3])
	plt.scatter(lengthTime[ind3:ind4]/1e3,UfTime[ind3:ind4]/365.25,[2],alpha=.25,color=colorList[3])
if(ind5>ind4):
	plt.plot(l/1e3,np.exp(m5[0])*(l + pwrShift)**(m5[1]),color=colorList[4],linestyle='--',label=nameList[4])
	plt.scatter(lengthTime[ind4:ind5]/1e3,UfTime[ind4:ind5]/365.25,[2],alpha=.25,color=colorList[4])
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
for i in range(len(UcList)):
	if(i == 0 or UcList[i] != UcList[0]):
		for b in np.linspace(.35,.6,15)*365.25:
			plt.plot(l/1e3,UcList[i]*Ht - b*l,color=colorList[i],alpha=.2,linewidth=.5)#label=f'Flux in - Melt (b = {b/365.25:.2f})')


plt.plot(l/1e3,np.exp(m1[0])*(l + pwrShift)**(m1[1])*365.25*Hf,color=colorList[0],linestyle='--',label=nameList[0])
plt.plot(l/1e3,np.exp(m2[0])*(l + pwrShift)**(m2[1])*365.25*Hf,color=colorList[1],linestyle='--',label=nameList[1])
plt.scatter(lengthTime[:ind1]/1e3,UfTime[:ind1]*Hf,[2],alpha=.25,color=colorList[0])
plt.scatter(lengthTime[ind1:ind2]/1e3,UfTime[ind1:ind2]*Hf,[2],alpha=.25,color=colorList[1])
if(ind3>ind2):
	plt.plot(l/1e3,np.exp(m3[0])*(l + pwrShift)**(m3[1])*365.25*Hf,color=colorList[2],linestyle='--',label=nameList[2])
	plt.scatter(lengthTime[ind2:ind3]/1e3,UfTime[ind2:ind3]*Hf,[2],alpha=.25,color=colorList[2])
if(ind4>ind3):
	plt.plot(l/1e3,np.exp(m4[0])*(l + pwrShift)**(m4[1])*365.25*Hf,color=colorList[3],linestyle='--',label=nameList[3])
	plt.scatter(lengthTime[ind3:ind4]/1e3,UfTime[ind3:ind4]*Hf,[2],alpha=.25,color=colorList[3])
if(ind5>ind4):
	plt.plot(l/1e3,np.exp(m5[0])*(l + pwrShift)**(m5[1])*365.25*Hf,color=colorList[4],linestyle='--',label=nameList[4])
	plt.scatter(lengthTime[ind4:ind5]/1e3,UfTime[ind4:ind5]*Hf,[2],alpha=.25,color=colorList[4])
# plt.plot(l/1e3,fluxFast,color='green',linestyle='--')
# plt.plot(l/1e3,fluxSlow,color='green',linestyle='--')

plt.ylim([0,1e7])
plt.title('Flux Diagram')
plt.xlabel('Length [km]')
plt.ylabel('Volume Flux [m^2/year]')

plt.legend()
plt.savefig("FluxDiagram.png", format='png', dpi=400)


l = np.linspace(1, 25000, 250)

def eq(b,length,m,calvSpd):
	return np.exp(m[0])*(length  + pwrShift )**(m[1])*365.25*(25) + b*365.25*length - calvSpd*600

def eqFun(b,length,fun,calvSpd):
	return fun(length)*365.25*(25) + b*365.25*length - calvSpd*600
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
		if(eq(b,l[li],m1,UcList[0]) * eq(b,l[li-1],m1,UcList[0]) < 0):
			rt = l[li-1] + eq(b,l[li],m1,UcList[0])/(eq(b,l[li-1],m1,UcList[0]) - eq(b,l[li],m1,UcList[0]))* (l[li]-l[li-1])
			eq1x.append(rt)
			eq1y.append(b)
		# if(eqFun(b,l[li],fun1,UcList[0]) * eqFun(b,l[li-1],fun1,UcList[0]) < 0):
		# 	rt = l[li-1] + eqFun(b,l[li],fun1,UcList[0])/(eqFun(b,l[li-1],fun1,UcList[0]) - eqFun(b,l[li],fun1,UcList[0]))* (l[li]-l[li-1])
		# 	eq1x.append(rt)
		# 	eq1y.append(b)
		if(eq(b,l[li],m2,UcList[1]) * eq(b,l[li-1],m2,UcList[1]) < 0):
			rt = l[li-1] + eq(b,l[li],m2,UcList[1])/(eq(b,l[li-1],m2,UcList[1]) - eq(b,l[li],m2,UcList[1]))* (l[li]-l[li-1])
			eq2x.append(rt)
			eq2y.append(b)
		if(ind3 > ind2 and eq(b,l[li],m3,UcList[2]) * eq(b,l[li-1],m3,UcList[2]) < 0):
			rt = l[li-1] + eq(b,l[li],m3,UcList[2])/(eq(b,l[li-1],m3,UcList[2]) - eq(b,l[li],m3,UcList[2]))* (l[li]-l[li-1])
			eq3x.append(rt)
			eq3y.append(b)
		if(ind4>ind3 and eq(b,l[li],m4,UcList[3]) * eq(b,l[li-1],m4,UcList[3]) < 0):
			rt = l[li-1] + eq(b,l[li],m4,UcList[3])/(eq(b,l[li-1],m4,UcList[3]) - eq(b,l[li],m4,UcList[3]))* (l[li]-l[li-1])
			eq4x.append(rt)
			eq4y.append(b)
		if(ind5 > ind4 and eq(b,l[li],m5,UcList[4]) * eq(b,l[li-1],m5,UcList[4]) < 0):
			rt = l[li-1] + eq(b,l[li],m5,UcList[4])/(eq(b,l[li-1],m5,UcList[4]) - eq(b,l[li],m5,UcList[4]))* (l[li]-l[li-1])
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

if(ind3 > ind2):
	sorted_indices = np.argsort(eq3x)
	eq3x = eq3x[sorted_indices]
	eq3y = eq3y[sorted_indices]
	maxInd3 = np.argmax(eq3y)

if(ind4 > ind3):
	sorted_indices = np.argsort(eq4x)
	eq4x = eq4x[sorted_indices]
	eq4y = eq4y[sorted_indices]
	maxInd4 = np.argmax(eq4y)

if(ind5 > ind4):
	sorted_indices = np.argsort(eq5x)
	eq5x = eq5x[sorted_indices]
	eq5y = eq5y[sorted_indices]
	maxInd5 = np.argmax(eq5y)

plt.figure(4,figsize=(8, 5))

# plt.plot(l,eq1(.3,l))

plt.plot([np.min(l),np.max(l)],[bTime[ind1-5], bTime[ind1-5]],color=colorList[0],linestyle=':')
plt.plot([np.min(l),np.max(l)],[bTime[ind2-5], bTime[ind2-5]],color=colorList[1],linestyle=':')
if(ind3>ind2):
	plt.plot([np.min(l),np.max(l)],[bTime[ind3-5], bTime[ind3-5]],color=colorList[2],linestyle=':')
if(ind4>ind3):
	plt.plot([np.min(l),np.max(l)],[bTime[ind4-5], bTime[ind4-5]],color=colorList[3],linestyle=':')
if(ind5>ind4):
	plt.plot([np.min(l),np.max(l)],[bTime[ind5-5], bTime[ind5-5]],color=colorList[4],linestyle=':')

plt.plot(eq1x[:maxInd1+1],eq1y[:maxInd1+1],color=colorList[0],linestyle='--')
plt.plot(eq2x[:maxInd2+1],eq2y[:maxInd2+1],color=colorList[1],linestyle='--')

plt.plot(eq1x[maxInd1:],eq1y[maxInd1:],color=colorList[0],linestyle='-',label=nameList[0])
plt.plot(eq2x[maxInd2:],eq2y[maxInd2:],color=colorList[1],linestyle='-',label=nameList[1])
if(ind3>ind2):
	plt.plot(eq3x[:maxInd3+1],eq3y[:maxInd3+1],color=colorList[2],linestyle='--')
	plt.plot(eq3x[maxInd3:],eq3y[maxInd3:],color=colorList[2],linestyle='-',label=nameList[2])
if(ind4>ind3):
	plt.plot(eq4x[:maxInd4+1],eq4y[:maxInd4+1],color=colorList[3],linestyle='--')
	plt.plot(eq4x[maxInd4:],eq4y[maxInd4:],color=colorList[3],linestyle='-',label=nameList[3])
if(ind5>ind4):
	plt.plot(eq5x[:maxInd5+1],eq5y[:maxInd5+1],color=colorList[4],linestyle='--')
	plt.plot(eq5x[maxInd5:],eq5y[maxInd5:],color=colorList[4],linestyle='-',label=nameList[4])

plt.legend()
plt.ylim([0,.95])
plt.title('Stability Diagram')
plt.xlabel('Length [km]')
plt.ylabel('Melt Rate [m/day]')
plt.savefig("Stability.png", format='png', dpi=400)
plt.show()
