import numpy as np
import os
import sys
from matplotlib import pyplot as plt
# from glaciome1D_dimensional import glaciome, basic_figure, plot_basic_figure, constants
sys.path.append('/Users/psummers8/Documents/glaciome1D')
from glaciome1D import glaciome, basic_figure, plot_basic_figure, constants
from scipy.integrate import trapz
import pickle
import time
import glob
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: 
    Paul Summers
    May 2025
    Example of running glaciome1d

"""

files = sorted(glob.glob('*.pickle'))

plt.figure()

for j in [0,100]:
    file = open(files[j], 'rb')
    data = pickle.load(file)
    file.close()
    # plot_basic_figure(data, axes, color_id, j)

    print(data)


    plt.subplot(2,2,1)
    plt.plot(data.X,data.U/365,linestyle='--',color='xkcd:sunflower')

    plt.subplot(2,2,2)
    plt.plot(data.X_,data.H,linestyle='--',color='xkcd:red')

    plt.subplot(2,2,3)
    plt.plot(data.X_,data.gg,linestyle='--',color='xkcd:apple')

    plt.subplot(2,2,4)
    plt.plot(data.X,data.muW,linestyle='--',color='xkcd:lavender')

    plt.suptitle(f"Mélange Uc: {data.Uc}")
    plt.tight_layout()
plt.savefig("fFTestFigure{data.Uc}.png",format='png',dpi=150)
plt.show()
plt.close()







