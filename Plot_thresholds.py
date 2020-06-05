# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 10:23:45 2020

@author: goodwela
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as col
import pickle


thresh_file='Spatial_thresholds_mm.pickle'                  
f_myfile = open(thresh_file,'rb')
threshvals = pickle.load(f_myfile)
f_myfile.close() 

string = ['DJF Threshold (mm)','MAM Threshold (mm)','JJA Threshold (mm)','SON Threshold (mm)']

fig = plt.figure(figsize=(6.5, 2.5))

for i in range(1,5):
    ax = fig.add_subplot(2, 2, i,frameon=False)
    
    Iplot = threshvals[:,:,i-1]
    Iplot[Iplot == 0] = np.nan
    pic=plt.imshow(np.flipud(Iplot),interpolation = 'none', vmin = 0.3,cmap='jet')

    ax.xaxis.set_visible(False)
    ax.yaxis.set_visible(False)
    plt.colorbar(pic)
    plt.title(string[i-1])


plt.show()
fig.savefig('ThresholdFigures.pdf', bbox_inches='tight')



#for i in range(1,5):
#    plt.subplot(2,2,i,frameon=False)
#    # add a subplot with no frame
#    Iplot = threshvals[:,:,i-1]
#    Iplot[Iplot == 0] = np.nan
#    pic=plt.imshow(np.flipud(Iplot),interpolation = 'none', vmin = 0.3,cmap='jet') 
#    plt.colorbar(pic)
#    pic.axes.get_xaxis().set_visible(False)
#    pic.axes.get_yaxis().set_visible(False)
#    plt.title(string[i-1])