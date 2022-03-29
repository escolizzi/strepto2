#!/usr/bin/python3

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
mpl.use('QT5Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import seaborn as sns
from matplotlib.lines import Line2D
from collections import OrderedDict
import numpy as np

import seaborn as sns

fig, ax = plt.subplots()
plt.subplots_adjust(wspace=0, hspace=0)

heatmap_colormap = sns.diverging_palette(240, 10, n=9,as_cmap=True) #plt.cm.RdYlBu_r
cmap = heatmap_colormap #plt.cm.bwr

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
# cax = fig.add_axes([0.8, 0.15, 0.05, 0.3])

n_colors = 9
bounds = [x for x in range(n_colors+1)] 

# print("bounds: ", bounds)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

my_colorbar = fig.colorbar(
                mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                cax=cbar_ax, orientation='vertical',
                label="who wins")
             
# tick_locs = (np.arange(len(bounds)) + 0.5)*(len(bounds)-1)/len(bounds)
tick_locs = [ 0.5+x for x in range(n_colors) ]
tick_labels=[ x-4 for x in range(n_colors) ] 
my_colorbar.set_ticks(tick_locs)
my_colorbar.set_ticklabels(tick_labels)

plt.show()