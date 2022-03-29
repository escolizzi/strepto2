#!/usr/bin/python3

'''
Gets division of labor data from break prob one season, plots it!

'''

import matplotlib as mpl
mpl.use('Qt5Agg')

import sys,math,os,string
from subprocess import Popen, PIPE
from collections import Counter

#from PIL import Image
# import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import seaborn as sns
import numpy as np

ddata={}

#data is written as
# filename (from which you get break prob and brinfl)
# X axis
# ab
# repl

with open(sys.argv[1]) as fin:
    lines=fin.readlines()
    for i in range(0, len(lines), 4):
        filename, lX,lab,lrepl = lines[i:i + 4]
        breakprob_pos = filename.find('breakprob')
        breakprob= filename[breakprob_pos+len('breakprob'):].split('_')[0]
        
        breakinfl_pos = filename.find('breakinfl')
        breakinfl= filename[breakinfl_pos+len('breakinfl'):].split('_')[0][:-4] # to remove .txt
        
        # if breakinfl == '0.0001':
        #     print("Got it")
        
        lX = [int(x) for x in lX[6:-2].split(', ')]
        lab = [float(x) for x in lab[16:-2].split(', ')]
        lrepl = [float(x) for x in lrepl[16:-2].split(', ')]
        if breakprob not in ddata:
            ddata[breakprob]={}
        if breakinfl not in ddata[breakprob]:
            ddata[breakprob][breakinfl]=[lX,lab,lrepl]
        else:
            print("Error, already got this")

# for i,breakprob in enumerate(sorted(ddata)):
#     for j,breakinfl in enumerate(sorted(ddata[breakprob])):
#         print(breakprob,breakinfl)

def find_median(ll):
        sum_lab = sum(ll)
        run_sum = 0.
        for i,ab in enumerate(ll):
            run_sum += ab/sum_lab
            if run_sum>0.5:
                return i

all_data = True
if all_data == True:
    fig = plt.figure()
    spec = gridspec.GridSpec(nrows=5, ncols=4, figure=fig)
    for i,breakprob in enumerate(sorted(ddata,reverse=True)):
        for j,breakinfl in enumerate(sorted(ddata[breakprob])):
            ax = fig.add_subplot(spec[i, j])
            lX= ddata[breakprob][breakinfl][0]
            lab = ddata[breakprob][breakinfl][1]
            lrepl = ddata[breakprob][breakinfl][2]
            ax.plot(lX,lab,label="AB produced")
            ax.plot(lX,lrepl,label="Replication")
            # ax.set_title("breakprob = "+breakprob+", breakinfl = "+breakinfl)
            if i == len(ddata)-1:
                ax.set_xlabel("Nr. of growth-promoting genes")
            if j == 0:
                ax.set_ylabel("Frequency")
            # ax.legend(frameon=False)
else:
    # plot condensed data
    import copy
    
    dcondata = copy.deepcopy(ddata)
    for breakprob in ddata:
        for breakinfl in ddata[breakprob]:
            lab = ddata[breakprob][breakinfl][1]
            lrepl = ddata[breakprob][breakinfl][2]
            dcondata[breakprob][breakinfl] = find_median(lrepl) - find_median(lab)
    lcondata= np.zeros((len(dcondata),len(dcondata[breakprob])), dtype=float)
    for i,breakprob in enumerate(sorted(dcondata)):
        for j,breakinfl in enumerate(sorted(dcondata[breakprob])):
            lcondata[i,j] = dcondata[breakprob][breakinfl]
    
    fig = plt.figure()
    ax = plt.gca()
    
    from matplotlib import colors
    
    vmax= 10
    divnorm=colors.DivergingNorm(vmin=0, vcenter=5, vmax=vmax)
    # my_cmap = sns.color_palette("vlag", as_cmap=True)
    # my_cmap = sns.diverging_palette(250, 30, l=65, as_cmap=True)
    my_cmap = sns.color_palette('RdYlBu' , as_cmap=True)

    bla = ax.pcolormesh(lcondata, vmax=vmax,cmap=my_cmap,norm=divnorm)
    
    ax.set_xticks([0.5+i for i in range(4)])
    ax.set_xticklabels([x for x in sorted(dcondata[breakprob])]) # rows are y axis
    ax.set_yticks([0.5+i for i in range(5)])
    ax.set_yticklabels([x for x in sorted(dcondata)])
    cbar = plt.colorbar(bla,extend="max")
    
    cbar.set_label('$\Delta$ growth genes (AB producers, replicating bacteria)')
plt.show()
