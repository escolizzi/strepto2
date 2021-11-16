#!/usr/bin/python3.8

'''
Adapted from the competition plotting of the MC experiments.
Plot the populationsize of both groups over time
'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
#mpl.use('QT5Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import seaborn as sns
from matplotlib.lines import Line2D
from collections import OrderedDict
import numpy as np

#manager = plt.get_current_fig_manager()
#print manager
#sys.exit(1)
#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################



colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""

#fig, (ax0,ax1) = plt.subplots(nrows=2, sharex=True)
fig, ax = plt.subplots(nrows=1)

if len(sys.argv) <4:
    print ("This is the program 'plot_competition_popsize_dist.py'")
    print ("Usage: ./plot_division_time.py <figure name> <season duration> <nr replicates> <datafile 1> <datafile 2> <...>")
    sys.exit(1)
else:
    figname=sys.argv[1]
    season=int(sys.argv[2])
    nrreplicates=int(sys.argv[3])
redcols=["firebrick", "orangered", "lightcoral"]
blucols=["mediumblue", "cornflowerblue", "skyblue"]
filenr=0

heatmapdata=np.zeros((2*5,),dtype=float)

for files in sys.argv[4:]:
    
    popsize=[ 0 for x in range(2)]
    print(files)
    print (filenr) 
    filetype=0

    with open(files,"r") as fin:

        #read first line first to get the first time point (could probably more clever, but hey)
        seasoncount=0
        prevtime=0
        registered=0
        index=int(filenr/nrreplicates)
        print (index)

        for line in reversed(fin.readlines()):
            
            line1=line.split(' ')
            if(len(line1)>3 or ',' in line):
                time = int(line1[0])
                nrab=int(line1[1])
                sporetype=int(line1[3])
                #count population size
                if sporetype:
                    popsize[sporetype-1]+=1

            else:
                filetype=1
                if (int(line1[1])>int(line1[2])):
                    heatmapdata[index]+=1./float(nrreplicates)
                    print("Evolved won")
                elif (int(line1[1])<int(line1[2])):
                    heatmapdata[index]-=1./float(nrreplicates)
                    print("Competitor won")
                break

          #  if time!=prevtime:
        if(not filetype):
            if (popsize[0]>popsize[1]):
                heatmapdata[index]+=1./float(nrreplicates)
                print("Evolved won")
            elif (popsize[0]<popsize[1]):
                heatmapdata[index]-=1./float(nrreplicates)
                print("Competitor won")

    filenr+=1

#display consensus image
imcons=np.reshape(heatmapdata,(2,5))
#fig=plt.figure() #!
#fig.set_size_inches(1, 1, forward = False) #!
#ax = plt.Axes(fig, [0., 0., 1., 1.]) #!
##ax.set_axis_off() #!
#fig.add_axes(ax) #!
x=[0,1,2,3,4]
xlabels=["1", "2", "3", "4", "5"]
ylabels=["0.075", "0.1"]
ax.set_xlabel("Kab (10^-5)")
ax.set_ylabel("Krep")
plt.xticks(x, xlabels)
plt.yticks(x, ylabels)
im=plt.imshow(imcons, cmap=plt.cm.bwr, origin='lower')
plt.clim(-1, 1)
plt.colorbar(im, ticks=[-1, -0.5, 0, 0.5, 1]) # adding the colorbar on the right
plt.savefig(figname, bbox_inches='tight')
#plt.close()
#plt.show()