#!/usr/bin/python3.8

'''
Adapted from the competition plotting of the MC experiments.
Plot the populationsize of both groups over time
'''

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
fig, ax0 = plt.subplots(nrows=1)

if len(sys.argv) <4:
    print ("This is the program 'plot_competition_popsize_dist.py'")
    print ("Usage: ./plot_division_time.py <figure name> <season duration> <datafile 1> <datafile 2> <...>")
    sys.exit(1)
else:
    figname=sys.argv[1]
    season=int(sys.argv[2])

redcols=["firebrick", "orangered", "lightcoral"]
blucols=["mediumblue", "cornflowerblue", "skyblue"]
filenr=0

for files in sys.argv[3:]:
    
    popsize=[ [0] for x in range(2)]
    times=[0]

    print(files)

    with open(files,"r") as fin:

        #read first line first to get the first time point (could probably more clever, but hey)
        seasoncount=0
        prevtime=0
        registered=1

        for line in fin:
            #print(line)
            line=line.split(' ')
            time = int(line[0])
            nrab=int(line[1])
            sporetype=int(line[3])
                
            if time!=prevtime:
                popsize[0].append(0)
                popsize[1].append(0)
                prevtime=time
                times.append(time)
                registered=0
            
            #count population size
            if sporetype:
                popsize[sporetype-1][-1]+=1


            #register end of season
            if not(prevtime%season) and not registered:
                seasoncount+=1
                registered=1
                print("season: ",seasoncount)

    if filenr==0:        
        ax0.plot(times, popsize[0], c=redcols[filenr], label="evolved", linewidth=2)
        ax0.plot(times, popsize[1], c=blucols[filenr], label="competitor", linewidth=2)
        ax0.legend()
    
    else:
        ax0.plot(times, popsize[0], c=redcols[filenr], linewidth=2)
        ax0.plot(times, popsize[1], c=blucols[filenr], linewidth=2)
        ax0.legend()
        #ax1.plot(times, dist[0], c=redcols[filenr], linewidth=1.5)
        #ax1.fill_between(times, bottom[0], top[0], color=redcols[filenr],alpha=0.25,linewidth=0.0)

        #ax1.plot(times, dist[1], c=blucols[filenr], linewidth=1.5)
        #ax1.fill_between(times, bottom[1], top[1], color=blucols[filenr],alpha=0.25,linewidth=0.0)
    filenr+=1

#ax0.set_ylim(0,season)
#ax0.set_xlabel('time (MCS)')
ax0.set_ylabel('group size')
ax0.set_xlabel('time')
ax0.set_title(figname)
plt.show()
