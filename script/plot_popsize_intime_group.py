#!/usr/bin/python3

"""
plots popsize of multiple data files in time
"""

import matplotlib as mpl
mpl.use('qt5Agg')

import sys,math,os,string
from subprocess import Popen, PIPE
#from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np

genome_pos_in_file = 4
fig,ax=plt.subplots(3,1,sharex=True)#,sharey=True)

if sys.argv[1] == '-maxtime': 
    maxtime=int(sys.argv[2])
    start_pos_filename=3
else: 
    maxtime =float("inf") 
    start_pos_filename=1

lfilenames = sys.argv[start_pos_filename:]
for filecounter,filename in enumerate(lfilenames):
    ltime=[]
    lpop=[]
    pop=0
    print("Opening file: ",filename)
    time = int( (Popen(["head","-n1",filename], stdout=PIPE).stdout.read()).split()[0])
    print( "Initial time =",time)
    with open(filename,"r") as fin:
        for line in fin:
            line =line.split()
            if len(line)<4: continue
            if "n" in line: continue
                        
            timenow=int(line[0])
            if timenow != time:
                if time%250000==0: print("Time: ", time)
                if time>maxtime: break
                time = timenow
                ltime.append(time)
                lpop.append(pop)
                pop=0
                
            pop+=1
    
    alpha= 0.4+0.6/float(len(lfilenames))
    linewidth = 2. if filecounter==1 else 0.5
    
    ax[0].plot(ltime,lpop,label = str(filecounter)+" population size (end of growth cycle)",lw=linewidth, alpha=alpha, color='tab:grey')
    
ax[0].set_xlim(xmin=0)
ax[0].set_xlabel("Time steps")
ax[0].set_ylabel("Population size")

logscale=False
if logscale:
    plt.yscale('log')
    plt.ylim(ymin=0)
else:
    plt.ylim(ymin=0)

ax[0].legend()
    
# plt.legend()
plt.title('\n'.join(lfilenames))
plt.show()
