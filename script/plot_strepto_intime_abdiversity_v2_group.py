#!/usr/bin/python3

"""
plots a multiple data files in time... only medians of course
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

ltime=[]
lF=[]
lA=[]
lB=[]
lav_F=[]
lav_A=[]
lav_B=[]
lstd_F=[]
lstd_A=[]
lstd_B=[]
lstd2_F=[]
lstd2_A=[]
lstd2_B=[]

genome_pos_in_file = 4
fig,ax=plt.subplots(4,1,sharex=True)#,sharey=True)

if sys.argv[1] == '-maxtime': 
    maxtime=int(sys.argv[2])
    start_pos_filename=3
else: 
    maxtime =float("inf") 
    start_pos_filename=1

lfilenames = sys.argv[start_pos_filename:]
for filecounter,filename in enumerate(lfilenames):
    ltime=[]
    lF=[]
    lA=[]
    lB=[]
    lav_F=[]
    lav_A=[]
    lav_B=[]
    time = int( (Popen(["head","-n1",filename], stdout=PIPE).stdout.read()).split()[0])
    print( "Initial time =",time)
    with open(filename,"r") as fin:
        for line in fin:
            line =line.split()
            if len(line)<4: continue
            if "n" in line: continue
            
            timenow=int(line[0])
            if timenow != time:
                if time>maxtime: break
                
                ltime.append(time)
                if time%100000==0: print(time)
                time = timenow
                lav_F.append(np.median(lF))
                # lstd_F.append(np.quantile(lF,0.25))
                # lstd2_F.append(np.quantile(lF,0.75))
                lav_A.append(np.median(lA))
                # lstd_A.append(np.quantile(lA,0.25))
                # lstd2_A.append(np.quantile(lA,0.75))
                lav_B.append(np.median(lB))
                # lstd_B.append(np.quantile(lB,0.25))
                # lstd2_B.append(np.quantile(lB,0.75))

                lF=[];lA=[];lB=[]
                
            try:
                genome=line[genome_pos_in_file]
                lF.append(genome.count("F"))
                lA.append(genome.count("A"))
                lB.append(genome.count("B"))
            except:
                pass
    
    alpha= 0.4+0.6/float(len(lfilenames))
    linewidth = 2. if filecounter==1 else 0.5
    
    ax[0].plot(ltime,lav_F,label = str(filecounter)+" median nr. growth-promoting genes +/- 25\%", alpha=alpha, color='tab:blue')
    ax[1].plot(ltime,lav_B,label = str(filecounter)+" median nt. fragile sites +/- 25\%", alpha=alpha, color='tab:green')
    ax[2].plot(ltime,lav_A,label = str(filecounter)+" median nr. antibiotic genes +/- 25\%", alpha=alpha, color='tab:orange')
    
    ax[3].plot(ltime,lav_F,label = str(filecounter)+" median nr. growth-promoting genes +/- 25\%", alpha=alpha, lw=linewidth, color='tab:blue')
    ax[3].plot(ltime,lav_A,label = str(filecounter)+" median nr. antibiotic genes +/- 25\%", alpha=alpha,lw=linewidth, color='tab:orange')
    ax[3].plot(ltime,lav_B,label = str(filecounter)+" median nt. fragile sites +/- 25\%", alpha=alpha, lw=linewidth,color='tab:green')

ax[-1].set_xlim(xmin=0)
ax[-1].set_xlabel("Time steps")
ax[-1].set_ylabel("Nr. genes")

logscale=True
if logscale:
    plt.yscale('log')
    plt.ylim(ymin=0)
else:
    plt.ylim(ymin=0)

ax[0].legend()
ax[1].legend()
ax[2].legend()
    
# plt.legend()
plt.title('\n'.join(lfilenames))
plt.show()
