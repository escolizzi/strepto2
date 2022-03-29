#!/usr/bin/python3
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

filename = sys.argv[1]
maxtime =float("inf") 
try :
    maxtime = int(sys.argv[2])
except:
    pass
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
            if time%50000==0: print(time)
            time = timenow
            lav_F.append(np.median(lF))
            lstd_F.append(np.quantile(lF,0.25))
            lstd2_F.append(np.quantile(lF,0.75))
            lav_A.append(np.median(lA))
            lstd_A.append(np.quantile(lA,0.25))
            lstd2_A.append(np.quantile(lA,0.75))
            lav_B.append(np.median(lB))
            lstd_B.append(np.quantile(lB,0.25))
            lstd2_B.append(np.quantile(lB,0.75))

            lF=[];lA=[];lB=[]
            
        try:
            genome=line[genome_pos_in_file]
            lF.append(genome.count("F"))
            lA.append(genome.count("A"))
            lB.append(genome.count("B"))
        except:
            pass

plt.plot(ltime,lav_F,label = "median nr. growth-promoting genes +/- 25\%")
plt.plot(ltime,lav_A,label = "median nr. antibiotic genes +/- 25\%")
plt.plot(ltime,lav_B,label = "median nt. fragile sites +/- 25\%")

plt.fill_between(ltime, lstd_F,lstd2_F , alpha=0.5)
plt.fill_between(ltime, lstd_A,lstd2_A,alpha=0.5)
plt.fill_between(ltime, lstd_B,lstd2_B,alpha=0.5)

plt.xlim(xmin=0)
plt.xlabel("Time steps")
plt.ylabel("Nr. genes")

logscale=True
if logscale:
    plt.yscale('log')
    plt.ylim(ymin=0)
else:
    plt.ylim(ymin=0)
    
plt.legend()
plt.title(filename)
plt.show()
