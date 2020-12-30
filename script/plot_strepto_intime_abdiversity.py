#!/usr/bin/python2.7


import sys,math,os,string
from subprocess import Popen, PIPE
#from PIL import Image
import matplotlib as mpl
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

genome_pos_in_file = 4

filename = sys.argv[1]
maxtime =sys.maxint 
try :
    maxtime = int(sys.argv[2])
except:
    pass
time = int( (Popen(["head","-n1",filename], stdout=PIPE).stdout.read()).split()[0])
print "Initial time =",time
with open(filename,"r") as fin:
    for line in fin.readlines():
        line =line.split()
        if len(line)<4: continue
        if "n" in line: continue
        
        timenow=int(line[0])
        if timenow != time:
            if time>maxtime: break
            
            ltime.append(time)
            if time%50000==0: print time
            time = timenow
            lav_F.append(np.mean(lF))
            lstd_F.append(np.std(lF))
            lav_A.append(np.mean(lA))
            lstd_A.append(np.std(lA))
            lav_B.append(np.mean(lB))
            lstd_B.append(np.std(lB))

            lF=[];lA=[];lB=[]
            
        try:
            genome=line[genome_pos_in_file]
            lF.append(genome.count("F"))
            lA.append(genome.count("A"))
            lB.append(genome.count("B"))
        except:
            pass

plt.plot(ltime,lav_F,label = "avrg #F")
plt.plot(ltime,lav_A,label = "avrg #A")
plt.plot(ltime,lav_B,label = "avrg #B")

plt.fill_between(ltime, [x-y for x,y  in zip(lav_F,lstd_F)], [x+y for x,y  in zip(lav_F,lstd_F)],alpha=0.5)
plt.fill_between(ltime, [x-y for x,y  in zip(lav_A,lstd_A)], [x+y for x,y  in zip(lav_A,lstd_A)],alpha=0.5)
plt.fill_between(ltime, [x-y for x,y  in zip(lav_B,lstd_B)], [x+y for x,y  in zip(lav_B,lstd_B)],alpha=0.5)
plt.ylim(ymin=0)
plt.legend()
plt.title(filename)
plt.show()
