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
lreg=[]
lstd_reg=[]
lav_reg=[]

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
        timenow=int(line[0])
        if timenow != time:
            if time>maxtime: break
            
            ltime.append(time)
            if time%50000==0: print time
            time = timenow
            
            lstd_reg.append(np.std(lreg))
            lav_reg.append(np.mean(lreg))
            
            lreg=[]
        try:
            lreg.append(float(line[-1]))
            
        except:
            pass

plt.plot(ltime,lav_reg,label = "avrg ref")

plt.fill_between(ltime, [x-y for x,y  in zip(lav_reg,lstd_reg)], [x+y for x,y  in zip(lav_reg,lstd_reg)],alpha=0.5)
plt.ylim(ymin=0)
plt.legend()
plt.title(filename)
plt.show()
