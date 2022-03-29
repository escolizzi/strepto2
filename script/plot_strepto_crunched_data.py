#!/usr/bin/python3
from ctypes.util import find_library
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

ldata=[]

with open(filename,"r") as fin:
    for line in fin:
        line =[int(x) if not '.' in x else float(x) for x in line.split()]
        ldata.append(line)
ldata=[*zip(*ldata)]
ltime,lstd_F,lav_F,lstd2_F,lstd_A,lav_A,lstd2_A,lstd_B,lav_B,lstd2_B=ldata

plt.plot(ltime,lav_F,label = "median nr. growth-promoting genes +/- 25\%")
plt.plot(ltime,lav_A,label = "median nr. antibiotic genes +/- 25\%")
plt.plot(ltime,lav_B,label = "median nt. fragile sites +/- 25\%")

plt.fill_between(ltime, lstd_F,lstd2_F , alpha=0.5)
plt.fill_between(ltime, lstd_A,lstd2_A,alpha=0.5)
plt.fill_between(ltime, lstd_B,lstd2_B,alpha=0.5)

plt.xlim(xmin=0)
plt.xlabel("Time steps")
plt.ylabel("Nr. genes")

logscale=False
if logscale:
    plt.yscale('log')
    plt.ylim(ymin=0.1)
else:
    plt.ylim(ymin=0)
    # plt.ylim(ymax=20)
    
plt.legend()
plt.title(filename)
plt.show()
