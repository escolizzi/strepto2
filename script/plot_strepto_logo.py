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

def Gene2Number(gene):
    if gene=='F':return 0
    if gene=='A':return 1
    if gene=='B':return 2


lgenome=[]

lF=[]
lA=[]
lB=[]
lav_F=[]
lav_A=[]
lav_B=[]
lstd_F=[]
lstd_A=[]
lstd_B=[]

filename = sys.argv[1]
genome_logo=[]
l_totpos=[]
l_totpos2=[]
max_totpos=0
maxlen=0
lgenome_length=[]
genome_logo2=[]
totpos2=0
max_totpos2=0

with open(filename,"r") as fin:
    for line in fin.readlines():
        line =line.split()
        lgenome.append(line[3])
        if maxlen<len(line[3]): maxlen = len(line[3])

for pos in range(maxlen):
    genome_logo.append([0,0,0])
    
    for genome in lgenome:
        # if len(genome)<10: continue # This has no visual effect on the final logo 
        lgenome_length.append(len(genome))
        try:
            gene = genome[pos]
            genome_logo[-1][Gene2Number(gene)]+=1
                            
        except:
            pass
    tot_thispos = sum(genome_logo[-1])
    max_totpos = tot_thispos if tot_thispos > max_totpos else max_totpos
    l_totpos.append(tot_thispos)
    genome_logo[-1] = [x/float(tot_thispos) for x in genome_logo[-1]]

genome_logo = map(list, zip(*genome_logo))

for pos in range(11):
    genome_logo2.append([0,0,0])

    for genome in lgenome:
        if len(genome) > 10: continue
        try:
            gene = genome[pos]
            genome_logo2[-1][Gene2Number(gene)]+=1
                            
        except:
            pass
    tot_thispos2 = sum(genome_logo2[-1])
    max_totpos2 = tot_thispos2 if tot_thispos2 > max_totpos2 else max_totpos2
    # l_totpos2.append(tot_thispos2)
    genome_logo2[-1] = [x/float(tot_thispos) for x in genome_logo2[-1]]

genome_logo2 = map(list, zip(*genome_logo2))

fig, ax = plt.subplots(ncols=1,nrows=3,sharex=True)
ax[0].plot(range(maxlen), genome_logo[0], label="F")
ax[0].plot(range(maxlen), genome_logo[1], label="A")
ax[0].plot(range(maxlen), genome_logo[2], label="B")

ax[0].plot(range(maxlen), [x/float(max_totpos) for x in l_totpos], lw=0.5, c='k', ls='--', label='support')
# ax[0].title(filename)
ax[0].legend()

# plt.hist()

ax[1].hist(lgenome_length,  bins=range(0, 1+max(lgenome_length)) )

ax[2].plot(range(11), genome_logo2[0], label="F")
ax[2].plot(range(11), genome_logo2[1], label="A")
ax[2].plot(range(11), genome_logo2[2], label="B")

plt.margins(0,0)
plt.show()
