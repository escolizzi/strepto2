#!/usr/bin/python3

'''
This is to check how I am actually initialising the genomes for Fig. SF2.
'''

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
from numpy.core.fromnumeric import prod

lgenome=[]
for filename in sys.argv[1:]:
    print("Opening: ",filename)
    ltime=[]
    lab_prod=[]

    with open(filename,"r") as fin: 
        for line in fin:
            line = line.split()
            genome = line[4]
            if line[4] =='n':continue
            if len(genome) != 20: 
                print("Error - what's wrong with this genome? Len =",len(genome))
                print(genome)
                continue
            lgenome.append(genome) 
                        
lgenome_T = [*zip(*lgenome)]
lcounts = [ [bla.count('F')/len(bla),bla.count('A')/len(bla),bla.count('B')/len(bla)] for bla in lgenome_T ]
lcounts_T = [*zip(*lcounts)]
plt.fill_between(range(20), 0, lcounts_T[0])

plt.fill_between(range(20), lcounts_T[0], [x+y for x,y in zip(lcounts_T[0],lcounts_T[1])])

plt.fill_between(range(20), [x+y for x,y in zip(lcounts_T[0],lcounts_T[1])], [x+y+z for x,y,z in zip(lcounts_T[0],lcounts_T[1],lcounts_T[2])])

print("mean F: ",np.mean(lcounts_T[0]),", mean A: ",np.mean(lcounts_T[1]),", mean B: ",np.mean(lcounts_T[2]))
print("mean F+A+B: ",np.mean(lcounts_T[0]) + np.mean(lcounts_T[1]) + np.mean(lcounts_T[2]))

# plt.fill_between(range(20), 0, lcounts_T[0])

plt.show()
