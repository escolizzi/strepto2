#!/usr/bin/python3

'''
Makes plots of number of mutants in time. 
'''
import matplotlib as mpl
mpl.use('Qt5Agg')

import sys,math,os,string,random,copy
from subprocess import Popen, PIPE
from collections import Counter

#from PIL import Image
# mpl.use('GTK3Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np


def Gene2Number(gene):
    if gene=='F':return 0
    if gene=='A':return 1
    if gene=='B':return 2

dg={}
count=0
broken_repl=0
filename = sys.argv[1]
# tot_people_alive_pergnm={}
line_number=0
maxtime = 2500 #1250 # 1000#2500

ltime = [0]
l_nr_mutants = [0]
l_nr_tot = [0]
l_anc=[]
avr_nr_F_in_anc=0.0

with open(filename,"r") as fin:
    
    for line in fin:
        line_number+=1
        line = line.split()
        time = int(line[0])
        
        if time > maxtime : break
        
        # if time%100==0: print("Time:",time)
        # print("line_number =", line_number)
        # print("Time: ", time)
        if time not in ltime:
            # if time%100==0: print("Time:",time)
            # print(l_nr_mutants)
            l_nr_mutants.append(l_nr_mutants[-1])
            l_nr_tot.append(l_nr_tot[-1])
            
            ltime.append(time)
            # print(ltime)
            # print(l_nr_mutants)
        if "I:" in line[1] and len(line)==4: 
            try:
                l_anc.append(line[3])
                print(l_anc)
                l1 = [x.count("F") for x in l_anc ]
                print(l1)
                avr_nr_F_in_anc = sum( l1 )/float( len(l_anc) )
                print("avr_nr_F_in_anc: ", avr_nr_F_in_anc)
            except:
                print("Error. line: ", line)
                # exit(1)
        if "R:" in line[1]: 
            l_nr_tot[-1] +=1
            try:
                if line[5] not in l_anc and line[5].count("F")<avr_nr_F_in_anc*1/4. and line[5].count("A")>0:
                    l_nr_mutants[-1] +=1
            except:
                # happens if B is at the beginning of seq
                if line[4]:
                    l_nr_mutants[-1] +=1
                # print("Error R. line: ", line)
                # exit(1)
        # May want to rm this if you want to know mut rate
        if "D:" in line[1]:
            l_nr_tot[-1] -=1
            try:
                if line[3] not in l_anc and line[3].count("F")<avr_nr_F_in_anc/4. and line[3].count("A")>0:
                    l_nr_mutants[-1] -=1
            except:
                # happens if B is at the beginning of seq
                if line[2]:
                    l_nr_mutants[-1] -=1
                # print("Error D. line: ", line)
                # exit(1)
# print(ltime)
# print(l_nr_mutants)
plt.plot( ltime, [ x/float(y) if y!=0 else 0. for x,y in zip(l_nr_mutants,l_nr_tot) ] )
plt.fill_between( ltime, 0, [ x/float(y) if y!=0 else 0. for x,y in zip(l_nr_mutants,l_nr_tot) ], alpha = 0.2 )
plt.title(filename)
plt.show()