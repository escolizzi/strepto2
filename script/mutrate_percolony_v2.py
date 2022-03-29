#!/usr/bin/python3

import sys,math
import matplotlib.pyplot as plt
import numpy as np

dc={}
tot_lines=0
tot_linesnoF=0
if len(sys.argv)<2:
  print("mutrate_percolony.py, calculates and plots the fraction of AB mutants *PER COLONY* in the system, feed it only one or few time steps")
filename=sys.argv[1]

l_abpr=[]

#filename="data_strpt_season1000_maxpAB0.01_betatradeoff3_bitstrABlen12_frnonproducing0._nomix_btwn_seasns_muttypeC.txt_t1000000.txt"
with open(filename, "r") as fin:
  for line in fin.readlines():
    line=line.split()
    if len(line)<5: continue
    if line[4]=='n': continue
    time=line[0]
    tot_lines+=1
    
    cln=line[3]
    #print "line:",line
    #print "seq:",seq
    seq=line[4]
    noF=0
    
    F=seq.count("F")
    A=seq.count("A")
    abpr=A/(A+3.)*(math.exp(-1.*F));
    growth=0.1*F/(F+10.)
    l_abpr.append(abpr)
    if abpr>0.1 and growth<0.01:
    #if seq.count("F")<=2:
      # print( abpr)
      noF=1
      tot_linesnoF+=1
      #print "noF true", seq
    key=time+'_'+cln
    if key in dc: 
      dc[key][0]+=1
      dc[key][1]+=noF
    else:
      dc[key]=[1,noF]

l_mutfr = [ x[1]/float(x[0]) for key, x in dc.items() ]
# step = 0.0005
step = 0.002
maxhist=0.15
hist,bins= np.histogram(l_mutfr, bins=np.arange(0., maxhist, step=step))
plt.plot(bins[:-1],hist,lw=2,label="distribution of AB mutants\nper colony")

# hist2,bins2= np.histogram(l_abpr, bins=np.arange(0., 1., step=step))
# plt.plot(bins2[10:-1],hist2[10:],lw=2,label="distribution of AB prod total")

print ("tot_lines",tot_lines,"tot_linesnoF",tot_linesnoF)
averagemutrate=tot_linesnoF/float(tot_lines)

plt.plot([averagemutrate,averagemutrate],[0,0.667*max(hist)],lw=2,linestyle='dashed',label="avrg. fraction of\nAB mutants")
plt.xlabel("fraction AB producing mutant per colony")
plt.ylabel("colonies")
plt.legend()
#print "No mut?"
#for key,x in dc.items():
#  if x[1]/float(x[0])==0: 
#    print key, x[0]

plt.show()
