#!/usr/bin/python2.7

'''
Plots the distribution of who produces antibiotics
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

ltime=[]
l_abprod=[]
lav_abprod=[]
ldistr_abprod=[]

lF=[]
lA=[]
lB=[]
lav_F=[]
lav_A=[]
lav_B=[]
lstd_F=[]
lstd_A=[]
lstd_B=[]
tot_nr_gnms_noF=0
tot_nr_gnms=0
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
        # if timenow>1000000: 
        #     print("WARNING: stopping at 500000")
        #     break
        if timenow != time:
            if time>maxtime: break
            
            ltime.append(time)
            if time%50000==0: print time
            time = timenow
            lav_abprod.append(np.mean(l_abprod))
            ldistr_abprod.append(l_abprod)# np.histogram(l_abprod,bins=np.linspace(0.,1.,num=20)) )
            l_abprod=[]
            
            lav_F.append(np.mean(lF))
            lstd_F.append(np.std(lF))
            lav_A.append(np.mean(lA))
            lstd_A.append(np.std(lA))
            lav_B.append(np.mean(lB))
            lstd_B.append(np.std(lB))

            lF=[];lA=[];lB=[]

        try:
            tot_nr_gnms+=1
            genome=line[genome_pos_in_file]
            nF = genome.count("F")
            nA = genome.count("A")
            ab_prod = nA/(float(3.+nA)) * math.exp(-nF)
            if nF==0:
                tot_nr_gnms_noF+=1
                # print("Kijk; ", genome)
            #repl_prod = 0.1*nF/(float(10.+nF))
            #ab_prod = ab_prod/(repl_prod + ab_prod)
            l_abprod.append( ab_prod )

            lF.append(genome.count("F"))
            lA.append(genome.count("A"))
            lB.append(genome.count("B"))
        except:
            pass
print("Tot nr gnms without F: ", tot_nr_gnms_noF, ", out of tot genomes = ", tot_nr_gnms)
fig,ax = plt.subplots(nrows=2,ncols=1)#,sharex=True)
ax[0].plot(ltime,lav_F,label = "avrg #F")
ax[0].plot(ltime,lav_A,label = "avrg #A")
ax[0].plot(ltime,lav_B,label = "avrg #B")
ax[0].set_ylim(ymin=0)
ax[0].legend()

# ax[1].plot(ltime,lav_abprod,label = "avrg ABprod")

#mybins=np.linspace(0.,1.02,103)#np.linspace(0.,1.,51)
mybins=np.linspace(0.,1.2,7)
print(mybins)

for thistime,thisdistr_abprod in zip(ltime[::13],ldistr_abprod[::13]):
    print("hisogramming time: ", thistime)
    hist,bins = np.histogram(thisdistr_abprod, bins= mybins)
    ax[1].plot(bins[:-1],hist, label='time = '+str(thistime))
ax[1].legend()
# ax[1].set_yscale('log')
# ax[1].violinplot(ldistr_abprod[::5], ltime[::5], widths=ltime[2]-ltime[0], showextrema=False, points=50)
print("t1-t0 = ",ltime[1]-ltime[0])

plt.title(filename)
plt.show()
