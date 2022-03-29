#!/usr/bin/python3

'''
Plots the distribution of who produces antibiotics
'''

import sys,math,os,string
from subprocess import Popen, PIPE
#from PIL import Image
import matplotlib as mpl
mpl.use('QT5Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np
from numpy.core.fromnumeric import prod

colors = ['darkgoldenrod', 'royalblue','crimson','seagreen','darkgrey','chartreuse']

for filenumber,filename in enumerate(sys.argv[1:]):
    print("Opening: ",filename)
    ltime=[]
    lab_prod=[]
    lpopsize=[]
    if "popsize" in filename:
        
        with open(filename,"r") as fin:
            for line in fin:
                line = line.split()
                time=int(line[0])
                popsize = float(line[1])
                ab_prod = float(line[2])#*317.
                
                ltime.append(time)
                lpopsize.append(popsize)
                # lab_prod.append(ab_prod/popsize)
                lab_prod.append(ab_prod)
        
        # l_deltaab_prod = [  (y-x)/(z*(1-z/90000.)) for x,y,z in zip(lab_prod[:-1],lab_prod[1:],lpopsize[:-1])  ]
        # lab_sliced = [ l_deltaab_prod[i:i+250] for i in range(0,len(lab_prod),250) ]
        # for lab_slice in lab_sliced[2:-1]:
        #     # print (lab_slice)
        #     N=10
        #     plt.plot([x for x in range(0,2500,10)][1:],np.convolve(lab_slice[1:],np.ones(N)/N, mode='valid'),lw=0.2,alpha = 0.7,c=colors[filenumber])
        
        # plt.hlines(0, xmin = 0, xmax=2500, lw=0.5,alpha = 0.7,colors='black')
        lab_prod = [ab_prod/(popsize*317.) for ab_prod,popsize in zip(lab_prod,lpopsize)]
        lab_sliced = [ lab_prod[i:i+250]  for i in range(0,len(lab_prod),250) ]
        #remove initial conditions because they are not representative
        # and last because there is only one data point before strepto breaks
        for lab_slice in lab_sliced[2:-1]:
            # print (lab_slice)
            plt.plot([x for x in range(0,2500,10)],lab_slice,lw=0.1,alpha = 0.3,c=colors[filenumber])
        
        # for lab_slice in lab_sliced[:1]:
        #     # print (lab_slice)
        #     plt.plot([x for x in range(0,2500,10)],lab_slice,lw=1.3,alpha = 0.99,c=colors[filenumber])
        # mean_abprod=np.mean(lab_prod[250:])
        # std_abprod=np.std(lab_prod[250:])
        # print("Mean: ", mean_abprod, ", std: ",std_abprod)
        
        # for bla in lab_sliced:
        #     print(len(bla))
        # sys.exit(1)
        
        lab_sliced_T = np.transpose(lab_sliced[2:-1])
        # print("len(lab_sliced_T)",len(lab_sliced_T))
        # print("len(lab_sliced[0])",len(lab_sliced[0]))
        lab_sliced_T_mean_per_point = [np.mean(bla) for bla in lab_sliced_T]
        # plt.plot([x for x in range(0,2500,10)],lab_sliced_T_mean_per_point,lw=0.5,alpha = 1.,c=colors[filenumber])
        mean_abprod=np.mean(lab_sliced_T_mean_per_point)
        std_abprod=np.std(lab_sliced_T_mean_per_point)
        print("Mean: ", mean_abprod, ", std: ",std_abprod)
        #plot mean and std
        plt.plot(2*[500*filenumber] , [mean_abprod-std_abprod,mean_abprod+std_abprod],c=colors[filenumber])
        plt.scatter([500*filenumber],[mean_abprod],c=colors[filenumber])
        
    # THIS MATCHES DATA PERFECTLY - UN-NEEDED
    # elif "avr" in filename:
    #     with open(filename,"r") as fin:
    #         for line in fin:
    #             line = line.split()
    #             ltime.append(int(line[0]))
    #             lab_prod.append(float(line[2]))
    #         # plt.plot(ltime,lab_prod)
    #         # plt.show()
    #         # sys.exit(1)
    #         lab_sliced = [ lab_prod[i:i+5]  for i in range(1,len(lab_prod),5) ]
    #         print("lab_sliced:",lab_sliced)
    #         for lab_slice in lab_sliced[3:-1]:
    #             plt.plot([x for x in range(500,2500+500,500)],lab_slice,linestyle="",marker = "o",alpha = 0.7,c=colors[filenumber])

# plt.vlines(range(0,ltime[-1],2500), ymin=0,ymax=0.1/317., linestyles='dashed', lw = 0.5,color = 'red')
# plt.yscale('log')
# plt.xticks(range(0,ltime[-1]+2500,500))
plt.xticks(range(0,2500+250,250))
plt.xlabel("time (within growth cycle)")
plt.ylabel("Average AB production")
plt.show()
