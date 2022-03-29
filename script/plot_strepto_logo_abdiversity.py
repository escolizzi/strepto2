#!/usr/bin/python3

'''
Reads a data file (best is one time step - but the script does not care)
Makes a bunch of plots that show genome structure
'''

import matplotlib as mpl
mpl.use('qt5agg')
import sys,math,os,string
from subprocess import Popen, PIPE
#from PIL import Image
import seaborn as sns

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
import numpy as np

def Gene2Number(gene):
    if gene=='F':return 0
    if gene=='A':return 1
    if gene=='B':return 2
    if gene=='H':return 3


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

genome_pos_in_file = 4

with open(filename,"r") as fin:
    for line in fin.readlines():
        line =line.split()
        if len(line)<5 or line[genome_pos_in_file]=='n' or line[genome_pos_in_file]==',': continue
        if len(  set(line[genome_pos_in_file]).difference( set(['F','A','B','H'])) )!=0:
            print("Something wrong here")
            print(line[genome_pos_in_file])
            
        try:
          lgenome.append(line[genome_pos_in_file])
          if maxlen<len(line[genome_pos_in_file]): maxlen = len(line[genome_pos_in_file])
        except:
          pass
          #print line
          #sys.exit(1)  

tot_letter_used = len(set().union(*lgenome) )
print("tot_letter_used", tot_letter_used)
print("Letters: ", set().union(*lgenome))

for pos in range(maxlen):
    genome_logo.append(tot_letter_used*[0])
    
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
    genome_logo2.append(tot_letter_used*[0])

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

# l_nF_avpos=[]
# for genome in lgenome:
#     len_genome = float(len(genome))
#     n_F = 0
#     av_pos=0.
#     for i,x in enumerate(genome):
#         if x == 'F': 
#             n_F +=1
#             av_pos += i/len_genome # pos goes from 0 to 1
#     if n_F>0:
#         av_pos = av_pos / float(n_F)
#     else:
#         if av_pos > 0.:
#             print("Error: n_f = 0 but av_pos>0?")
#             print("Genome: ", genome)
#     l_nF_avpos.append([ len_genome , av_pos ])
# l_nF_avpos = map(list, zip(*l_nF_avpos))

#Oh position in the genome, what's your frequency of fitness genes?
bin_gnm=25
print("Maxlen=",maxlen)
l_gnmsize_posF=np.zeros(shape=(bin_gnm,maxlen+1), dtype=float)
l_gnmsize_posA=np.zeros(shape=(bin_gnm,maxlen+1), dtype=float)
l_gnmsize_posB=np.zeros(shape=(bin_gnm,maxlen+1), dtype=float)
if tot_letter_used==4:
    l_gnmsize_posH=np.zeros(shape=(bin_gnm,maxlen+1), dtype=float)
# print("numpy emptyF is:",l_gnmsize_posF)
# print("numpy emptyA is:",l_gnmsize_posA)
# print("numpy emptyB is:",l_gnmsize_posB)
if tot_letter_used==4:
    l_pos_10th_H = []

for ignm,genome in enumerate(lgenome):
    # l_thisgnm_posF=bin_gnm*[0.]
    len_genome = len(genome)
    # if len_genome>60:continue
    flen_genome = float(len_genome)
    n_F = float(genome.count('F'))
    n_A = float(genome.count('A'))
    n_B = float(genome.count('B'))
    if tot_letter_used==4:
        n_H = float(genome.count('H'))
        
    if n_F>0:
        for i,x in enumerate(genome):
            if x == 'F': 
                try:
                    l_gnmsize_posF[ int(bin_gnm*(i/(flen_genome))) , len_genome ] += 1./n_F
                except:
                    print("Error?")
                    print(int(bin_gnm*(i/(flen_genome))) )
                    print(len_genome)
                    sys.exit(1)
    if n_A>0:
        for i,x in enumerate(genome):
            if x == 'A': 
                try:
                    l_gnmsize_posA[ int(bin_gnm*(i/(flen_genome))) , len_genome ] += 1./n_A
                except:
                    print("Error?")
                    print(int(bin_gnm*(i/(flen_genome))) )
                    print(len_genome)
                    sys.exit(1)
    if n_B>0:
        # print("Genome: ",genome, ", nr. B: ", n_B)
        # av_n_B += n_B
        for i,x in enumerate(genome):
            if x == 'B': 
                try:
                    #print("Before:",l_gnmsize_posB[ int(bin_gnm*(i/(flen_genome))) , len_genome ])
                    l_gnmsize_posB[ int(bin_gnm*(i/(flen_genome))) , len_genome ] += 1./n_B
                    #print("After:",l_gnmsize_posB[ int(bin_gnm*(i/(flen_genome))) , len_genome ])
                    # if l_gnmsize_posB[ int(bin_gnm*(i/(flen_genome))) , len_genome ] > 10000:
                    #     print(l_gnmsize_posB[ int(bin_gnm*(i/(flen_genome))) , len_genome ])
                    #     print("genome nr ",ignm,":", genome)
                    #     print("nB=",n_B, "1/n_B=",1./n_B,"i:",i,"int(bin_gnm*(i/(flen_genome))): ",int(bin_gnm*(i/(flen_genome))))
                    #     sys.exit(1)
                except:
                    print("Error?")
                    print(int(bin_gnm*(i/(flen_genome))) )
                    print(len_genome)
                    sys.exit(1)
    if tot_letter_used==4:
        if n_H>0:
            for i,x in enumerate(genome):
                if x == 'H': 
                    try:
                        l_gnmsize_posH[ int(bin_gnm*(i/(flen_genome))) , len_genome ] += 1./n_H
                    except:
                        print("Error?")
                        print(int(bin_gnm*(i/(flen_genome))) )
                        print(len_genome)
        if n_H>=10:
            l_pos_10th_H.append( [Hindex for Hindex,gene in enumerate(genome) if gene=='H' ][9] / float(len(genome)) ) # index of 10th occurrence from the 5' of gene H
#Find min adn max, this is so that we can share a color map across them
min_all_array = np.min( [np.min(l_gnmsize_posF), np.min(l_gnmsize_posA), np.min(l_gnmsize_posB)] )
max_all_array = np.max( [np.max(l_gnmsize_posF), np.max(l_gnmsize_posA), np.max(l_gnmsize_posB)] )

if tot_letter_used==4:
    min_all_array = np.min([min_all_array, np.min(l_gnmsize_posH)])
    max_all_array = np.max([max_all_array, np.max(l_gnmsize_posH)])

# fig, ax = plt.subplots(ncols=1,nrows=3)
fig = plt.figure(constrained_layout=True)
spec2 = gridspec.GridSpec(nrows=5, ncols=3*tot_letter_used+1 , figure=fig)
ax10 = fig.add_subplot(spec2[1:4 , 0:3])
ax00 = fig.add_subplot(spec2[0   , 0:3], sharex=ax10)
ax20 = fig.add_subplot(spec2[4   , 0:3])
# ax_cbarF = fig.add_subplot(spec2[1:4, -1])

ax11 = fig.add_subplot(spec2[1:4 , 3:6], sharey = ax10)
ax01 = fig.add_subplot(spec2[0   , 3:6], sharex = ax11, sharey = ax00) 
ax21 = fig.add_subplot(spec2[4   , 3:6])
# ax_cbarA = fig.add_subplot(spec2[5:8, -1])

ax12 = fig.add_subplot(spec2[1:4 , 6:9], sharey = ax10)
ax02 = fig.add_subplot(spec2[0   , 6:9], sharex = ax12, sharey = ax00)
ax22 = fig.add_subplot(spec2[4   , 6:9])
# ax_cbarB = fig.add_subplot(spec2[9:, -1])

if tot_letter_used==4:
    ax13 = fig.add_subplot(spec2[1:4 , 9:12], sharey = ax10)
    ax03 = fig.add_subplot(spec2[0   , 9:12], sharex = ax12, sharey = ax00)
    ax23 = fig.add_subplot(spec2[4   , 9:12])

ax1last = fig.add_subplot(spec2[1:4 , -1], sharey = ax10)

# my_cmap = sns.color_palette('afmhot_r', as_cmap=True)
# my_cmap = ListedColormap(sns.color_palette('Blues').as_hex())
my_cmapF = sns.color_palette("gist_earth_r", as_cmap=True)
# F genes
X, Y = np.meshgrid( range(maxlen+1) , np.linspace(0.,1.,num=bin_gnm))
F_2dplot = ax10.pcolormesh(Y,X, l_gnmsize_posF, cmap = my_cmapF, vmin=min_all_array, vmax=max_all_array )
fig.colorbar(F_2dplot, ax=ax20, orientation='horizontal')

l_gnmsize_posF_sumgnmsize = [ sum(x) for x in  l_gnmsize_posF ]
ax00.plot(np.linspace(0,1,bin_gnm), [x/float(sum(l_gnmsize_posF_sumgnmsize)) for x in l_gnmsize_posF_sumgnmsize], c='dodgerblue', label = 'Distribution of\ngrowth genes', lw=4)
ax00.fill_between(np.linspace(0,1,bin_gnm), [x/float(sum(l_gnmsize_posF_sumgnmsize)) for x in l_gnmsize_posF_sumgnmsize], color='aliceblue')
ax00.legend()
ax00.margins(0,0)   
ax10.margins(0,0) 

# A genes
my_cmapA = sns.color_palette("gist_earth_r", as_cmap=True)
A_2dplot = ax11.pcolormesh(Y,X, l_gnmsize_posA, cmap = my_cmapA, vmin=min_all_array, vmax=max_all_array)
fig.colorbar(A_2dplot, ax=ax21, orientation='horizontal')

l_gnmsize_posA_sumgnmsize = [ sum(x) for x in  l_gnmsize_posA ]
ax01.plot(np.linspace(0,1,bin_gnm), [x/float(sum(l_gnmsize_posA_sumgnmsize)) for x in l_gnmsize_posA_sumgnmsize], c='darkorange', label = 'Distribution of\nantibiot. genes', lw=4)
ax01.fill_between(np.linspace(0,1,bin_gnm), [x/float(sum(l_gnmsize_posA_sumgnmsize)) for x in l_gnmsize_posA_sumgnmsize], color='moccasin')
ax01.legend()
ax11.margins(0,0)
ax01.margins(0,0)

# B genes
my_cmapB = sns.color_palette("gist_earth_r", as_cmap=True)
B_2dplot = ax12.pcolormesh(Y,X, l_gnmsize_posB, cmap = my_cmapB, vmin=min_all_array, vmax=max_all_array)
fig.colorbar(B_2dplot, ax=ax22, orientation='horizontal')

l_gnmsize_posB_sumgnmsize = [ sum(x) for x in  l_gnmsize_posB ]
ax02.plot(np.linspace(0,1,bin_gnm), [x/float(sum(l_gnmsize_posB_sumgnmsize)) for x in l_gnmsize_posB_sumgnmsize],c='seagreen', label = 'Distribution of\nfragile sites', lw=4)
ax02.fill_between(np.linspace(0,1,bin_gnm), [x/float(sum(l_gnmsize_posB_sumgnmsize)) for x in l_gnmsize_posB_sumgnmsize], color='honeydew')
ax02.legend()
ax02.margins(0,0)
ax12.margins(0,0)

# H genes
if tot_letter_used==4: 
    my_cmapH = sns.color_palette("gist_earth_r", as_cmap=True)
    H_2dplot = ax13.pcolormesh(Y,X, l_gnmsize_posH, cmap = my_cmapH, vmin=min_all_array, vmax=max_all_array)
    fig.colorbar(H_2dplot, ax=ax22, orientation='horizontal')
    
    l_gnmsize_posH_sumgnmsize = [ sum(x) for x in  l_gnmsize_posH ]
    # ax03.plot(np.linspace(0,1,bin_gnm), [x/float(sum(l_gnmsize_posH_sumgnmsize)) for x in l_gnmsize_posH_sumgnmsize],c='seagreen', label = 'Distribution of\nhomeost. genes', lw=4)
    # ax03.fill_between(np.linspace(0,1,bin_gnm), [x/float(sum(l_gnmsize_posH_sumgnmsize)) for x in l_gnmsize_posH_sumgnmsize], color='honeydew')
    
    
    Hhist,Hbins = np.histogram( l_pos_10th_H , bins = np.linspace(0,1,bin_gnm), normed=True )
    scale_for_match_with_other_plots=0.02
    ax03.plot(Hbins[:-1], [scale_for_match_with_other_plots*bla for bla in Hhist], c='darkgoldenrod', label = 'Distribution of\n10th homeost. genes', lw=4)
    
    ax03.legend()
    
    ax03.margins(0,0)
    ax13.margins(0,0)

# Distribution of genome length
hist, bins = np.histogram( lgenome_length, bins=range(0, 1+max(lgenome_length)), density=True )
ax1last.plot( hist, bins[:-1], c='darkgrey', lw=4)
ax1last.fill_betweenx(bins[:-1],hist, color='whitesmoke')
ax1last.set_ylim([0.,max(bins[:-1])])
ax1last.margins(0,0)
# ax13.legend()
plt.show()
