#!/usr/bin/python3

'''
Makes plots of genome and genes distribution.
Puts on the x axis a parameter of your choice. 
Puts on a mini axis within the x axis a second parameter of your choice. 
The parameters you want to use must be in the file name, and you have to specify its position or the subtring

Usage: 
./plot_genomelendistr_2parameters.py [pos1] [pos2] [list of files]

'''

import matplotlib as mpl
from numpy.core.fromnumeric import _all_dispatcher
mpl.use('Qt5Agg')

import sys,math,os,string
from subprocess import Popen, PIPE
from collections import Counter

#from PIL import Image
# import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import seaborn as sns

from mpl_toolkits.mplot3d import Axes3D
import numpy as np


lgenome=[]
llgenome=[]
lgenome_noFlen=[]
llgenome_noFlen=[]

llF=[]
llA=[]
llB=[]
llABfield=[]

lav_F=[]
lav_A=[]
lav_B=[]
lstd_F=[]
lstd_A=[]
lstd_B=[]
lABfield=[]

# filename = sys.argv[1]
genome_logo=[]
l_totpos=[]
l_totpos2=[]
max_totpos=0
maxlen=0
lgenome_length=[]
llgenome_length=[]

ll_howmany_ab=[]
ll_susceptibility=[]

genome_logo2=[]
totpos2=0
max_totpos2=0

nr_ab_pos_in_file = 1
ab_list_pos_in_file = 2
genome_pos_in_file = 4
gene_ab_pos_in_file = 5

dpos_violinplot={}

look_for_matching_string1 = False
look_for_matching_string2 = False

calculate_susceptibility=False
#  Wishful thinking without better input management
# if "-nosusceptbility" in sys.argv: 
#     calculate_susceptibility==False
if calculate_susceptibility==True: 
    print("WARNING, NOT CALCULATING SUSCEPTIBILITY")
    sys.exit(1)

try:
    where_to_look_in_filename1 = int(sys.argv[2])
except ValueError:
    #maybe it is a string, in which case it has to match one substring of the lfilenames
    print("WARNING: Could not convert second argument into an int, for variable where_to_look_in_filename")
    print("\nWill try to match to something in the filename\nThis can break real bad")
    substring_to_look_for1 = sys.argv[2]
    look_for_matching_string1 = True
    
    # Also just small check that there is a substring:
    if '_' in substring_to_look_for1:
        print("Error: there is something wrong (an underscore _ ) in this substring:")
        print(substring_to_look_for1)
    # sys.exit(1)

try:
    where_to_look_in_filename2 = int(sys.argv[3])
except ValueError:
    #maybe it is a string, in which case it has to match one substring of the lfilenames
    print("WARNING: Could not convert second argument into an int, for variable where_to_look_in_filename")
    print("\nWill try to match to something in the filename\nThis can break real bad")
    substring_to_look_for2 = sys.argv[3]
    look_for_matching_string2 = True
    
    # Also just small check that there is a substring:
    if '_' in substring_to_look_for2:
        print("Error: there is something wrong (an underscore _ ) in this substring:")
        print(substring_to_look_for2)
    # sys.exit(1)
    
lfilename = sys.argv[4:]

for i,filename in enumerate(lfilename):
    print( "Opening: ", filename)
    if look_for_matching_string1 == True:
        filename_split = filename.split('_')
        filename_split_nonumbers = [ x.rstrip(string.digits+'.') for x in filename_split]
        print("filename_split_nonumbers:",filename_split_nonumbers)
        if substring_to_look_for1 in filename_split_nonumbers:
            where_to_look_in_filename1 = filename_split_nonumbers.index(substring_to_look_for1)
        else:
            print("Error: Could not find substring", substring_to_look_for1,"matching substring")
            sys.exit(1)
    if look_for_matching_string2 == True:
        filename_split = filename.split('_')
        filename_split_nonumbers= [ x.rstrip(string.digits+'.') for x in filename_split]
        print("filename_split_nonumbers:",filename_split_nonumbers)
        if substring_to_look_for2 in filename_split_nonumbers:
            where_to_look_in_filename2 = filename_split_nonumbers.index(substring_to_look_for2)
        else:
            print("Error: Could not find substring", substring_to_look_for2,"matching substring")
            sys.exit(1)
    # POS BASED ON BITSTRING LENGTH
    # pos_violinplot = list(filter(lambda x:'bitstrABlen' in x, filename.split('_')))
    # pos_violinplot = int( pos_violinplot[0][len('bitstrABlen'):] )
    
    if where_to_look_in_filename1>=0: 
        pos_violinplot1 = filename.split('_')[where_to_look_in_filename1] 
        print("Got pos_violinplot1 =",pos_violinplot1)
        if ".txt" in pos_violinplot1:
            pos_violinplot1 = pos_violinplot1.replace('.txt', '')
        pos_violinplot1 = float( ''.join(c for c in pos_violinplot1 if c in string.digits+'.' ) )
        # pos_violinplot = int( pos_violinplot[1:] )
    else: 
        pos_violinplot1 = i
    
    if where_to_look_in_filename2>=0: 
        pos_violinplot2 = filename.split('_')[where_to_look_in_filename2] 
        print("Got pos_violinplot2 =",pos_violinplot2)
        if ".txt" in pos_violinplot2:
            pos_violinplot2 = pos_violinplot2.replace('.txt', '')
        pos_violinplot2 = float( ''.join(c for c in pos_violinplot2 if c in string.digits+'.' ) )
        # pos_violinplot = int( pos_violinplot[1:] )
    else: 
        pos_violinplot2 = i

    if pos_violinplot1 in dpos_violinplot:
        # pos_violinplot+=0.5
        # lpos_violinplot[index()]
        # lpos_violinplot[lpos_violinplot.index(pos_violinplot)] -=0.5
        dpos_violinplot[pos_violinplot1][pos_violinplot2] = dict()
    else:
        dpos_violinplot[pos_violinplot1]=dict()
        dpos_violinplot[pos_violinplot1][pos_violinplot2] = dict()
    print ("Got pos1 = ",pos_violinplot1, ", pos2 =",pos_violinplot2)
    #sys.exit(1)
    
    # lpos_violinplot.append(pos_violinplot)
    with open(filename,"r") as fin:
        lgenome_length=[]
        lgenome_noFlen=[]
        lF=[]
        lA=[]
        lB=[]
        lABfield=[]
        dpop={}
        l_antib=[]
        l_percentiles=[]
        # l_nr_abs=[]
        for line in fin:
            line =line.split()
            if "n" not in line and len(line)>genome_pos_in_file:
                try:
                    # lgenome.append(line[genome_pos_in_file])
                    genome = line[genome_pos_in_file]
                    genlen = len(genome)
                    lgenome_length.append(genlen)
                    if maxlen<genlen: maxlen = genlen

                    # if 'F' not in line[genome_pos_in_file]:
                    #     lgenome_noFlen.append(genlen)
                    lF.append(genome.count('F'))
                    lA.append(genome.count('A'))
                    lB.append(genome.count('B'))
                    
                except:
                    print("Error in the line:", line)
                    pass
            #print line
            #sys.exit(1)  
        # Populate dictionary
        dpos_violinplot[pos_violinplot1][pos_violinplot2]['lgenome_length'] = lgenome_length
        dpos_violinplot[pos_violinplot1][pos_violinplot2]['lF'] = lF
        dpos_violinplot[pos_violinplot1][pos_violinplot2]['lA'] = lA
        dpos_violinplot[pos_violinplot1][pos_violinplot2]['lB'] = lB
        
        
#For now just a bad trick for pos, later we make a mini x axis to put them properly, for now equally spaced - sorry.
lpos = [ pos1 + pos1*y for pos1 in sorted(dpos_violinplot) for y in range(len(dpos_violinplot[pos1]))]
print (lpos)
#lpos, llgenome_length_sorted,llF_sorted,llA_sorted,llB_sorted = (list(t) for t in zip(*sorted(zip(lpos_violinplot, llgenome_length,llF,llA,llB))))
# sys.exit(1)    


from mpl_toolkits.mplot3d.art3d import Poly3DCollection # New import
# make 3d axes
fig = plt.figure()
ax = fig.gca(projection='3d')
def CalculateMutRate(pc_mutrate,nb):
    totp=0.
    for x in range(int(nb)):
        totp+=pc_mutrate*(1.-float(pc_mutrate))**x
    return totp

ll_data =[[]]
tag='lB'
for i,pos1 in enumerate(sorted(dpos_violinplot)):
    for pos2 in sorted(dpos_violinplot[pos1]):
        # print(dpos_violinplot[pos1][pos2][tag])
        nb = np.median(dpos_violinplot[pos1][pos2][tag])
        # pc_mutrate = float(pos1)
        # ll_data[i].append( 0. if nb == 0 else CalculateMutRate(pc_mutrate,nb) / float(nb) )
        ll_data[i].append(CalculateMutRate(float(pos1),nb) )
    if len(sorted(dpos_violinplot)) > i+1:
        ll_data.append([])
print("ll_data: ",ll_data)

# for x,ldata in enumerate(ll_data):
#     X = np.array([x for _ in range(len(ldata))])
#     Y = np.array([y for y in range(len(ldata)) ])
#     Z = np.array(ldata)
#     print(X.shape, Y.shape, Z.shape)
#     ax.plot(Y,X,Z,lw=1,color=u'#2ca02c')
#     # verts = [(Y[i],X[i],Z[i]) for i in range(len(Y))] + [(Y.max(),x,0),(Y.min(),x,0)]
#     # ax.add_collection3d(Poly3DCollection([verts],color='#2ca02c',alpha=0.3)) # Add a polygon instead of fill_between

ll_data = np.array(ll_data)

ax.set_xticks([x for x in range(len(dpos_violinplot[pos1])) ])
ax.set_yticks([y for y in range(len(dpos_violinplot)) ])
ax.set_xticklabels([x for x in sorted(dpos_violinplot[pos1])]) # column makes x-axis
ax.set_yticklabels([x for x in sorted(dpos_violinplot)])
ax.set_xlabel("break_infl") # column makes x-axis - column is break_infl
ax.set_ylabel("break_prob") # <- breakprob is rows, rows make y axis

ll_dataT = np.transpose(ll_data)
for x,ldata in enumerate(ll_dataT):
    X = np.array([x for _ in range(len(ldata))])
    Y = np.array([y for y in range(len(ldata)) ])
    Z = np.array(ldata)
    print(X.shape, Y.shape, Z.shape)
    ax.plot(X,Y,Z,lw=4,color='black',zorder = 50)
    verts = [(X[i],Y[i],Z[i]) for i in range(len(X))] + [(X.max(),Y.max(),0),(X.min(),Y.min(),0)]
    ax.add_collection3d(Poly3DCollection([verts],color='black',alpha=0.3,zorder = 15)) # Add a polygon instead of fill_between
# ax.set_zlim(0, 30)
# plt.title("Normalised!")

# Plot vertical lines
for x,ldata in enumerate(ll_dataT):
    for y,z in enumerate(ldata):
        ax.plot([x,x],[y,y],[0,z],lw=1,color='black',alpha=0.5)


# my_cmap = sns.color_palette('RdYlBu' , as_cmap=True)
# norm = mpl.colors.Normalize(vmin=ll_data.min(), vmax=ll_data.max() )
# colors = my_cmap(norm(ll_data))



ddata={}

#data is written as
# filename (from which you get break prob and brinfl)
# X axis
# ab
# repl

print("Opening: ",sys.argv[1])

with open(sys.argv[1]) as fin:
    lines=fin.readlines()
    for i in range(0, len(lines), 4):
        filename, lX,lab,lrepl = lines[i:i + 4]
        breakprob_pos = filename.find('breakprob')
        breakprob= filename[breakprob_pos+len('breakprob'):].split('_')[0]
        
        breakinfl_pos = filename.find('breakinfl')
        breakinfl= filename[breakinfl_pos+len('breakinfl'):].split('_')[0][:-4] # to remove .txt
        
        # if breakinfl == '0.0001':
        #     print("Got it")
        
        lX = [int(x) for x in lX[6:-2].split(', ')]
        lab = [float(x) for x in lab[16:-2].split(', ')]
        lrepl = [float(x) for x in lrepl[16:-2].split(', ')]
        if breakprob not in ddata:
            ddata[breakprob]={}
        if breakinfl not in ddata[breakprob]:
            ddata[breakprob][breakinfl]=[lX,lab,lrepl]
        else:
            print("Error, already got this")

# for i,breakprob in enumerate(sorted(ddata)):
#     for j,breakinfl in enumerate(sorted(ddata[breakprob])):
#         print(breakprob,breakinfl)

# plot condensed data
import copy
def find_median(ll):
    sum_lab = sum(ll)
    run_sum = 0.
    for i,ab in enumerate(ll):
        run_sum += ab/sum_lab
        if run_sum>0.5:
            return i

dcondata = copy.deepcopy(ddata)
for breakprob in ddata:
    for breakinfl in ddata[breakprob]:
        lab = ddata[breakprob][breakinfl][1]
        lrepl = ddata[breakprob][breakinfl][2]
        dcondata[breakprob][breakinfl] = find_median(lrepl) - find_median(lab)
lcondata= np.zeros((len(dcondata),len(dcondata[breakprob])), dtype=float)
for i,breakprob in enumerate(sorted(dcondata)):
    for j,breakinfl in enumerate(sorted(dcondata[breakprob])):
        lcondata[i,j] = dcondata[breakprob][breakinfl]

# fig = plt.figure()
# ax = plt.gca()


vmax= 10
divnorm=colors.DivergingNorm(vmin=0, vcenter=5, vmax=vmax)
my_cmap = sns.color_palette('RdYlBu' , as_cmap=True)
colors = my_cmap(divnorm(lcondata))


#to plot a flat 3d surface as a kind of pcolormesh we need to centre the data points
# so if we have, say, 4x5 data points, we'll need 5x6 corners for the squares, 
# at the center of which there are our data points
X, Y = np.meshgrid([x-0.5 for x in range(1+len(dpos_violinplot[pos1]))], [y-0.5 for y in range(1+len(dpos_violinplot))])
surf = ax.plot_surface(X, Y, np.zeros_like(X), facecolors=colors,zorder = 1)
# heatmap_data= ax.pcolormesh( ll_mutrate, cmap = sns.cubehelix_palette(start=2, rot=0, dark=0, light=.95, reverse=True, as_cmap=True) )





# from matplotlib import colors

# vmax= 10
# divnorm=colors.DivergingNorm(vmin=0, vcenter=5, vmax=vmax)
# # my_cmap = sns.color_palette("vlag", as_cmap=True)
# # my_cmap = sns.diverging_palette(250, 30, l=65, as_cmap=True)
# my_cmap = sns.color_palette('RdYlBu' , as_cmap=True)

# bla = ax.pcolormesh(lcondata, vmax=vmax,cmap=my_cmap,norm=divnorm)

# ax.set_xticks([0.5+i for i in range(4)])
# ax.set_xticklabels([x for x in sorted(dcondata[breakprob])]) # rows are y axis
# ax.set_yticks([0.5+i for i in range(5)])
# ax.set_yticklabels([x for x in sorted(dcondata)])
# cbar = plt.colorbar(bla,extend="max")

# cbar.set_label('$\Delta$ growth genes (AB producers, replicating bacteria)')




plt.show()

sys.exit(1)
