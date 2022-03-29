#!/usr/bin/python3

'''
Adapted from the competition plotting of the MC experiments.
Plot whether the evolved guy (red) or the generalist (blue) wins
This version will need all the data to make a complete plot
'''

import sys,math,os,subprocess,random
#from PIL import Image
import matplotlib as mpl
mpl.use('QT5Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import seaborn as sns
from matplotlib.lines import Line2D
from collections import OrderedDict
import matplotlib.gridspec as gridspec
import numpy as np

import seaborn as sns
#manager = plt.get_current_fig_manager()
#print manager
#sys.exit(1)
#########################
###                   ###
### ---   BEGIN   --- ###
###                   ###
#########################



colours=["firebrick","royalblue", "darkgoldenrod", "green", "salmon", "lightskyblue","orchid"]
filename=""

#fig, (ax0,ax1) = plt.subplots(nrows=2, sharex=True)
# fig, ax = plt.subplots(nrows=1,ncols=2)
# plt.subplots_adjust(wspace=0, hspace=0)

fig = plt.figure()
spec = gridspec.GridSpec(nrows=2, ncols=40, figure=fig, wspace=0 )
ax=[fig.add_subplot(spec[0, 0:7]), 
    fig.add_subplot(spec[0, 7:40]), 
    fig.add_subplot(spec[1, :]) ]

if len(sys.argv) <4:
    print ("This is the program 'plot_competition_popsize_dist.py'")
    print ("Usage: ./plot_division_time.py <figure name> <season duration> <nr replicates> <datafile 1> <datafile 2> <...>")
    sys.exit(1)
else:
    figname=sys.argv[1]
    season=int(sys.argv[2])
    nrreplicates=int(sys.argv[3])
redcols=["firebrick", "orangered", "lightcoral"]
blucols=["mediumblue", "cornflowerblue", "skyblue"]
filenr=0

heatmapdata=np.zeros((6,10,),dtype=float)    #np.zeros((5*7,),dtype=float)
heatmapdata[:] = np.nan


kab=["0.00001", "0.00002", "0.00003", "0.00004", "0.00005", "0.00006", "0.00007", "0.0001", "0.0002", "0.0004"]
krep=["0.01", "0.025", "0.05", "0.0625", "0.075", "0.1"]

for files in sys.argv[4:]:
    
    popsize=[ 0 for x in range(2)]
    print(files)
    print (filenr) 
    filetype=0
    

    fileparts=files.split('_')
    for el in fileparts:
        if "kab" in el:
            kabval=el[3:]
            if ".txt" in el:
                kabval=kabval[:-4]
            kabind=kab.index(kabval)
        elif "krep" in el:
            krepval=el[4:]
            krepind=krep.index(krepval)
    print("Got kab = ",kabval,", krep = ",krepval)
    print("Indexes: ",kabind,krepind)
    if krepind==2 and kabind==3:
        print("    *******      ")
        print(" ")
    print("This value is currently:",heatmapdata[krepind][kabind])
    
    with open(files,"r") as fin:

        #read first line first to get the first time point (could probably more clever, but hey)
        seasoncount=0
        prevtime=0
        registered=0
        index=int(filenr/nrreplicates)
        print (index)

        for line in reversed(fin.readlines()):
            
            line1=line.split(' ')
            if(len(line1)>3 or ',' in line):
                time = int(line1[0])
                nrab=int(line1[1])
                sporetype=int(line1[3])
                #count population size
                if sporetype:
                    popsize[sporetype-1]+=1

            else:
                filetype=1
                if (int(line1[1])>10*int(line1[2])):
                    if np.isnan(heatmapdata[krepind][kabind]):
                        heatmapdata[krepind][kabind] = 1./float(nrreplicates)
                    else:
                        heatmapdata[krepind][kabind] += 1./float(nrreplicates)
                    print("Evolved won")
                elif (10*int(line1[1])<int(line1[2])):
                    if np.isnan(heatmapdata[krepind][kabind]):
                        heatmapdata[krepind][kabind] = -1./float(nrreplicates)
                    else:
                        heatmapdata[krepind][kabind] -= 1./float(nrreplicates)
                    print("Competitor won")
                else:
                    if np.isnan(heatmapdata[krepind][kabind]):
                        heatmapdata[krepind][kabind] = 0
                    
                break

          #  if time!=prevtime:
        if(not filetype):
            if (popsize[0]>10*popsize[1]):
                if np.isnan(heatmapdata[krepind][kabind]):
                    heatmapdata[krepind][kabind] = 1./float(nrreplicates)
                    
                else:
                    heatmapdata[krepind][kabind] += 1./float(nrreplicates)
                print("Evolved won")
            elif (10*popsize[0]<popsize[1]):
                if np.isnan(heatmapdata[krepind][kabind]):
                    heatmapdata[krepind][kabind] = -1./float(nrreplicates)
                    
                else:
                        heatmapdata[krepind][kabind] -= 1./float(nrreplicates)
                print("Competitor won")
            else:
                if np.isnan(heatmapdata[krepind][kabind]):
                    heatmapdata[krepind][kabind] = 0

    filenr+=1

# defined above
# kab=["0.00001", "0.00002", "0.00003", "0.00004", "0.00005", "0.00006", "0.00007", "0.0001", "0.0002", "0.0004"]
# krep=["0.01", "0.025", "0.05", "0.0625", "0.075", "0.1"]

# first part of heat map uses values 
krep_0 = ["0.01", "0.05", "0.0625", "0.075", "0.1"]
X_0 = [0.]+  [(float(x)+float(y))/2. for x,y in zip(krep_0[:-1],krep_0[1:])] + [0.1 + (0.1-0.075)/2.] #the last elements exceed a bit, to make it look nicer
kab_0= [ "0.00001" , "0.00002" , "0.00003" , "0.00004" , "0.00005" , "0.00006" , "0.00007" ]
Y_0 = [0.00001/2.] + [(float(x)+float(y))/2. for x,y in zip(kab_0[:-1],kab_0[1:])] + [0.00007+0.00001/2.] #the last elements exceed a bit, to make it look nicer

print("HEllo, this is the whole heatmap")
print(heatmapdata)
print("Hello this hetamp[2,3]: ", heatmapdata[2,3])

heatmapdata_0 = heatmapdata[[0,2,3,4,5],:7]
print("Going to plot the following heatmapdata_0: ")
print(heatmapdata_0)
# sys.exit(1)

heatmap_colormap=sns.diverging_palette(240, 10, n=9,as_cmap=True) #plt.cm.RdYlBu_r
X_mesh0 = ax[0].pcolormesh(Y_0, X_0, heatmapdata_0,edgecolor='white',cmap=heatmap_colormap)
ax[0].set_yticks([float(x) for x in krep])
ax[0].set_xticks([float(x) for x in kab_0])
ax_0_xticklabels = [x if i%2==0 else "" for i,x in enumerate(kab_0)] # skips some tick labels for clarity
ax[0].set_xticklabels(ax_0_xticklabels)

heatmapdata_1 = heatmapdata[:,7:]
print("Going to plot the following heatmapdata_1: ")
print(heatmapdata_1)
krep_1 = krep[:]
X_1 = [0.]+  [(float(x)+float(y))/2. for x,y in zip(krep_1[:-1],krep_1[1:])] + [0.1 + (0.1-0.075)/2.] #the last elements exceed a bit, to make it look nicer
kab_1= [ "0.0001" , "0.0002" , "0.0004" ]
Y_1 = [0.00007+0.00001/2.] + [(float(x)+float(y))/2. for x,y in zip(kab_1[:-1],kab_1[1:])] + [0.0004 + 0.0001/2.] #the last elements exceed a bit, to make it look nicer
X_mesh1 = ax[1].pcolormesh(Y_1, X_1, heatmapdata_1,edgecolor='white',cmap=heatmap_colormap)
ax[1].yaxis.set_ticklabels([])
ax[1].set_yticks([])
ax[1].set_yticks([], minor=True) # for minor ticks
ax[1].spines['left'].set_visible(False) 
ax[1].set_xticks([float(x) for x in kab_1])
# ax[1].set_xticklabels([float(x) for x in kab_1])


cmap = heatmap_colormap#plt.cm.bwr

#genome1
# print()
# print()
# print("USING GENOME 1\n\n")
# krepl = 0.058
# mean_ab_prod_wt = 0.00011729
# std_ab_prod_wt = 0.000064293

#genome2
# print()
# print()
# print("USING GENOME 2\n\n")
# krepl = 0.068
# mean_ab_prod_wt = 0.00010468
# std_ab_prod_wt = 0.000039623

#genome3
print()
print()
print("USING GENOME 3\n\n")
krepl = 0.06
mean_ab_prod_wt = 0.000067862
std_ab_prod_wt = 0.000030915

#adds average and std dev in a separate plot with equal dimensions, to overlay with inkscape or so.
ax[2].plot( [ mean_ab_prod_wt - std_ab_prod_wt , mean_ab_prod_wt + std_ab_prod_wt ],[krepl,krepl],lw=2,c='darkgoldenrod')
ax[2].scatter( [ mean_ab_prod_wt ],[krepl], s = 10)
ax[2].set_ylim([0 , 0.1 + (0.1-0.075)/2.])
ax[2].set_xlim([0.00001/2. , 0.0004 + 0.0001/2.])
ax[2].set_xticks([float(x) for x in kab])
ax[2].set_yticks([float(x) for x in krep])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
# cax = fig.add_axes([0.8, 0.15, 0.05, 0.3])

n_colors = 9
bounds = [x for x in range(n_colors+1)] 

# print("bounds: ", bounds)
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

my_colorbar = fig.colorbar(
                mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
                cax=cbar_ax, orientation='vertical',
                label="who wins")
             
# tick_locs = (np.arange(len(bounds)) + 0.5)*(len(bounds)-1)/len(bounds)
tick_locs = [ 0.5+x for x in range(n_colors) ]
tick_labels=[ x-4 for x in range(n_colors) ] 
my_colorbar.set_ticks(tick_locs)
my_colorbar.set_ticklabels(tick_labels)
plt.show()
sys.exit(1)
# ax[0].pcolormesh([0].append()

#display consensus image
#imcons=np.reshape(heatmapdata,(5,7))
#fig=plt.figure() #!
#fig.set_size_inches(1, 1, forward = False) #!
#ax = plt.Axes(fig, [0., 0., 1., 1.]) #!
##ax.set_axis_off() #!
#fig.add_axes(ax) #!
x=[0,1,2,3,4,5,6]
xlabels=["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "1", "2", "4"]
#ylabels=["0.01", "0.05", "0.0625", "0.075", "0.1"]
ax.set_xlabel("Kab (10^-4)")
ax.set_ylabel("Krep")
plt.xticks(x, xlabels)
plt.yticks(x, krep)
im=plt.imshow(heatmapdata, cmap=plt.cm.bwr, origin='lower')


plt.clim(-1, 1)
plt.colorbar(im, ticks=[-1, -0.5, 0, 0.5, 1]) # adding the colorbar on the right
plt.savefig(figname, bbox_inches='tight')
#plt.close()
#plt.show()
