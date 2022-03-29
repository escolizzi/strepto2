#!/usr/bin/python3

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
llgenome=[]
lgenome_noFlen=[]
llgenome_noFlen=[]

llF=[]
llA=[]
llB=[]


lav_F=[]
lav_A=[]
lav_B=[]
lstd_F=[]
lstd_A=[]
lstd_B=[]

# filename = sys.argv[1]
genome_logo=[]
l_totpos=[]
l_totpos2=[]
max_totpos=0
maxlen=0
lgenome_length=[]
llgenome_length=[]

genome_logo2=[]
totpos2=0
max_totpos2=0

genome_pos_in_file = 4

lpos_violinplot=[]

lfilename = sys.argv[1:]
for filename in lfilename:
    print( "Opening: ", filename)
    # pos_violinplot = filename.split('_')
    pos_violinplot = list(filter(lambda x:'bitstrABlen' in x, filename.split('_')))
    # pos_violinplot = int( pos_violinplot[0][len('bitstrABlen24'):-len('.txt')] )
    pos_violinplot = int( pos_violinplot[0][len('bitstrABlen'):] )
    if pos_violinplot in lpos_violinplot:
        # pos_violinplot+=0.5
        # lpos_violinplot[index()]
        # lpos_violinplot[lpos_violinplot.index(pos_violinplot)] -=0.5
        lpos_violinplot.append( pos_violinplot+1 )
    else:
        lpos_violinplot.append(pos_violinplot)
    print ("Got bistrlen = ",pos_violinplot)
    
    # lpos_violinplot.append(pos_violinplot)
    with open(filename,"r") as fin:
        lgenome_length=[]
        lgenome_noFlen=[]
        lF=[]
        lA=[]
        lB=[]
        for line in fin.readlines():
            line =line.split()
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
                pass
            #print line
            #sys.exit(1)  
        llgenome_length.append(lgenome_length)
        # llgenome_noFlen.append(lgenome_noFlen)
        llF.append(lF)
        llA.append(lA)
        llB.append(lB)
        
print( "Length", len(lpos_violinplot),len(llgenome_length))
fig, ax = plt.subplots(1)
vln1 = ax.violinplot(llgenome_length,positions = [x-0.375 for x in lpos_violinplot],  widths=1., showextrema=False, showmedians=True)
vlnF = ax.violinplot(llF,            positions = [x-0.125 for x in lpos_violinplot], bw_method=1, widths=1., showextrema=False, showmedians=True)
vlnA = ax.violinplot(llA,            positions = [x+0.125 for x in lpos_violinplot], bw_method=1, widths=1., showextrema=False, showmedians=True)
vlnB = ax.violinplot(llB,            positions = [x+0.375 for x in lpos_violinplot], bw_method=1, widths=1., showextrema=False, showmedians=True)

for pc in vln1['bodies']:
    pc.set_facecolor(u'#808080')
    pc.set_alpha(0.5)
vln1['cmedians'].set_edgecolor(u'#808080')

for pc in vlnF['bodies']:
    pc.set_facecolor(u'#1f77b4')
    pc.set_alpha(0.5)
vlnF['cmedians'].set_edgecolor(u'#1f77b4')

for pc in vlnA['bodies']:
    pc.set_facecolor(u'#ff7f0e')
    pc.set_alpha(0.5)
vlnA['cmedians'].set_edgecolor(u'#ff7f0e')

for pc in vlnB['bodies']:
    pc.set_facecolor(u'#2ca02c')
    pc.set_alpha(0.5)
vlnB['cmedians'].set_edgecolor(u'#2ca02c')

ax.set_xticks(lpos_violinplot)
ax.set_xticklabels([str(x) for x in lpos_violinplot])
ax.set_ylim([0,None])
# ax.set_yscale('symlog')
plt.show()
sys.exit(1)



#######################################      THE END     #####################################







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
