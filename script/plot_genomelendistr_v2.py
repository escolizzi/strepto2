#!/usr/bin/python3.6

'''
Makes plots of genome and genes distribution.
Positions must be specified in some way.. for now hard coded
-also added antib production
'''


import sys,math,os,string
from subprocess import Popen, PIPE
from collections import Counter

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

genome_logo2=[]
totpos2=0
max_totpos2=0

genome_pos_in_file = 4

lpos_violinplot=[]

try:
    where_to_look_in_filename = int(sys.argv[1])
except:
    print("Error: Could not convert second argument into an int, for variable where_to_look_in_filename")
    sys.exit(1)

lfilename = sys.argv[2:]
for i,filename in enumerate(lfilename):
    print( "Opening: ", filename)
    
    #  POS BASED ON BITSTRING LENGTH
    # pos_violinplot = list(filter(lambda x:'bitstrABlen' in x, filename.split('_')))
    # pos_violinplot = int( pos_violinplot[0][len('bitstrABlen'):] )
    
    if where_to_look_in_filename>=0: 
        pos_violinplot = filename.split('_')[where_to_look_in_filename] 
        print("Got pos_violinplot =",pos_violinplot)
        pos_violinplot = int( ''.join(c for c in pos_violinplot if c.isdigit() ) )
        # pos_violinplot = int( pos_violinplot[1:] )
    else: 
        pos_violinplot = i


    if pos_violinplot in lpos_violinplot:
        # pos_violinplot+=0.5
        # lpos_violinplot[index()]
        # lpos_violinplot[lpos_violinplot.index(pos_violinplot)] -=0.5
        lpos_violinplot.append( pos_violinplot+1 )
    else:
        lpos_violinplot.append(pos_violinplot)
    print ("Got pos = ",pos_violinplot)
    
    # lpos_violinplot.append(pos_violinplot)
    with open(filename,"r") as fin:
        lgenome_length=[]
        lgenome_noFlen=[]
        lF=[]
        lA=[]
        lB=[]
        lABfield=[]
        for line in fin.readlines():
            line =line.split()
            lABfield.extend(line[2].split(",")[:-1]) #[:-1] because empty after last comma ends up in the list
            if "n" not in line:
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
        llABfield.append( len(set(lABfield))/float(len(lABfield)) ) # AB diversity (unique/tot)
        llABfield.append(Counter(lABfield))

print( "Length", len(lpos_violinplot),len(llgenome_length))
fig, ax = plt.subplots(1)
# widht = 100.

lpos, llgenome_length_sorted,llF_sorted,llA_sorted,llB_sorted,llABfield_sorted = (list(t) for t in zip(*sorted(zip(lpos_violinplot, llgenome_length,llF,llA,llB,llABfield))))


ll_data=llgenome_length_sorted
# ax.plot(lpos,[np.median(x) for x in ll_data],marker='o',color=u'#808080')
# ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=u'#808080')

ll_data=llF_sorted
ax.plot(lpos,[np.median(x) for x in ll_data],marker='o',color=u'#1f77b4',label="nr. growth-promoting genes")
ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=u'#1f77b4')

ll_data=llA_sorted
ax.plot(lpos,[np.median(x) for x in ll_data],marker='o',color=u'#ff7f0e',label="nr. antibiotic genes")
ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=u'#ff7f0e')

ll_data=llABfield_sorted
#for later, notice that llABfield_sorted now is a list of dictionaries, in each dict, key = each antib, value = how many times it appeared

# ax.plot(lpos,[x for x in ll_data],marker='o',color=u'#aa7f0e',label="diversity antibiotic in field")
# ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=u'#aa7f0e')

ll_data=llB_sorted
ax.plot(lpos,[np.median(x) for x in ll_data],marker='o',color=u'#2ca02c')
ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=u'#2ca02c')


# w = 1.
# width = lambda p, w: 10**(np.log10(p)+w/2.)-10**(np.log10(p)-w/2.)
# vln1 = ax.violinplot(llgenome_length,positions = lpos_violinplot,   widths=width(lpos_violinplot,w), showextrema=False, showmedians=True)

# vln1 = ax.violinplot(llgenome_length,positions = [x-0.375 for x in lpos_violinplot],   widths=width([x-0.375 for x in lpos_violinplot],w), showextrema=False, showmedians=True)
# vlnF = ax.violinplot(llF,            positions = [x-0.125 for x in lpos_violinplot], bw_method=1, widths=width([x-0.125 for x in lpos_violinplot],w), showextrema=False, showmedians=True)
# vlnA = ax.violinplot(llA,            positions = [x+0.125 for x in lpos_violinplot], bw_method=1, widths=width([x+0.125 for x in lpos_violinplot],w), showextrema=False, showmedians=True)
# vlnB = ax.violinplot(llB,            positions = [x+0.375 for x in lpos_violinplot], bw_method=1, widths=width([x+0.375 for x in lpos_violinplot],w), showextrema=False, showmedians=True)

# for pc in vln1['bodies']:
#     pc.set_facecolor(u'#808080')
#     pc.set_alpha(0.5)
# vln1['cmedians'].set_edgecolor(u'#808080')

# for pc in vlnF['bodies']:
#     pc.set_facecolor(u'#1f77b4')
#     pc.set_alpha(0.5)
# vlnF['cmedians'].set_edgecolor(u'#1f77b4')

# for pc in vlnA['bodies']:
#     pc.set_facecolor(u'#ff7f0e')
#     pc.set_alpha(0.5)
# vlnA['cmedians'].set_edgecolor(u'#ff7f0e')

# for pc in vlnB['bodies']:
#     pc.set_facecolor(u'#2ca02c')
#     pc.set_alpha(0.5)
# vlnB['cmedians'].set_edgecolor(u'#2ca02c')

ax.set_ylim([0,None])

# ax.set_yscale('symlog')
ax.set_xticks(lpos_violinplot)

if where_to_look_in_filename>=0:
    ax.set_xscale('log')
    ax.set_xticklabels([str(x) for x in lpos_violinplot])
else:
    #remove stuff in front that's always the same
    stripped_lfilename=[  '_'.join(x for x in fn.split('_')[3:]) for fn in lfilename]
    #add \n once in a while to make label readable
    stripped_lfilename2=[]
    for fn in stripped_lfilename:
        # print("old fn: ",fn)
        newfn=""
        count=1
        for i,a in enumerate(fn):
            if a == '_':
                count +=1 
                if count%2 == 1:
                    newfn+='\n'
                else:
                    newfn+=a
            else: 
                newfn+=a
        print("newfn: ",newfn)
        stripped_lfilename2.append(newfn)
    
    print("stripped_lfilename2: ")
    for x in stripped_lfilename2:
      print(x)
    print("xticks labelled with: ",lfilename)
    ax.set_xticklabels( stripped_lfilename2  ,rotation=0,ha='center')
    print("ticks at pos: ",lpos_violinplot)

# ax.set_xlabel("Colony development time")
ax.legend()

plt.tight_layout()
plt.show()
sys.exit(1)



#######################################      THE END     #####################################


