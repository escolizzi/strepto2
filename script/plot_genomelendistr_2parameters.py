#!/usr/bin/python3

'''
Makes plots of genome and genes distribution.
Puts on the x axis a parameter of your choice. 
Puts on a mini axis within the x axis a second parameter of your choice. 
The parameters you want to use must be in the file name, and you have to specify its position or the subtring

Usage: 
./plot_genomelendistr_2parameters.py [pos1] [pos2] [list of files]

-also added antib production
'''

import matplotlib as mpl
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
import numpy as np

def Birthrate(a,b):
    hd = bin(a^b).count('1')
    return math.exp(-hd**2.)

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


try:
    where_to_look_in_filename1 = int(sys.argv[1])
except ValueError:
    #maybe it is a string, in which case it has to match one substring of the lfilenames
    print("WARNING: Could not convert second argument into an int, for variable where_to_look_in_filename")
    print("\nWill try to match to something in the filename\nThis can break real bad")
    substring_to_look_for1 = sys.argv[1]
    look_for_matching_string1 = True
    
    # Also just small check that there is a substring:
    if '_' in substring_to_look_for1:
        print("Error: there is something wrong (an underscore _ ) in this substring:")
        print(substring_to_look_for1)
    # sys.exit(1)

try:
    where_to_look_in_filename2 = int(sys.argv[2])
except ValueError:
    #maybe it is a string, in which case it has to match one substring of the lfilenames
    print("WARNING: Could not convert second argument into an int, for variable where_to_look_in_filename")
    print("\nWill try to match to something in the filename\nThis can break real bad")
    substring_to_look_for2 = sys.argv[1]
    look_for_matching_string2 = True
    
    # Also just small check that there is a substring:
    if '_' in substring_to_look_for2:
        print("Error: there is something wrong (an underscore _ ) in this substring:")
        print(substring_to_look_for2)
    # sys.exit(1)
    
lfilename = sys.argv[3:]

if calculate_susceptibility==True: 
    print("No support for calculating susceptibility now, please disable this option")
    sys.exit(1)

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
            if calculate_susceptibility==True: 
                    lABfield.extend(line[2].split(",")[:-1]) #[:-1] because empty after last comma ends up in the list
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
                    
                    if calculate_susceptibility==True: 
                        # Antibiotics in the field
                        if len(line[nr_ab_pos_in_file])>0:
                            antib_at_pos = line[ab_list_pos_in_file]
                            l_antib.extend( antib_at_pos.split(',')[:-1] ) #get ab deposited in the field, remove the last bc it's empty after last comma
                        # Antibiotics in the genome
                        ab_genome = list( set( line[gene_ab_pos_in_file].split(',')[:-1] ) )
                        if len(ab_genome):
                            if  len(ab_genome)==1 and '' in ab_genome: continue #empty list with just a comma fails, this catches it
                            colony_val = int(line[3])
                            # this is the colony antibiotic pan-genome
                            # we will exclude the very rare abs later - if they - happen, because they are not representative of the colony
                            if colony_val in dpop:
                                dpop[colony_val].extend(ab_genome)
                                #print("Extended: dpop[colony_val] = ",dpop[colony_val])
                            else:
                                dpop[colony_val] = ab_genome
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
        
        if calculate_susceptibility==True: 
            llABfield.append( len(set(lABfield))/float(len(lABfield)) ) # AB diversity (unique/tot)
            llABfield.append(Counter(lABfield))
        
        l_susceptibility=[]
        
        if calculate_susceptibility==True:
            l_antib=list(set(l_antib))
            for dcounter, colony in enumerate(dpop):
                # counted_abs = Counter(dpop[colony])
                set_abs = set(dpop[colony])
                # print("counted abs: ",counted_abs)
                #dpop[colony] = [ x for x in set_abs if counted_abs[x]>0]
                dpop[colony]=list(set(dpop[colony]))  # we should exclude rare genes, because they are not representative - but how many are there?
                                                # FOR NOW EXCLUDE NOBODY, THEN WE'LL SEE
                #for each colony, check it min birthrate to each ab in the field
                # so go through all antib in the field
                l_susceptibility.append(0)
                if len(dpop[colony]) != 0:
                    for field_ab in l_antib:
                        try:
                            l_icel_birthrates=[Birthrate( int(i_ab),int(field_ab) ) for i_ab in dpop[colony]] # This has to be kept in sync with simulations
                        except:
                            print(int(field_ab))
                            print(dpop[colony])
                            print([int(i_ab) for i_ab in dpop[colony]])
                            sys.exit(1)
                        # now l_icel_birthrates contains a bunch of growth rates, for each AB in the field.
                        # given the AB genes in the bacterium, it grows with the rate given by the max resistance for this particular AB
                        # so if this max growth rate is less than 0.25 (or the number you want), then we say it is susceptible
                        if max(l_icel_birthrates) < 0.5 : l_susceptibility[-1] += 1
        # print("l_suscept = ", l_susceptibility)
        if calculate_susceptibility==True:
            try:
                quantiles = np.quantile(l_susceptibility, [0.05,0.25,0.5,0.75,0.95] )
            except Exception as error: 
                print("Got some error, it probably means that there are no antibiotics, you might want to check this if you think there should be ABs")
                print("Or that you set calculate_susceptibility to False")
                print("this will just go on with: quantiles = [0.,0.,0.,0.,0.] ")
                print("The error message is:")
                print(error)
                quantiles = [0.,0.,0.,0.,0.]
                
            ll_howmany_ab.append( len(l_antib) )
            ll_susceptibility.append( quantiles )
        
if calculate_susceptibility==True: 
    lpos, llgenome_length_sorted,llF_sorted,llA_sorted,llB_sorted,llABfield_sorted,ll_howmany_ab_sorted,ll_susceptibility_sorted = (list(t) for t in zip(*sorted(zip(lpos_violinplot, llgenome_length,llF,llA,llB,llABfield,ll_howmany_ab,ll_susceptibility))))
else:
    #For now just a bad trick for pos, later we make a mini x axis to put them properly, for now equally spaced - sorry.
    lpos = [ pos1 + pos1*y for pos1 in sorted(dpos_violinplot) for y in range(len(dpos_violinplot[pos1]))]
    print (lpos)
    #lpos, llgenome_length_sorted,llF_sorted,llA_sorted,llB_sorted = (list(t) for t in zip(*sorted(zip(lpos_violinplot, llgenome_length,llF,llA,llB))))
    # sys.exit(1)    

# fig, ax = plt.subplots(1,3)
# plt.subplots_adjust(wspace=0.2)
fig = plt.figure()
spec = gridspec.GridSpec(nrows=1, ncols=6*3, figure=fig)
ax = [ fig.add_subplot(spec[:, 0:3]), 
       fig.add_subplot(spec[:, 4]), 
       fig.add_subplot(spec[:, 6:9]), 
       fig.add_subplot(spec[:, 10]), 
       fig.add_subplot(spec[:, 12:15]),
       fig.add_subplot(spec[:, 16])]


from mpl_toolkits.axes_grid1 import make_axes_locatable
def Make2dplot(tag, ax_nr):
# tag='lF'
    ll_data =[[]]

    for i,pos1 in enumerate(sorted(dpos_violinplot)):
        for pos2 in sorted(dpos_violinplot[pos1]):
            # print(dpos_violinplot[pos1][pos2][tag])
            ll_data[i].append( np.median(dpos_violinplot[pos1][pos2][tag]) )
        if len(sorted(dpos_violinplot)) > i+1:
            ll_data.append([])
    # ll_data=[ dpos_violinplot[pos1][pos2][tag] for pos1 in sorted(dpos_violinplot) for pos2 in sorted(dpos_violinplot[pos1])]
    print("Got tag ", tag)
    print("ddpos keys:", [x for x in sorted(dpos_violinplot)])
    
    for i,pos1 in enumerate(sorted(dpos_violinplot)):
        print("ddpos[",pos1,"] keys:", [x for x in sorted(dpos_violinplot[pos1])])
    print("lldata")
    print(ll_data)
    # sys.exit(1)
    mycmap=sns.color_palette("hot", as_cmap=True)
    ltags = ["lF","lB","lA"]
    lvminvmax=[(0,10),(0,30),(40,300)]
    vmin,vmax= lvminvmax[ltags.index(tag)]
    
    heatmap_data= ax[ax_nr].pcolormesh( np.array(ll_data), cmap=mycmap, vmin=vmin,vmax=vmax )
    ax[ax_nr].set_aspect('equal')
    # bla= ax[ax_nr].imshow( np.array(ll_data) )
    
    # bla = [[0,1,2],[3,4,5],[6,7,8]]
    
    # bla= ax[ax_nr].pcolor( np.transpose(bla) )
    # bla= ax[ax_nr].pcolor( bla )
    # fig.colorbar(bla, cax=ax[ax_nr+1])
    ax[ax_nr].title.set_text(tag)
    # ax[ax_nr].plot([0,1,2],[0,1,2])
    # plt.show()
    # ax[ax_nr].set_xticklabels([0,1,2])
    # ax[ax_nr].set_yticklabels([3,4,5])
    # plt.show()
    ax[ax_nr].set_xticks([0.5+x for x in range(len(dpos_violinplot[pos1])) ])
    ax[ax_nr].set_yticks([0.5+x for x in range(len(dpos_violinplot)) ])
    # print(ax[ax_nr].get_xticks())
    # sys.exit(1)
    ax[ax_nr].set_xticklabels([x for x in sorted(dpos_violinplot[pos1])]) # column makes x-axis
    ax[ax_nr].set_yticklabels([x for x in sorted(dpos_violinplot)])
    ax[ax_nr].set_xlabel("break_infl") # column makes x-axis - column is break_infl
    ax[ax_nr].set_ylabel("break_prob") # <- breakprob is rows, rows make y axis
    
    # ax[ax_nr].set_xscale('log')
    # ax[ax_nr].set_yscale('log')
    # cax = fig.add_axes([ax[ax_nr].get_position().x1+0.01,ax[ax_nr].get_position().y0,0.02,ax[ax_nr].get_position().height])
    fig.colorbar(heatmap_data, cax=ax[ax_nr+1])
    ax[ax_nr].title.set_text(tag)
    # cax.set_axis_off()
    # fig.colorbar(bla,cax=ax[ax_nr+1],fraction=0.046, pad=0.04)
    
Make2dplot('lF', 0)
Make2dplot('lB', 2)
Make2dplot('lA', 4)

plt.show()
sys.exit(1)

split_axes=False     
if split_axes==True:
    print( "Length", len(dpos_violinplot),len(llgenome_length))
    print("Warning splitting y-axis btwn 80 and 120")
    print("Warning not printing susceptibility")
    
    # fig, (ax,ax2) = plt.subplots(2,1,sharex=True)
    fig2 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=5, figure=fig2)
    ax = fig2.add_subplot(spec2[0, 0])
    ax2 = fig2.add_subplot(spec2[1:, 0],sharex=ax)
    
    
    ax.plot(lpos,[np.median(x) for x in llF_sorted],marker='o',color=u'#1f77b4',label="nr. growth-promoting genes")
    ax.fill_between(lpos, [np.quantile(x, .25) for x in llF_sorted], [np.quantile(x, .75) for x in llF_sorted], alpha=0.5,color=u'#1f77b4')
    ax.plot(lpos,[np.median(x) for x in llA_sorted],marker='o',color=u'#ff7f0e',label="nr. antibiotic genes")
    ax.fill_between(lpos, [np.quantile(x, .25) for x in llA_sorted], [np.quantile(x, .75) for x in llA_sorted], alpha=0.5,color=u'#ff7f0e')
    ax.plot(lpos,[np.median(x) for x in llB_sorted],marker='o',color=u'#2ca02c', label="nr. fragile sites")
    ax.fill_between(lpos, [np.quantile(x, .25) for x in llB_sorted], [np.quantile(x, .75) for x in llB_sorted], alpha=0.5,color=u'#2ca02c')
    # plot the same data on both axes
    ax2.plot(lpos,[np.median(x) for x in llF_sorted],marker='o',color=u'#1f77b4',label="nr. growth-promoting genes")
    ax2.fill_between(lpos, [np.quantile(x, .25) for x in llF_sorted], [np.quantile(x, .75) for x in llF_sorted], alpha=0.5,color=u'#1f77b4')
    ax2.plot(lpos,[np.median(x) for x in llA_sorted],marker='o',color=u'#ff7f0e',label="nr. antibiotic genes")
    ax2.fill_between(lpos, [np.quantile(x, .25) for x in llA_sorted], [np.quantile(x, .75) for x in llA_sorted], alpha=0.5,color=u'#ff7f0e')
    ax2.plot(lpos,[np.median(x) for x in llB_sorted],marker='o',color=u'#2ca02c', label="nr. fragile sites")
    ax2.fill_between(lpos, [np.quantile(x, .25) for x in llB_sorted], [np.quantile(x, .75) for x in llB_sorted], alpha=0.5,color=u'#2ca02c')
    # zoom-in / limit the view to different portions of the data
    
    ax.set_ylim(60, 120)  # outliers only
    ax2.set_ylim(0, 60)  # most of the data

    # hide the spines between ax and ax2
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.xaxis.tick_top()
    ax.tick_params(labeltop='off')  # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    plt.show()
    sys.exit(1)
    
    
    ax.set_ylim([0,None])

    ax.set_xticks(lpos_violinplot)

    if where_to_look_in_filename>=0:
        # ax.set_xscale('log')
        ax.set_xticklabels([str(x) for x in lpos_violinplot])
        print("x Positions: ", lpos_violinplot)
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

    
    
else:
    
    # print( "Length", len(lpos_violinplot),len(llgenome_length))
    fig, ax = plt.subplots(1)
    # widht = 100.

    # lpos, llgenome_length_sorted,llF_sorted,llA_sorted,llB_sorted,llABfield_sorted,ll_howmany_ab_sorted,ll_susceptibility_sorted = (list(t) for t in zip(*sorted(zip(lpos_violinplot, llgenome_length,llF,llA,llB,llABfield,ll_howmany_ab,ll_susceptibility))))


    # ll_data=llgenome_length_sorted
    # ax.plot(lpos,[np.median(x) for x in ll_data],marker='o',color=u'#808080')
    # ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=u'#808080')
    
    def Plotquartiles(tag, mycolor,mylabel):
        ll_data=[ dpos_violinplot[pos1][pos2][tag] for pos1 in sorted(dpos_violinplot) for pos2 in sorted(dpos_violinplot[pos1])]
        ax.plot(lpos,[np.median(x) for x in ll_data],marker='o',color=mycolor,label=mylabel)
        ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=mycolor)
    Plotquartiles('lF',u'#1f77b4',"nr. growth-promoting genes")
    Plotquartiles('lA',u'#ff7f0e',"nr. antibiotic genes")
    Plotquartiles('lB',u'#2ca02c',"nr. fragile sites")
    
    # ll_data=[ dpos[pos1][pos2]['lF'] for pos1 in sorted(dpos_violinplot) for pos2 in sorted(dpos_violinplot[pos1])]
    # ax.plot(lpos,[np.median(x) for x in ll_data],marker='o',color=u'#1f77b4',label="nr. growth-promoting genes")
    # # print(lpos,[np.median(x) for x in ll_data])
    # ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=u'#1f77b4')

    # ll_data=llA_sorted
    # ax.plot(lpos,[np.median(x) for x in ll_data],marker='o',color=u'#ff7f0e',label="nr. antibiotic genes")
    # ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=u'#ff7f0e')
    
    # # ll_data=llABfield_sorted
    # #for later, notice that llABfield_sorted now is a list of dictionaries, in each dict, key = each antib, value = how many times it appeared

    # # ax.plot(lpos,[x for x in ll_data],marker='o',color=u'#aa7f0e',label="diversity antibiotic in field")
    # # ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=u'#aa7f0e')

    # ll_data=llB_sorted
    # ax.plot(lpos,[np.median(x) for x in ll_data],marker='o',color=u'#2ca02c', label="nr. fragile sites")
    # ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=u'#2ca02c')
    
    if calculate_susceptibility==True:     
        
        # ax.plot(lpos,ll_howmany_ab,color='dimgrey',label="Tot. nr. different antibiotics produced")
        ll_susceptibility_sorted = np.array(ll_susceptibility_sorted).T
        ax1 = ax.twinx()
        ax1.plot(lpos,[x/float(y) for x,y in zip(ll_susceptibility_sorted[2],ll_howmany_ab_sorted)],c='darkred')
        ax1.fill_between(lpos, [x/float(y) for x,y in zip(ll_susceptibility_sorted[1],ll_howmany_ab_sorted)],[x/float(y) for x,y in zip(ll_susceptibility_sorted[3],ll_howmany_ab_sorted)],color='darkred',alpha=0.5)
        ax1.set_ylabel("Susceptibility", color='darkred')
        ax1.tick_params(axis='y', labelcolor='darkred')
        ax1.set_ylim(0,None)

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
    
    logxaxis=True
    if logxaxis==True:
        ax.set_xscale('log')
        # ax.set_xticks(lpos_violinplot)
    # else:
    
    ax.set_xticks(lpos)
    
    
    # if where_to_look_in_filename>=0:
    #     # ax.set_xscale('log')
    #     ax.set_xticklabels([str(x) for x in lpos_violinplot])
    #     print("x Positions: ", lpos_violinplot)
    # else:
    #     #remove stuff in front that's always the same
    #     stripped_lfilename=[  '_'.join(x for x in fn.split('_')[3:]) for fn in lfilename]
    #     #add \n once in a while to make label readable
    #     stripped_lfilename2=[]
    #     for fn in stripped_lfilename:
    #         # print("old fn: ",fn)
    #         newfn=""
    #         count=1
    #         for i,a in enumerate(fn):
    #             if a == '_':
    #                 count +=1 
    #                 if count%2 == 1:
    #                     newfn+='\n'
    #                 else:
    #                     newfn+=a
    #             else: 
    #                 newfn+=a
    #         print("newfn: ",newfn)
    #         stripped_lfilename2.append(newfn)
        
    #     print("stripped_lfilename2: ")
    #     for x in stripped_lfilename2:
    #         print(x)
    #     print("xticks labelled with: ",lfilename)
    #     ax.set_xticklabels( stripped_lfilename2  ,rotation=0,ha='center')
    #     print("ticks at pos: ",lpos_violinplot)

    # ax.set_xlabel("Colony development time")
    ax.legend()

    plt.tight_layout()
    plt.show()
    sys.exit(1)



    #######################################      THE END     #####################################


