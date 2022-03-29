#!/usr/bin/python3

'''
Makes plots of genome and genes distribution.
Puts on the x axis a parameter of your choice. 
The parameter you want to use must be in the file name, and you have to specify its position or the subtring

Usage: 
./plot_genomelendistr_2parameters.py [pos1] [pos2] [list of files]

-also added antib production
'''

import matplotlib as mpl
mpl.use('Qt5Agg')

import sys,math,os,string,random
from subprocess import Popen, PIPE
from collections import Counter

#from PIL import Image
# import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
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

dpos_violinplot=dict()

look_for_matching_string = False

calculate_susceptibility=False
#  Wishful thinking without better input management
# if "-nosusceptbility" in sys.argv: 
#     calculate_susceptibility==False
if calculate_susceptibility==False: 
    print("WARNING, NOT CALCULATING SUSCEPTIBILITY")

try:
    where_to_look_in_filename = int(sys.argv[1])
except ValueError:
    #maybe it is a string, in which case it has to match one substring of the lfilenames
    print("WARNING: Could not convert second argument into an int, for variable where_to_look_in_filename")
    print("\nWill try to match to something in the filename\nThis can break real bad")
    substring_to_look_for = sys.argv[1]
    look_for_matching_string = True
    
    # Also just small check that there is a substring:
    if '_' in substring_to_look_for:
        print("Error: there is something wrong (an underscore _ ) in this substring:")
        print(substring_to_look_for)
    # sys.exit(1)

lfilename = sys.argv[2:]

for i,filename in enumerate(lfilename):
    print( "Opening: ", filename)
    if look_for_matching_string == True:
        filename_split = filename.split('_')
        filename_split_nonumbers = [ x.rstrip(string.digits+'.') for x in filename_split]
        print("filename_split_nonumbers:",filename_split_nonumbers)
        if substring_to_look_for in filename_split_nonumbers:
            where_to_look_in_filename = filename_split_nonumbers.index(substring_to_look_for)
        else:
            print("Error: Could not find substring", substring_to_look_for,"matching substring")
            sys.exit(1)
    
    # POS BASED ON BITSTRING LENGTH
    # pos_violinplot = list(filter(lambda x:'bitstrABlen' in x, filename.split('_')))
    # pos_violinplot = int( pos_violinplot[0][len('bitstrABlen'):] )
    
    if where_to_look_in_filename>=0: 
        pos_violinplot = filename.split('_')[where_to_look_in_filename] 
        print("Got pos_violinplot =",pos_violinplot)
        if ".txt" in pos_violinplot:
            pos_violinplot = pos_violinplot.replace('.txt', '') #removes trailing .txt if present
        pos_violinplot = ''.join(c for c in pos_violinplot if c in string.digits+'.' )
        # pos_violinplot = int( pos_violinplot[1:] )
    else: 
        pos_violinplot = i

    # ATTEMPTS TO DO SOMETHING WITH DUPLICATE X positions
    # if pos_violinplot in lpos_violinplot:
    #     # pos_violinplot+=0.5
    #     # lpos_violinplot[index()]
    #     # lpos_violinplot[lpos_violinplot.index(pos_violinplot)] -=0.5
    #     lpos_violinplot.append( pos_violinplot+1 )
    # else:
    #     lpos_violinplot.append(pos_violinplot)
    n=0
    while pos_violinplot in dpos_violinplot:
        n+=1
        pos_violinplot = str(pos_violinplot)+ '_'+str(n)
    
    print ("Got pos =",pos_violinplot)
    dpos_violinplot[pos_violinplot]=dict()
    
    # lpos_violinplot.append(pos_violinplot)
    with open(filename,"r") as fin:
        lgenome_length=[]
        lgenome_noFlen=[]
        lF=[]
        lA=[]
        lB=[]
        lH=[]
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
                    lH.append(genome.count('H'))
                    
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
                    pass
            #print line
            #sys.exit(1)  
        dpos_violinplot[pos_violinplot]['lgenome_length'] = lgenome_length
        dpos_violinplot[pos_violinplot]['lF'] = lF
        dpos_violinplot[pos_violinplot]['lA'] = lA
        dpos_violinplot[pos_violinplot]['lB'] = lB
        dpos_violinplot[pos_violinplot]['lH'] = lH
        
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
    print("Error, current update does not support calculating susceptibility, come again later and fix me")
    sys.exit(1)
    lpos, llgenome_length_sorted,llF_sorted,llA_sorted,llB_sorted,llABfield_sorted,ll_howmany_ab_sorted,ll_susceptibility_sorted = (list(t) for t in zip(*sorted(zip(lpos_violinplot, llgenome_length,llF,llA,llB,llABfield,ll_howmany_ab,ll_susceptibility))))
else:
    lpos = [ pos for pos in sorted(dpos_violinplot)]
    
    # lpos, llgenome_length_sorted,llF_sorted,llA_sorted,llB_sorted = (list(t) for t in zip(*sorted(zip(lpos_violinplot, llgenome_length,llF,llA,llB))))

print (lpos)
# sys.exit(1)

# HANDLING REPETITIONS and placements
# create dpos_without_repetition
dpos_without_repetition = {}
scaling_factor = 10.
for pos in dpos_violinplot:
    thispos = float(pos.split('_')[0])
    if thispos in dpos_without_repetition:
        dpos_without_repetition[thispos] += 1
    else:
        dpos_without_repetition[thispos] = 1
xaxis_scale = scaling_factor * ( max([float(x) for x in dpos_without_repetition]) - min([float(x) for x in dpos_without_repetition]) ) / len(dpos_without_repetition)
print ("dpos_without_repetition",dpos_without_repetition)

#create list of float(pos) that is sorted, but keep a list that says how it is sorted
lpos_tags=[x for x in dpos_violinplot] #these are the tags
lpos_float=[float(pos.split('_')[0]) for pos in lpos_tags ] #these are the tags converted to numbers for sorting
#now we sort lpos_tags using lpos_float's value (there may be repetitions such as [1,1,...] but maybe it does not matter for now)
lpos_tags_sorted = [x for _, x in sorted(zip(lpos_float, lpos_tags))] #actually this suffices to call the array in the right order
# lpos_tags_sorted_positions = [lpos_tags.index(x) for x in lpos_tags] #this says where the original elements of lpos_tags are.
print ("lpos_tags_sorted",lpos_tags_sorted)

lpos=[]
for pos in sorted(dpos_without_repetition):
    print (type( pos))
    nrep=dpos_without_repetition[pos]
    if nrep>1:
        val_pos = float(pos)
        to_add = [(x-(nrep-1)/2)/xaxis_scale for x in range(nrep)]
        lpos.extend( [val_pos+x for x in to_add] )
    else: 
        lpos.append(float(pos))

print("lpos",lpos)
print("IF THESE lpos don't look right, please check the parameter scaling_factor, which is now = ",scaling_factor)
print("Chances are you can make the categories neatly separated by playing with it")
print()

# sys.exit(1)
# At this point lpos is still sorted alphanumerically, so 1, 10, 2 ...
#we need to sort this but also rember the order, 

# sorted_key_dpos = [x for x in sorted(dpos_violinplot)]  #these are the handles, they are char
# lpos_sorted = [x for x in sorted(lpos)] #these are the positions on the plot, in the same order as sorted_key_dpos
# lpos_sorted_pos = [lpos.index(x) for x in lpos_sorted]  #this is the order in which they should be called, so and lpos_sorted[ lpos_sorted_pos] and dpos_violinplot[sorted_key_dpos[ lpos_sorted_pos ]] give x and y for plots

# print("sorted_key_dpos", sorted_key_dpos)
# print("lpos", lpos)
# print ("lpos_sorted_pos",lpos_sorted_pos)

# sys.exit(1 )
split_axes=False     
if split_axes==True:
    print( "Length", len(lpos_violinplot),len(llgenome_length))
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
    
    fig, ax = plt.subplots(1)
    
    def PlotData(tag, color, label):
        # this_lpos = [lpos[ x] for x in lpos_sorted_pos]
        # print()
        # print("this_lpos", [lpos[ x] for x in lpos_sorted_pos])
        # print("[sorted_key_dpos[ x ] for x in lpos_sorted_pos]: ",[sorted_key_dpos[ x ] for x in lpos_sorted_pos])
        ll_data=[ dpos_violinplot[x][tag] for x in lpos_tags_sorted]
        ax.plot(lpos,[np.median(x) for x in ll_data],marker='o',color=color,label=label)
        # ax.fill_between(lpos, [np.quantile(x, .25) for x in ll_data], [np.quantile(x, .75) for x in ll_data], alpha=0.5,color=color)
        for pos,x in zip(lpos,ll_data):
           ax.plot([pos,pos],[np.quantile(x, .25),np.quantile(x, .75)], color=color, alpha=0.7)
        
        # sys.exit(1)
    PlotData('lF' , u'#1f77b4' , "nr. growth-promoting genes")
    PlotData('lA' , u'#ff7f0e' , "nr. antibiotic genes")
    PlotData('lB' , u'#2ca02c' , "nr. fragile sites")
    if sum(lB)!=0:
        PlotData('lH' , 'darkgoldenrod' , "nr. homeostatic genes")
    
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
    
    logxaxis=False
    if logxaxis==True:
        ax.set_xscale('log')
        
        # ax.set_xticks(lpos_violinplot)
    # else:
    
    ax.set_xticks([float(x) for x in sorted(dpos_without_repetition)])
    ax.set_xticklabels([str(x) for x in sorted(dpos_without_repetition)])
    
        
    # if where_to_look_in_filename>=0:
    #     # ax.set_xscale('log')
    #     ax.set_xticklabels([str(x) for x in sorted(dpos_violinplot)])
    #     print("x Positions: ", dpos_violinplot)
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
    #     print("ticks at pos: ",[x for x in sorted(dpos_violinplot)])

    # ax.set_xlabel("Colony development time")
    ax.legend()

    plt.tight_layout()
    plt.show()
    sys.exit(1)



    #######################################      THE END     #####################################


