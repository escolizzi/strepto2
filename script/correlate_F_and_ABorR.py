#!/usr/bin/python3

'''
Makes plots of genome and genes distribution.
Positions must be specified in some way.. for now hard coded
-also added antib production
'''
import matplotlib as mpl
mpl.use('Qt5Agg')

import sys,math,os,string,random,copy
from subprocess import Popen, PIPE
from collections import Counter

#from PIL import Image
# mpl.use('GTK3Agg')
from matplotlib import pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np


def Gene2Number(gene):
    if gene=='F':return 0
    if gene=='A':return 1
    if gene=='B':return 2

dg={}
count=0
broken_repl=0
filename = sys.argv[1]
# tot_people_alive_pergnm={}
line_number=0
maxtime = 1250 # 1000#2500

with open(filename,"r") as fin:
    
    for line in fin:
        line_number+=1
        line = line.split()
        time = int(line[0])
        if time >= maxtime : break
        
        # E = extinction
        if "E:" in line[1]:
            #we be done early!
            for key in dg:
                if dg[key]["death"]==maxtime:
                    dg[key]["death"]=time+1 # if someone is already died we don't artificially extend its life time
            break
        tag=line[2]
        if tag not in dg:
            #set line 2
            # I = initial 
            if "I:" in line[1]:
                try:
                    dg[tag]={"gnm":line[3], "birth":time, "death":maxtime, "A":0, "R":0, "offs_tag":[], "parent":"n"} #set death to 2500, this will be true unless death is recorded in log
                except:
                    print("Error, here is line: ", line)
                    print("line_number = ", line_number)
                    # sys.exit(1)
                    pass
            else:
                if len(line)>3:
                    print("Error, this is not an initial guy.")
                    print(line)
                    sys.exit(1)
                else:
                    continue # it's a guy with zero length genome
        else:
            #update life
            # R = replication
            if "R:" in line[1]:
                #it replicated
                dg[tag]["R"]+=1
                if len(line)>5:
                    offspr_tag = line[4]
                    dg[tag]["offs_tag"].append(offspr_tag)
                    if offspr_tag in dg:
                        print("Error, this offspring tag is already in dg")
                    else:
                        if len(line)>5:
                            dg[offspr_tag]={"gnm":line[5], "birth":time, "death":maxtime, "A":0, "R":0, "offs_tag":[], "parent":tag} #set death to maxtime, this will be true unless death is recorded in log
            # A = antibiotic production
            elif "A:" in line[1]:
                #it made an antibiotic
                dg[tag]["A"]+=1
            # D = death
            elif "D:" in line[1]:
                if tag in dg:
                    dg[tag]["death"]=time
        

print("Done, reading data. In total so many people lived:", len(dg))

print("Successful genomes - in order:")
l_nr_offs=[]
d_genomes={} #dictionary with genomes a keys; values: total abundance, earliest found
d_gnm_and_offs={} #dict with parent offspring relation. key: genome; values: descending genomes
for key in dg:
    gnm = dg[key]["gnm"]
    if gnm in d_genomes:
        d_genomes[gnm]["howmany_offs"]+=len(dg[key]["offs_tag"])
        d_genomes[gnm]["earliest"] = min( dg[key]["birth"] , d_genomes[gnm]["earliest"] )
        d_genomes[gnm]["howmany"]+=1
        d_genomes[gnm]["ab_prod"]+=dg[key]["A"]
    else: d_genomes[gnm]={"ab_prod": dg[key]["A"],  "howmany": 1 , "howmany_offs":len(dg[key]["offs_tag"]) , "earliest": dg[key]["birth"] }

    if gnm in d_gnm_and_offs:
        d_gnm_and_offs[gnm].extend( [dg[ot]["gnm"] for ot in dg[key]["offs_tag"]] )
    else:
        d_gnm_and_offs[gnm] = [ dg[ot]["gnm"] for ot in dg[key]["offs_tag"] ]

#just mutational relations
for gnm in d_gnm_and_offs:
    d_gnm_and_offs[gnm]=set(d_gnm_and_offs)


#Get succesful ancestors
l_ancestor = []
nr_offspring_to_count_as_successful = 500 # when is an ancestor successful? when it makes so many offsprings
for gnm in d_genomes:
    # if d_genomes[gnm]["howmany"]>50: print(gnm,": ",d_genomes[gnm])
    if d_genomes[gnm]["earliest"]==0 and d_genomes[gnm]["howmany_offs"]>nr_offspring_to_count_as_successful:
        l_ancestor.append(gnm)
# Now go get these guys from sys.argv[2] = data_bla.txt
if len(sys.argv)<=2:
    print("No ancestor genomes are searched")
else:
    print("Going to look for founder bacteria")
    with open(sys.argv[2],"r") as fin2:
        for line in fin2:
            line=line.split()
            if int(line[0]) != 0: break
            if line[4] in l_ancestor:
                print( line[4],line[5])

def set_size(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)


# GENOME PLOTTING OR LOLLYPOP PLOTTING REPRODUCTION vs AB.
# ALSO ANCESTOR OFFSPRING RELATION deltaF and deltab
l_repr_ab_for_lollyplot=[]
ldeltaf_po=[]
ldeltaab=[]
tot_deltaab=0
print_genomes = False
if print_genomes:
    import dnaplotlib as dpl
    dr = dpl.DNARenderer()
    part_renderers = dr.SBOL_part_renderers()
sp = {'type':'EmptySpace', 'name':'S1', 'fwd':True, 'opts':{'x_extent':3}} #spacer
Fgn = {'type':'UserDefined', 'name':'ud' , 'fwd':True, 'opts':{'color':(0.0,0.5,0.9),'x_extent':2, 'y_extent':7 } } #blue gene - fitness
Agn = {'type':'UserDefined', 'name':'ud' , 'fwd':True, 'opts':{'color':(0.9,0.5,0.1),'x_extent':2, 'y_extent':7 } } #orange gene - antibiotic
Bp = {'type':'Ribozyme', 'name':'bla', 'fwd':True, 'opts':{'color':(0.05,0.6,0.4),'x_extent':7, 'y_extent':11 } }  #break point
l_drawgenes=[Fgn,Agn,Bp]
l_genes=['F','A','B']
designs=[]
print("Ancestor list: ", l_ancestor)
import copy
# print some genomes and ther offspring
def MakeDesign(genomeseq):
    # print ("Seq:", genomeseq)
    design = [sp]
    for pos,x in enumerate(genomeseq):
        # if pos == 0:
        #     # len_block=1 
        #     # continue
        #     design.append( l_drawgenes[l_genes.index(x)] )
        #     # design.append(sp)
        if pos>0 and x != genomeseq[pos-1] : 
            design.append(sp)
        design.append( l_drawgenes[l_genes.index(x)] )
        # else:
        #     which_gene = l_genes.index(genomeseq[pos-len_block])
        #     draw_gene = copy.deepcopy( l_drawgenes[ which_gene ] ) # copies dict
        #     print ("len_block = ", len_block)
        #     if which_gene == 'B':
        #         for _ in len_block:
        #             design.append( draw_gene )

        #     else:
        #         draw_gene['opts']['x_extent'] = len_block*draw_gene['opts']['x_extent']
        #         #   print(draw_gene['opts']['x_extent'])
        #         design.append( draw_gene )
        #     design.append(sp)
        #     len_block=1
    design.append(sp)
    return design

print("Keys d_genome", [key for key in d_genomes])

for ancestor in l_ancestor:
    if d_genomes[ancestor]["howmany"]>2500:
        print("Ancestor:") 
        print(ancestor, d_genomes[ancestor]["howmany"],d_genomes[ancestor]["ab_prod"])
    l_repr_ab_for_lollyplot.append( [d_genomes[ancestor]["howmany"],d_genomes[ancestor]["ab_prod"], ancestor ] )
    
    # print("Offspring",) 
    designs.append(MakeDesign(ancestor))
    tot_deltaab += d_genomes[ancestor]["ab_prod"]
    for offs_seq in sorted(d_gnm_and_offs[ancestor],reverse=True):
        if offs_seq == ancestor: continue
        # print(offs_seq,d_genomes[offs_seq]["howmany"],d_genomes[offs_seq]["ab_prod"])
        l_repr_ab_for_lollyplot.append( [d_genomes[offs_seq]["howmany"],d_genomes[offs_seq]["ab_prod"],offs_seq] )
        ldeltaf_po.append( ancestor.count("F") - offs_seq.count("F") )
        ldeltaab.append( d_genomes[ancestor]["ab_prod"] - d_genomes[offs_seq]["ab_prod"])
        tot_deltaab += d_genomes[offs_seq]["ab_prod"]
        designs.append(MakeDesign(offs_seq))
        
        # design = [sp]
        # for x in offs_seq:
        #     design.append( l_drawgenes[l_genes.index(x)] )
        #     design.append(sp)
        # designs.append(design)
    ldeltaab=[x/float(tot_deltaab) for x in ldeltaab]
    # print("Genome: ", designs[-1])
if print_genomes:
    fig = plt.figure(figsize=(8,6))
    gs = gridspec.GridSpec( len(designs), 2+ sum([ x['opts']['x_extent'] for x in designs[0] ]) )

    ax_dna=[]
    for i in range( len(designs) ):
        if i==0:
            ax_dna.append( plt.subplot( gs[i,1:] ) )
        else:
            ax_dna.append( plt.subplot( gs[i, 1:  sum([ x['opts']['x_extent'] for x in designs[i] ])   ] ) )
        # ax_dna2 = plt.subplot(gs[1])

        start, end = dr.renderDNA(ax_dna[i], designs[i], part_renderers)
        print (i, start,end)
        ax_dna[i].set_xlim([start, end])
        ax_dna[i].set_ylim([-15,28])
        ax_dna[i].set_aspect('equal')
        ax_dna[i].set_xticks([])
        ax_dna[i].set_yticks([])
        ax_dna[i].axis('off')

    for i,my_ax in enumerate(ax_dna): 
        set_size(len(designs[i]),5, ax=my_ax)
    # plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.01)

    # fig.tight_layout()
    plt.savefig("pretty.pdf")
    sys.exit(1)
# else:
    #lollyplots

print()

# PARENT-OFFSPRING DIFFERENCES, AS A FUNCTION OF Delta(F)
print ("Hello 1")
max_nr_F=0
# sys.exit(1)
# ldeltaf_po=[]
# ldeltaab=[]
# tot_deltaab=0
# print("Warning - making only for ancestor - rest of colony relation")
# if False and len(l_ancestor)==1:
#     ancestor=l_ancestor[0]
for key in dg:
    max_nr_F = max(dg[key]["gnm"].count("F") , max_nr_F)
#         l_this_deltaab=[]
#         tot_deltaab = dg[key]["A"]
#         howmany_offs1=0
#         for offs_tag in dg[key]["offs_tag"]:
#             if offs_tag in dg:
#                 # print(dg[key]["A"] , dg[offs_gnm]["A"])
                
#                 # deltaab = dg[key]["A"]/float(dg[key]["death"] - dg[key]["birth"]) - dg[offs_tag]["A"]/float(dg[key]["death"] - dg[key]["birth"])
#                 # Just a try: divide by total AB prod by the colony- BELOW
#                 deltaab = dg[key]["A"] - dg[offs_tag]["A"]
#                 tot_deltaab += dg[offs_tag]["A"]
#                 howmany_offs1+=1
#                 # if deltaab> 50000: 
#                 #     print("Mother gnm, ab: ",key,dg[key]["A"], " Offsp gnm, ab", offs_tag, dg[offs_tag]["A"] )
#                 l_this_deltaab.append( deltaab )
#                 # ldeltagr.append( dg[key]["R"] - dg[offs_gnm]["R"] )
#                 ldeltaf_po.append( dg[key]["gnm"].count("F") - dg[offs_tag]["gnm"].count("F") )
#                 if ldeltaf_po[-1]==1 and deltaab>0:
#                     print(dg[key]["gnm"])
#                     print(dg[offs_tag]["gnm"])
#                     print()
#             #else:
#             #    # else this guy did nothing - what a whimp!
#         #normalise Delta(AB prod) by tot ab production
#         # we have to normalise the last howmany_offs1
#         # print ("TOT AB: ", tot_deltaab)
#         # print("")
#         # l_this_deltaab=[ x/float(tot_deltaab ) if tot_deltaab>0 else 0. for x in l_this_deltaab ]
#         ldeltaab.extend( l_this_deltaab ) # Now ldeltaf_po will be in sync with ldeltaab

print ("Hello 2")
# Make averages and std
max_deltaf=max(ldeltaf_po)
print("max_deltaf parent-offspr: ",max_deltaf)
l1 = [[] for _ in range(1+max(ldeltaf_po))] # thanks past me - WTF is l1?
# print("l1: ",l1)
for df,ab in zip(ldeltaf_po,ldeltaab):
    l1[df].append(ab)
#Sometimes some positions in l1 are not filled

print ("Hello 2.1")
# print ("l1 = ", l1)
lavr_deltaab = [np.mean(x) if x!=[] else 0 for x in l1 ]
lstd_deltaab = [np.std(x) if x!=[] else 0 for x in l1 ]
print ("Hello 2.2")
lmedian_deltaab = [np.median(x) if x!=[] else 0 for x in l1 ]
print ("Hello 2.3")
l25perc_deltaab = [np.percentile(x,25) if x!=[] else 0 for x in l1 ]
print ("Hello 2.4")
l75perc_deltaab = [np.percentile(x,75)  if x!=[] else 0 for x in l1 ]

print ("Hello 3")
#first thing to check: Who makes AB and who reproduces:
lab=[0 for _ in range(max_nr_F+1)]
lrepl=[0 for _ in range(max_nr_F+1)]
# print(lab)
ldeltaF=[]
ldeltaF2=[]
# d_genomes={}
for key in dg:
    me_nr_F = dg[key]["gnm"].count("F")
    lab[ me_nr_F ]+= dg[key]["A"] 
    lrepl[ me_nr_F ]+= dg[key]["R"]
    
    # me_genome = dg[key]["gnm"]
    # if me_genome in d_genomes:
    #     d_genomes[me_genome]["howmany"]+=1
    # else:
    #     d_genomes[me_genome]={"howmany":1, "offs_gnm":[]}
    
    if len(dg[key]["offs_tag"])>0:
        this_deltaf=[]
        for offs_tag in dg[key]["offs_tag"]:
            # print( dg[key]["gnm"] )
            # print( dg[offs_tag]["gnm"] )
            deltaFfract = ( dg[key]["gnm"].count("F") - dg[offs_tag]["gnm"].count("F") )/(float( dg[key]["gnm"].count("F") ))
            ldeltaF.append( deltaFfract )
            # this_deltaf.append( deltaFfract )
            
            # d_genomes[me_genome]["offs_gnm"].append( dg[offs_tag]["gnm"] )
            # if dg[offs_tag]["gnm"] in d_genomes:
            #     d_genomes[ dg[offs_tag]["gnm"] ]["howmany"]+=1
            # else:
            #     d_genomes[ dg[offs_tag]["gnm"] ]={"howmany":1, "offs_gnm":[]}
        # len_thisdeltaf = float(len(this_deltaf))
        # hist_thisdeltaf,bins = np.histogram( this_deltaf, bins = np.linspace(0.,1.1,num=11) )
        # ldeltaF2.append([x/len_thisdeltaf for x in hist_thisdeltaf])

hist_df, bin_edgesdf = np.histogram( ldeltaF, bins = np.linspace(0.,1.05,num=21) )
hist_df = [x/float(len(dg)) for x in hist_df]

# ldeltaF2 = [list(i) for i in zip(*ldeltaF2)] # so we have it per-bin - bins defined by the linespace 3 lines above
# print("Here0")
# lmedian_deltaf2 = [np.median(x) for x in ldeltaF2]
# lmean_deltaf2 = [np.mean(x) for x in ldeltaF2]
# l25per_deltaf2 = [np.percentile(x,25) for x in ldeltaF2]
# l75per_deltaf2 = [np.percentile(x,75) for x in ldeltaF2]

print("Plotting")

# fig, ax = plt.subplots(2,2)
fig = plt.figure()
gs = gridspec.GridSpec(4, 4)
ax=[]
ax.append([fig.add_subplot(gs[0:2, 0:2]), fig.add_subplot(gs[0:2,2:4])] )
ax.append([fig.add_subplot(gs[2:4, 0:2]),fig.add_subplot(gs[2:4,2]),fig.add_subplot(gs[2:4,3]) ])

print(ax)

sum_lab=sum(lab)
sum_repl=sum(lrepl)

# print("WARNING TURN THESE TWO LINES BACK ON WHEN YOU ARE DONE!!!")
#normalised
ax[0][0].plot(range(max_nr_F+1),[x/float(sum_lab) for x in lab],label="AB produced")
ax[0][0].plot(range(max_nr_F+1),[x/float(sum_repl) for x in lrepl],label="Replication")

# # d_genomes[gnm]={"ab_prod": dg[key]["A"],  "howmany": 1 , "howmany_offs":len(dg[key]["offs_tag"]) , "earliest": dg[key]["birth"] }
# l_gnm_ab_nrg = (max_nr_F+1)*[0]
# l_gnm_howmany_nrg = (max_nr_F+1)*[0]
# l_ab_vs_howmany = []
# for genome in d_genomes:
#     l_ab_vs_howmany.extend( [[genome.count("F"), math.log10( d_genomes[genome]["ab_prod"] ) if d_genomes[genome]["ab_prod"]!=0 else 0 ] for _ in range(d_genomes[genome]["howmany"])] )
#     #nr_g = genome.count("F")
#     #l_gnm_ab_nrg[nr_g] += d_genomes[genome]["ab_prod"]
#     #l_gnm_howmany_nrg[nr_g] += d_genomes[genome]["howmany"]

# l_ab_vs_howmany = [*zip(*l_ab_vs_howmany)]
# H,xedges,yedges = np.histogram2d( l_ab_vs_howmany[1],l_ab_vs_howmany[0], bins=(20, range(-1,20,1)) )
# xx, yy = np.meshgrid(xedges, yedges)
# cmap = copy.copy(mpl.cm.get_cmap("gist_earth_r"))
# pcm = ax[0][0].pcolormesh(xx[:-1], yy[:-1], H.T,cmap=cmap)

# ax[0][0].scatter(l_ab_vs_howmany[1],l_ab_vs_howmany[0], s= l_ab_vs_howmany[2], alpha=0.2)
# ax[0][0].set_xscale('log')
# ax[0][0].set_yscale('log')

#ax[0][0].plot(range(max_nr_F+1), l_gnm_ab_nrg ,label="AB produced gnm")
#ax[0][0].plot(range(max_nr_F+1), l_gnm_howmany_nrg ,label="how many gnm")

# non-normalised
# ax[0][0].plot(range(max_nr_F+1),[x for x in lab],label="AB produced")
# ax[0][0].plot(range(max_nr_F+1),[x for x in lrepl],label="Replication")

ax[0][0].set_xlabel("Nr. of Growth-promoting genes")
ax[0][0].legend()

ax[0][1].plot(bin_edgesdf[:-1],hist_df, label="density")
ax[0][1].set_xlabel("Fraction of gr-promoting genes after mutations")
# ax[1].set_yscale('log')
# print("Here1")
# ax[1].plot(bin_edgesdf[:-1],lmedian_deltaf2, label="median density 2")
# ax[1].plot(bin_edgesdf[:-1],lmean_deltaf2, label="mean density 2")
# ax[1].plot(bin_edgesdf[:-1],l25per_deltaf2, label="25perc density 2")
# ax[1].plot(bin_edgesdf[:-1],l75per_deltaf2, label="75 per 2")
# print("Here2")
ax[0][1].legend()

ax[0][1].fill_between(bin_edgesdf[:-1],hist_df,[0. for x in bin_edgesdf[:-1]], alpha=0.3)
ax[0][1].plot( [ 0 , 1 ] , [0,0], lw=0.5, c='red', linestyle = 'dotted' ) #a zero line to facilitate comparison
ax[0][1].set_xlim([0.,1.])


# print("Making Hist")
# H, xedges, yedges = np.histogram2d(ldeltaf_po, ldeltaab, bins=(range(20), range(-10,10,1)) )
# print("Transposing")
# H = H.T
# print("Meshgrid")
# xx, yy = np.meshgrid(xedges, yedges)
# # mask some 'bad' data, in your case you would have: data == 0
# H = np.ma.masked_where(H < 0.0000001, H)
# # cmap = plt.cm.magma_r
# cmap = copy.copy(mpl.cm.get_cmap("magma_r"))
# cmap.set_bad(color='white')
# print("plot")
# pcm = ax[2].pcolormesh(xx, yy, H,cmap=cmap)

# ax[2].plot(range(1+max(ldeltaf_po)), lavr_deltaab )
# ax[2].fill_between(range(1+max(ldeltaf_po)), [x+y for x,y in zip(lavr_deltaab,lstd_deltaab)],[x-y for x,y in zip(lavr_deltaab,lstd_deltaab)],alpha=0.3 )
# ax[2].plot(range(1+max(ldeltaf_po)), [x-y for x,y in zip(lavr_deltaab,lstd_deltaab)] )
ax[1][0].plot(range(1+max(ldeltaf_po)), lmedian_deltaab, label=r'$\Delta$'+" antibiotic production")
ax[1][0].fill_between(range(1+max(ldeltaf_po)), l25perc_deltaab, l75perc_deltaab, alpha=0.3 )
ax[1][0].plot( [ 0 , max(ldeltaf_po) ] , [0,0], lw=0.5, c='red', linestyle = 'dotted' ) #a zero line to facilitate comparison
ax[1][0].set_xlabel(r'$\Delta$'+" nr. growth-promoting genes")
ax[1][0].legend()
# fig.colorbar(pcm)
# c = ax.imshow(H,interpolation ='none', origin ='lower') 
# plt.tight_layout()
print("l_repr_ab_for_lollyplot",l_repr_ab_for_lollyplot)
#makes sense to make lollyplots only if we got founders, otherwise it's a bit messy
if len(sys.argv)>2:
    #lollyplot
    # l_repr_ab_for_lollyplot= sorted( l_repr_ab_for_lollyplot, key=lambda x: (x[1],-x[0]),reverse=True ) 
    l_repr_ab_for_lollyplot=l_repr_ab_for_lollyplot[::-1]
    my_range = range( len(l_repr_ab_for_lollyplot) )
    I_want_log10 = False
    I_want_fractions = False
    I_want_absolute_numbers=False
    I_am_hacking_stuff = False
    I_want_scatterplot = False
    I_want_separate_plots = True
    if I_want_log10:
        ax[1][1].hlines(y=my_range, xmin=[ -math.log10(x[0]) if x[0]>0 else 0. for x in l_repr_ab_for_lollyplot ], 
                                                    xmax=[  math.log10(x[1]) if x[1]>0 else 0. for x in l_repr_ab_for_lollyplot ], 
                                                    color='grey', alpha=.8)
        ax[1][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][1].scatter([ -math.log10(x[0]) if x[0]>0 else 0. for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[1][1].scatter([  math.log10(x[1]) if x[1]>0 else 0. for x in l_repr_ab_for_lollyplot ] , my_range, color='orange', alpha=1. , label='AB')
        ax[1][1].legend()
    elif I_want_absolute_numbers:
        ax[1][1].hlines(y=my_range, xmin=[ -x[0] for x in l_repr_ab_for_lollyplot ], 
                                    xmax=[  x[1] for x in l_repr_ab_for_lollyplot ], 
                                    color='grey', alpha=.8)
        ax[1][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][1].scatter([ -x[0] for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[1][1].scatter([  x[1] for x in l_repr_ab_for_lollyplot ], my_range, color='orange', alpha=1. , label='AB')
        ax[1][1].legend()
    elif I_want_fractions:
        sum_l_repr_ab_for_lollyplot_0 = sum(x[0] for x in l_repr_ab_for_lollyplot)
        sum_l_repr_ab_for_lollyplot_1 = sum(x[1] for x in l_repr_ab_for_lollyplot)
        l_repr_ab_for_lollyplot = [ [x[0]/float(sum_l_repr_ab_for_lollyplot_0), x[1]/float(sum_l_repr_ab_for_lollyplot_1) ] for x in l_repr_ab_for_lollyplot ]
        ax[1][1].hlines(y=my_range, xmin=[ -x[0] for x in l_repr_ab_for_lollyplot ], 
                                    xmax=[  x[1] for x in l_repr_ab_for_lollyplot ], 
                                    color='grey', alpha=.8)
        ax[1][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][1].scatter([ -x[0] for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[1][1].scatter([  x[1] for x in l_repr_ab_for_lollyplot ], my_range, color='orange', alpha=1. , label='AB')
        ax[1][1].legend()
    elif I_want_scatterplot:
        sum_l_repr_ab_for_lollyplot_0 = sum(x[0] for x in l_repr_ab_for_lollyplot)
        sum_l_repr_ab_for_lollyplot_1 = sum(x[1] for x in l_repr_ab_for_lollyplot)
        l_repr_ab_for_lollyplot = [ [x[0]/float(sum_l_repr_ab_for_lollyplot_0), x[1]/float(sum_l_repr_ab_for_lollyplot_1),x[2] ] for x in l_repr_ab_for_lollyplot ]
        
        # ax[1][1].scatter([ x[1] if x[1]>0 else 0. for x in l_repr_ab_for_lollyplot ], 
        #                 [ x[0] if x[0]>0 else 0. for x in l_repr_ab_for_lollyplot ] )
        # sorted_to_plot = l_repr_ab_for_lollyplot.sort(key=lambda x: x[1]) # sort on AB
        bla = ax[1][1].scatter([ x[0] for x in l_repr_ab_for_lollyplot ], 
                        [ x[1] for x in l_repr_ab_for_lollyplot ])
        print("l_repr_ab_for_lollyplot",l_repr_ab_for_lollyplot)
        for i,x in enumerate(l_repr_ab_for_lollyplot):
            print("x",x)
            ax[1][1].annotate(x[2], (x[0],x[1]),fontsize=2 ,rotation=45)
        # ax[1][1].invert_xaxis()
        # ax[1][1].invert_yaxis()
        ax[1][1].set_ylabel("AB")
        ax[1][1].set_xlabel("repl")
    if I_want_separate_plots:
        #we makes two separate plots in ax[1][1] and ax[1][2]
        sum_l_repr_ab_for_lollyplot_0 = sum(x[0] for x in l_repr_ab_for_lollyplot)
        sum_l_repr_ab_for_lollyplot_1 = sum(x[1] for x in l_repr_ab_for_lollyplot)
        l_repr_ab_for_lollyplot = [ [x[0]/float(sum_l_repr_ab_for_lollyplot_0), x[1]/float(sum_l_repr_ab_for_lollyplot_1) ] for x in l_repr_ab_for_lollyplot ]
        ax[1][1].hlines(y=my_range, xmin=[ 0 for x in l_repr_ab_for_lollyplot ], 
                                    xmax=[  x[0] for x in l_repr_ab_for_lollyplot ], 
                                    color='grey', alpha=.8)
        ax[1][2].hlines(y=my_range, xmin=[ 0 for x in l_repr_ab_for_lollyplot ], 
                                    xmax=[  x[1] for x in l_repr_ab_for_lollyplot ], 
                                    color='grey', alpha=.8)
                                    
        ax[1][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][2].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][1].scatter([ x[0] for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[1][2].scatter([ x[1] for x in l_repr_ab_for_lollyplot ], my_range, color='orange', alpha=1. , label='AB')
        ax[1][1].legend()
        ax[1][2].legend()
    else:
        print()
        print("For lollyplot choose one between I_want_log10,I_want_absolute,I_want_fractions")
        ax[1][1].text(0,0, "For lollyplot choose one between I_want_log10,I_want_absolute,I_want_fractions\n you chose nothing")
        print()
    if I_am_hacking_stuff:
        print("WARNING: this is just a little test, switch me off later!")
        ax[0][1].clear()
        ax[1][1].clear()
        ax[0][1].hlines(y=my_range, xmin=[ -x[0] for x in l_repr_ab_for_lollyplot ], 
                                    xmax=[  x[1] for x in l_repr_ab_for_lollyplot ], 
                                    color='grey', alpha=.8)
        ax[0][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[0][1].scatter([ -x[0] for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[0][1].scatter([  x[1] for x in l_repr_ab_for_lollyplot ], my_range, color='orange', alpha=1. , label='AB')
        ax[0][1].legend()
    
        sum_l_repr_ab_for_lollyplot_0 = sum(x[0] for x in l_repr_ab_for_lollyplot)
        sum_l_repr_ab_for_lollyplot_1 = sum(x[1] for x in l_repr_ab_for_lollyplot)
        print("sum_l_repr_ab_for_lollyplot_0,sum_l_repr_ab_for_lollyplot_1",sum_l_repr_ab_for_lollyplot_0,sum_l_repr_ab_for_lollyplot_1)
        print( [x[0] for x in l_repr_ab_for_lollyplot] )
        print( [x[1] for x in l_repr_ab_for_lollyplot] )
        # sys.exit(1)
        ax[0][1].clear()
        ax[1][1].clear()
        l_repr_ab_for_lollyplot = [ [x[0]/float(sum_l_repr_ab_for_lollyplot_0), x[1]/float(sum_l_repr_ab_for_lollyplot_1) ] for x in l_repr_ab_for_lollyplot ]
        # ax[1][1].hlines(y=my_range, xmin=[ -x[0] for x in l_repr_ab_for_lollyplot ], 
                                    # xmax=[  x[1] for x in l_repr_ab_for_lollyplot ], 
                                    # color='grey', alpha=.8)
        ax[1][1].plot([0,0],[min(my_range),max(my_range)], lw=0.5, c='red', linestyle = 'dotted')
        ax[1][1].scatter([ -x[0] for x in l_repr_ab_for_lollyplot ], my_range, color='blue', alpha=1., label='repl')
        ax[0][1].scatter([  x[1] for x in l_repr_ab_for_lollyplot ], my_range, color='orange', alpha=1. , label='AB')
        ax[1][1].legend()
        plt.show()
else:
    #Collect anyone that is alive now, sort by howmany ABs they made
    #then make cumulative plot
    l_ab_prod=[]
    l_ab_predicted=[]
    l_hist =[]
    for key in dg:
        if dg[key]['death'] >= maxtime:
            #this guy is alive now
            ab_prod = dg[key]['A']
            #birthdate = dg[key]['birth']
            #l_ab_prod.append(ab_prod/float(maxtime-birthdate))
            l_ab_prod.append(ab_prod)
            #gnm = dg[key]['gnm']
            #nA = gnm.count('A') 
            #nF = gnm.count('F')
            #l_hist.append(nF)
            #l_ab_predicted.append( nA/(nA + 3.)*math.exp(-nF) )
    tot_people_here = len(l_ab_prod)
    print("tot people: ", tot_people_here)
    sum_ABs_here = sum(l_ab_prod)
    print("Sum ABs = ", sum_ABs_here)
    sorted_l_ab_prod = sorted(l_ab_prod, reverse=True)
    # print(len(l_ab_prod), sorted_l_ab_prod[0],sorted_l_ab_prod[-1] )
    from itertools import accumulate
    l_sum_sorted = list( accumulate(sorted_l_ab_prod) ) # cumulative sum A LOT faster than what I wrote!
    X = [x/float(tot_people_here) for x in range(tot_people_here)]
    Y = [y/float(sum_ABs_here) for y in l_sum_sorted]
    ax[1][1].plot( X , Y , label = 'cumulative antibiotic production' )
    
    # Attempt at log transforming the data
    # data is in [0,1], mapped to -log10(1-y)
    # result is nice and clear, but not really suited for anyting
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    #Create an inset in the lower right corner (loc=4) with borderpad=1, i.e.
    # 10 points padding (as 10pt is the default fontsize) to the parent axes
    axins = inset_axes(ax[1][1], width="100%", height="100%", bbox_to_anchor=(.6, .2, .4, .4), bbox_transform=ax[1][1].transAxes)

    lX=[]
    lY=[]
    for x,y in zip(X,Y):
        if y<0.999:
            try:
                lY.append(-math.log10(1.-y))
                lX.append(x)
            except:
                break
    axins.plot( lX , lY )

    for x,y in zip(X,Y):
        if y>=0.99:
            print("HEllo x,y", x,y)
            ax[1][1].scatter([x],[y],s=20, c='r')
            ax[1][1].plot([0,x,x],[y,y,0],ls='--', lw=0.5, c='r', label='99%')
            try:
                logy = -math.log10(1.-y)
                axins.scatter([x],[logy],s=10,c='r')
                axins.plot([0,x,x],[logy,logy,0],ls='--', lw=0.5, c='r')
            except ValueError:
                print("Math domain error, not plotting inset")
            break
    
    ax[1][1].legend()
    ax[1][1].set_xlim([0,1])
    ax[1][1].set_ylim([0,1])
    ax[1][1].set_xlabel('Cumulative fraction of population')
    ax[1][1].set_ylabel('Fraction of total antibiotic production')
    
    axins_yaxis = [0., 0.9,0.99,0.999]
    axins.set_yticks([-math.log10(1.-x) for x in axins_yaxis])
    axins.set_yticklabels([str(x) for x in axins_yaxis])
    # axins.set_y ([str(x) for x in axins_yaxis])
    
    #ax[1][1].plot( [x/float(tot_people_here) for x in range(tot_people_here)] ,[y/float(sum_ABs_here) for y in sorted_l_ab_prod])
    # ax[1][1].set_yscale('log')
    # mybins = np.linspace(0,1,21)
    # hist,bins  = np.histogram(l_ab_prod)
    # hist2,bins = np.histogram(l_ab_predicted, bins=mybins)
    
    # ax[1][1].plot(bins[:-1],hist, lw=1)
    # ax[1][1].plot(bins[:-1],hist2, lw=1, label='hist Ab predicted')
    #ax[1][1].plot( [ 0 , bins[-1] ] , [0,0], lw=0.1, c='red', linestyle = 'dotted' ) #a zero line to facilitate comparison
plt.show()
