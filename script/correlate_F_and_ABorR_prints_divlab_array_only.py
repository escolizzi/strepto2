#!/usr/bin/python3

'''
This is a much truncated version of correlate_F_and_AB...
It only prints the arrays of who makes replications and who does ab as a function of nr. of F

Makes plots of genome and genes distribution.
Positions must be specified in some way.. for now hard coded
-also added antib production
'''
# import matplotlib as mpl
# mpl.use('Qt5Agg')

import sys,math,os,string,random,copy
from subprocess import Popen, PIPE
from collections import Counter

#from PIL import Image
# mpl.use('GTK3Agg')
# from matplotlib import pyplot as plt
# import matplotlib.colors as colors
# import matplotlib.cm as cm
# from matplotlib.lines import Line2D
# import matplotlib.gridspec as gridspec
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
                    # print("Error, here is line: ", line)
                    #print("line_number = ", line_number)
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
        

#print("Done, reading data. In total so many people lived:", len(dg))

#print("Successful genomes - in order:")
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
    # if d_genomes[gnm]["howmany"]>50: #print(gnm,": ",d_genomes[gnm])
    if d_genomes[gnm]["earliest"]==0 and d_genomes[gnm]["howmany_offs"]>nr_offspring_to_count_as_successful:
        l_ancestor.append(gnm)
# Now go get these guys from sys.argv[2] = data_bla.txt
if len(sys.argv)<=2:
    I_am_doing_nothing_here=True
    #print("No ancestor genomes are searched")
else:
    print("Going to look for founder bacteria")
    print("This shouldn't happen here")
    sys.exit(1)
    with open(sys.argv[2],"r") as fin2:
        for line in fin2:
            line=line.split()
            if int(line[0]) != 0: break
            # if line[4] in l_ancestor:
            #     #print( line[4],line[5])

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

for ancestor in l_ancestor:
    # if d_genomes[ancestor]["howmany"]>2500:
        #print("Ancestor:") 
        #print(ancestor, d_genomes[ancestor]["howmany"],d_genomes[ancestor]["ab_prod"])
    l_repr_ab_for_lollyplot.append( [d_genomes[ancestor]["howmany"],d_genomes[ancestor]["ab_prod"], ancestor ] )
    
    # print("Offspring",) 
    # designs.append(MakeDesign(ancestor))
    tot_deltaab += d_genomes[ancestor]["ab_prod"]
    for offs_seq in sorted(d_gnm_and_offs[ancestor],reverse=True):
        if offs_seq == ancestor: continue
        # print(offs_seq,d_genomes[offs_seq]["howmany"],d_genomes[offs_seq]["ab_prod"])
        l_repr_ab_for_lollyplot.append( [d_genomes[offs_seq]["howmany"],d_genomes[offs_seq]["ab_prod"],offs_seq] )
        ldeltaf_po.append( ancestor.count("F") - offs_seq.count("F") )
        ldeltaab.append( d_genomes[ancestor]["ab_prod"] - d_genomes[offs_seq]["ab_prod"])
        tot_deltaab += d_genomes[offs_seq]["ab_prod"]
        # designs.append(MakeDesign(offs_seq))
        
        # design = [sp]
        # for x in offs_seq:
        #     design.append( l_drawgenes[l_genes.index(x)] )
        #     design.append(sp)
        # designs.append(design)
    ldeltaab=[x/float(tot_deltaab) for x in ldeltaab]
    # print("Genome: ", designs[-1])
# if print_genomes:
#     fig = plt.figure(figsize=(8,6))
#     gs = gridspec.GridSpec( len(designs), 2+ sum([ x['opts']['x_extent'] for x in designs[0] ]) )

#     ax_dna=[]
#     for i in range( len(designs) ):
#         if i==0:
#             ax_dna.append( plt.subplot( gs[i,1:] ) )
#         else:
#             ax_dna.append( plt.subplot( gs[i, 1:  sum([ x['opts']['x_extent'] for x in designs[i] ])   ] ) )
#         # ax_dna2 = plt.subplot(gs[1])

#         start, end = dr.renderDNA(ax_dna[i], designs[i], part_renderers)
#         #print (i, start,end)
#         ax_dna[i].set_xlim([start, end])
#         ax_dna[i].set_ylim([-15,28])
#         ax_dna[i].set_aspect('equal')
#         ax_dna[i].set_xticks([])
#         ax_dna[i].set_yticks([])
#         ax_dna[i].axis('off')

#     for i,my_ax in enumerate(ax_dna): 
#         set_size(len(designs[i]),5, ax=my_ax)
#     # plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, wspace=0.01)

#     # fig.tight_layout()
#     plt.savefig("pretty.pdf")
#     sys.exit(1)
# # else:
#     #lollyplots

# PARENT-OFFSPRING DIFFERENCES, AS A FUNCTION OF Delta(F)
max_nr_F=0

for key in dg:
    max_nr_F = max(dg[key]["gnm"].count("F") , max_nr_F)

# Make averages and std
max_deltaf=max(ldeltaf_po)
# print("max_deltaf parent-offspr: ",max_deltaf)
l1 = [[] for _ in range(1+max(ldeltaf_po))] # thanks past me - WTF is l1?
# print("l1: ",l1)
for df,ab in zip(ldeltaf_po,ldeltaab):
    l1[df].append(ab)
#Sometimes some positions in l1 are not filled


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
        
    if len(dg[key]["offs_tag"])>0:
        this_deltaf=[]
        for offs_tag in dg[key]["offs_tag"]:
            # print( dg[key]["gnm"] )
            # print( dg[offs_tag]["gnm"] )
            deltaFfract = ( dg[key]["gnm"].count("F") - dg[offs_tag]["gnm"].count("F") )/(float( dg[key]["gnm"].count("F") ))
            ldeltaF.append( deltaFfract )
            # this_deltaf.append( deltaFfract )
            
hist_df, bin_edgesdf = np.histogram( ldeltaF, bins = np.linspace(0.,1.05,num=21) )
hist_df = [x/float(len(dg)) for x in hist_df]

sum_lab=sum(lab)
sum_repl=sum(lrepl)

print("X = ",[x for x in range(max_nr_F+1)])
print("AB produced = ", [x/float(sum_lab) for x in lab])
print("Replication = ",[x/float(sum_repl) for x in lrepl])
