#!/usr/bin/python3

'''
This is hopefully a serious speed up from plot_ab_susc...py
We are doing things per-colony, so no single genome considered.
'''
import matplotlib as mpl
mpl.use('qt5agg')
import sys,math,os,string,random
from subprocess import Popen, PIPE
#from PIL import Image
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
from itertools import chain
from collections import Counter
import numpy as np

def Birthrate(a,b):
    hd = bin(a^b).count('1')
    return math.exp(-hd**2.)

ltime=[]

l_antib=[]
lab_resistance=[]
l_l_susceptibility=[]
l_l_howmany_unique_ab_genes=[]
l_l_howmany_unique_ab_genes_actually_used=[]
l_num_antib=[]
l_avr=[]
l_percentiles=[]
l_median=[]
l_avr_pairwise_susc=[]

dpop={}
l_nr_colonies = []

ab_list_pos_in_file = 2
genome_pos_in_file = 4
gene_ab_pos_in_file = 5

filename = sys.argv[1]
maxtime =float("inf")
if len(sys.argv)>2:
    try :
        maxtime = int( sys.argv[2] )
    except:
        print("Warning: did not get maxtime (argv[2])")
        sys.exit(1)

if os.path.isfile("bla_avrgsuscept_intime.txt"):
    print ("Warning: output file exists, re-calculate and overwrite? y/n?")
    answer =  sys.stdin.read(1)
    if(answer=='y'): 
        print ("Alright, let's go!")
    else:
        print ("Not go, program exits now, bye!")
        exit(0)

time = int( (Popen(["head","-n1",filename], stdout=PIPE).stdout.read()).split()[0])
print ("Initial time =",time)

i=0
with open(filename,"r") as fin:
    for line in fin:
        line =line.split()
        timenow=int(line[0])
        if timenow % 50000 != 0: continue
        # if timenow>maxtime: break
        
        if timenow != time:
            if time>maxtime: break
            if time%100000==0: 
                print()
                print( "Time:",time)

            ltime.append(time)
            time = timenow

            # Make a set from l_antib, you need antib diversity
            l_antib = list(set(l_antib))
            print ("So many different antibiotics: ", len(l_antib))
            l_num_antib.append(len(l_antib))
            # compare against lab_resistance, which has as many elements as there are individuals
            
            #TO SAMPLE UNCOMMENT THE FOLLOWING LINE. It's like 30000 bact -> 1000 will do
            # lab_resistance = random.sample(lab_resistance,min(len(lab_resistance),100))
            #l_susceptibility=len(lab_resistance)*[0.] #contains the cum susceptibility to antib
            #loop over all (sampled) bacteria
            # for i,icel in enumerate(lab_resistance):
            l_susceptibility=[]
            sum_pairwise_susc=0.
            
            lsorted_keys=sorted( dpop.keys() ) # this is too make sure that iteration is always in order
                                               # which is guaranteed in python 3.7+ but... well
            print("Got ",len(lsorted_keys)," different colonies")
            for counter1,colony1 in enumerate(lsorted_keys):
                set_dcolony1 = set( dpop[colony1] )

                for colony2 in lsorted_keys[counter1+1:]:
                    set_dcolony2 = set(dpop[colony2])
                    flen_union = float( len( ( set_dcolony1 | set_dcolony2 ) ) )
                    # sum_pairwise_susc += len( (set_dcolony1 - set_dcolony2 ) ) + len( (set_dcolony2 - set_dcolony1 )  ) / len_union
                    # the above is the sum of the differences, which is also the difference between union and intersection, 
                    # which has the advantage that we can do one operation less
                    sum_pairwise_susc += 1. - len( ( set_dcolony1 & set_dcolony2 ) )/flen_union

            l_avr_pairwise_susc.append(   sum_pairwise_susc / float( (len(dpop))**2. - len(dpop) )   )
            #actually looping over colonies
            n_to_exclude=2
            for dcounter, colony in enumerate(dpop):
                counted_abs = Counter(dpop[colony])
                set_abs = set(dpop[colony])
                
                dpop[colony] = [ x for x in set_abs if counted_abs[x]>n_to_exclude]
                #dpop[colony]=list(set(dpop[colony]))  # we should exclude rare genes, because they are not representative - but how many are there?
                                                # FOR NOW EXCLUDE NOBODY, THEN WE'LL SEE
                #print("dpop[colony]:",dpop[colony])
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
            
            quantiles = np.quantile(l_susceptibility, [0.05,0.25,0.5,0.75,0.95] )
            l_percentiles.append(quantiles)
            l_nr_colonies.append(len(dpop))
            #find unique ab genes, per colony
            l_union_ab_genes =   [ [ dpop[key] ,  chain(*[ dpop[key2] for key2 in dpop if key2!=key ]) ] for key in dpop ]  #unionize everybody except the list at that position
            l_unique_ab_genes = [ set(x[0]) - set(x[1]) for x in l_union_ab_genes]
            l_howmany_unique_ab_genes = [ len(x) for x in l_unique_ab_genes ]
            l_l_howmany_unique_ab_genes.append(l_howmany_unique_ab_genes)
            s_antib= set(l_antib)
            l_howmany_unique_ab_genes_actually_used = [ len(x-s_antib) for x in l_unique_ab_genes ]
            l_l_howmany_unique_ab_genes_actually_used.append(l_howmany_unique_ab_genes_actually_used)
            print("Time = ",ltime[-1]," mean(unique ab) = ", np.mean(l_howmany_unique_ab_genes), 
                      ", mean actually used: ", np.mean(l_howmany_unique_ab_genes_actually_used) )

            # print("Quantiles", quantiles)

            # #At this point l_susceptibility can be turned into a histogram    
            # hist,bins = np.histogram( l_susceptibility, bins = np.linspace(0, 1000, num=101) )
            # mean_susc = np.mean(l_susceptibility)
            # l_avr.append( np.mean(l_susceptibility) )
            # print ("Mean susceptibility: ", mean_susc)
            # #print l_susceptibility
            # l_l_susceptibility.append(hist)


            lF=[];lA=[];lB=[];l_antib=[];lab_resistance=[];
            dpop={}
            print ("Next time step, it's now: ", time,"\n")
        try:
            # Antibiotics in the field
            antib_at_pos = line[ab_list_pos_in_file]
            if antib_at_pos != '0,':
                l_antib.extend( antib_at_pos.split(',')[:-1] ) #get ab deposited in the field, remove the last bc it's empty after last comma
            # Genome, and if there is something after
            if line[genome_pos_in_file] != 'n' and len(line)>genome_pos_in_file and line[gene_ab_pos_in_file]!=',':
                ab_genome = list( set( line[gene_ab_pos_in_file].split(',')[:-1] ) )
                
                if len(ab_genome):
                    
                    colony_val = int(line[3])
                    # this is the colony antibiotic pan-genome
                    # we will exclude the very rare abs later - if they - happen, because they are not representative of the colony
                    if colony_val in dpop:
                        dpop[colony_val].extend(ab_genome)
                        #print("Extended: dpop[colony_val] = ",dpop[colony_val])
                    else:
                        dpop[colony_val] = ab_genome
        except Exception as e: 
            print(e)
            print ("Faulty line number", i)
            print (line)
            #sys.exit(1)
            pass
        i+=1


# PLOTTING 
print()


#fig, (ax0, ax1,ax2) = plt.subplots(3,sharex=True)
fig, (ax2,ax3) = plt.subplots(nrows=2,ncols=1,sharex=True)
ax0 = ax2.twinx()
#tot number AB

# Distribution of nr. of AB to which indiv are susc.
# l_l_susceptibility = np.array(l_l_susceptibility)
# print("Shape: ",l_l_susceptibility.shape)
# a = np.ma.masked_where(l_l_susceptibility < 0.000000001, l_l_susceptibility)
#ax0.pcolor(ltime, bins[:-1], [x[::-1] for x in zip(*a)])
#ax1.pcolormesh(ltime, bins[:-1], zip(*a))

# ax1.pcolormesh(ltime, bins[:-1], a.T)
# ax1.set_title("Distr of how many bact are susc to how many ABs")
# print l_num_antib
#plt.plot( [int(x)+500 for x in ltime] , l_num_antib)


l_percentiles = np.array(l_percentiles).T
#ltime_forplotting = [int(x)+500 for x in ltime] # not sure why +500
ltime_forplotting = [int(x) for x in ltime]
# print("len ltime_forplotting", len(ltime_forplotting),ltime_forplotting)
# print("len l_percentiles[1]",len(l_percentiles[1]),l_percentiles[1])
# print("len l_num_antib",len(l_num_antib),l_num_antib)

# print("Zip list: ",   [ [x,y] for x,y in zip(l_percentiles,l_num_antib) ] )
# print("Conditioned list: ",   [ x/float(y) if y!=0 else 0 for x,y in zip(l_percentiles[2],l_num_antib) ] )

ax2.plot( ltime_forplotting , [ x/float(y) if y!=0 else 0 for x,y in zip(l_percentiles[2],l_num_antib) ],color = u'#ee7f0e',zorder=10 )
ax2.fill_between( ltime_forplotting, 
                  [ x/float(y) if y!=0 else 0 for x,y in zip(l_percentiles[1],l_num_antib) ], 
                  [ x/float(y) if y!=0 else 0 for x,y in zip(l_percentiles[3],l_num_antib)] ,alpha=0.5, color = u'#ee7f0e',zorder=9)
ax2.fill_between( ltime_forplotting, 
                  [ x/float(y) if y!=0 else 0 for x,y in zip(l_percentiles[0],l_num_antib) ], 
                  [ x/float(y) if y!=0 else 0 for x,y in zip(l_percentiles[4],l_num_antib)] ,alpha=0.2,color = u'#ee7f0e',zorder=8)
ax2.set_ylabel('susceptibility', color= u'#ee7f0e')
ax2.tick_params(axis='y', labelcolor = u'#ee7f0e')

#plot( [int(x)+500 for x in ltime] , [ x[1]/float(y) if y!=0 else 0 for x,y in zip(l_percentiles,l_num_antib) ]
#ax2.plot( [int(x)+500 for x in ltime] , [ x[1]/float(y) if y!=0 else 0 for x,y in zip(l_percentiles,l_num_antib) ]

# # avrg. Fraction of AB susc.
# ax2.plot( [int(x)+500 for x in ltime] , [ x/float(y) if y!=0 else 0 for x,y in zip(l_avr,l_num_antib) ] )
ax2.set_ylim(0,1)
# ax2.set_title("Avrg. fraction of AB susceptibility")
# ax2.set_xlabel("Time")

ax0.plot([int(x) for x in ltime] , l_num_antib ,color='black', alpha=0.5,lw=3,zorder=0)
ax0.plot([int(x) for x in ltime] , l_nr_colonies, ls='--' ,color='black', alpha=0.5,lw=3,zorder=0)
# ax0.set_title("Total nr. of ABs in the field", color='dimgrey')
ax0.set_ylabel("Total nr. of ABs in the field", color='dimgrey')
ax0.tick_params(axis='y', labelcolor='dimgrey')
ax0.set_ylim(0,None)

plt.margins(0,0)

print("Histogramming unique genes and actually used")
q1=[]
q2=[]
# print("Lengths")
# print(len(l_howmany_unique_ab_genes))
# print(len(l_howmany_unique_ab_genes_actually_used))
# print("Arrays:")
# print(l_howmany_unique_ab_genes)
# print(l_howmany_unique_ab_genes_actually_used)
# print("End arrays;")

for x,y, in zip(l_l_howmany_unique_ab_genes, l_l_howmany_unique_ab_genes_actually_used):
    bins = range(100)
    q1.append( np.quantile(x,[0.05, 0.25, 0.5, 0.75, 0.95]) )
    q2.append( np.quantile(y,[0.05, 0.25, 0.5, 0.75, 0.95]) )

q1T=np.array(q1).T
q2T=np.array(q2).T
ax3.plot(ltime_forplotting,q1T[2],c='firebrick', label='unique antib genes')
ax3.fill_between(ltime_forplotting,q1T[1],q1T[3],color='firebrick',alpha=0.5)
ax3.fill_between(ltime_forplotting,q1T[0],q1T[4],color='firebrick',alpha=0.2)

ax3.plot(ltime_forplotting,q2T[2],c='royalblue', label='unique antib genes actualy used')
ax3.fill_between(ltime_forplotting,q2T[1],q2T[3],color='royalblue',alpha=0.5)
ax3.fill_between(ltime_forplotting,q2T[0],q2T[4],color='royalblue',alpha=0.2)

ax3.plot( ltime_forplotting, l_avr_pairwise_susc, label='Avr pairwise susc' )

ax3.set_ylim([0,None])
ax3.set_xlim([0,None])
ax3.legend()

if False:
    with open(filename+"_avrgsuscept_intime.txt", "w") as fout:
        fout.write("### Tot. AB in the field\n")
        for i,j in zip( ltime , l_num_antib ):
            fout.write(str(i)+" "+str(j)+"\n")
        
        fout.write("### Distr of how many bact are susc to how many ABs - [begin,end,step]: [0,1000,10]\n")
        for i,lnab in zip( ltime , l_l_susceptibility ):
            fout.write(str(i))
            for x in lnab: 
                fout.write(" "+str(x))
            fout.write("\n")

        fout.write("### Avrg. fraction of AB susceptibiltiy\n")
        for i,j in zip( ltime , [ x/float(y) if y!=0 else 0 for x,y in zip(l_avr,l_num_antib) ] ):
            fout.write(str(i)+" "+str(j)+"\n")

plt.show()
