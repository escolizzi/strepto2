#!/usr/bin/python2.7

'''
This script attempts to find how many AB clusters are really there
'''

import sys,math,os,string,random
from subprocess import Popen, PIPE
#from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec
import numpy as np

def Birthrate(a,b):
    hd = bin(a^b).count('1')
    return math.exp(-0.3*hd**2.)

ltime=[]
lF=[]
lA=[]
lB=[]
lav_F=[]
lav_A=[]
lav_B=[]
lstd_F=[]
lstd_A=[]
lstd_B=[]
l_antib=[]
lab_resistance=[]
l_l_susceptibility=[]
l_num_antib=[]
l_avr=[]

ab_list_pos_in_file = 2
genome_pos_in_file = 4
gene_ab_pos_in_file = 5

filename = sys.argv[1]
maxtime =sys.maxint 
try :
    maxtime = int(sys.argv[2])
except:
    pass

if os.path.isfile(filename+"_avrgsuscept_intime.txt"):
    print "Warning: output file exists, re-calculate and overwrite? y/n?"
    answer =  sys.stdin.read(1)
    if(answer=='y'): 
        print "Alright, let's go!"
    else:
        print "Not go, program exits now, bye!"
        exit(0)

time = int( (Popen(["head","-n1",filename], stdout=PIPE).stdout.read()).split()[0])
print "Initial time =",time

i=0
with open(filename,"r") as fin:
    for line in fin:
        line =line.split()
        timenow=int(line[0])
        if timenow % 10000 != 0: continue
        
        if timenow != time:
            if time>maxtime: break
            if time%1000==0: 
                print
                print "Time:",time

            ltime.append(time)
            time = timenow

            # Make a set from l_antib, you need antib diversity
            l_antib = list(set(l_antib))
            print "So many antibiotics: ", len(l_antib)
            l_num_antib.append(len(l_antib))
            # compare against lab_resistance, which has as many elements as there are individuals
            
            #TO SAMPLE UNCOMMENT THE FOLLOWING LINE. It's like 30000 bact -> 1000 will do
            # lab_resistance = random.sample(lab_resistance,min(len(lab_resistance),100))
            l_susceptibility=len(lab_resistance)*[0.] #contains the cum susceptibility to antib
            #loop over all (sampled) bacteria
            for i,icel in enumerate(lab_resistance):
                #for each bact, check it min birthrate to each ab in the field
                #so go through all antib in the field
                for field_ab in l_antib:
                    l_icel_birthrates=[Birthrate( int(i_ab),int(field_ab) ) for i_ab in icel] # This has to be kept in sync with simulations
                    # now l_icel_birthrates contains a bunch of growth rates, for each AB in the field.
                    # given the AB genes in the bacterium, it grows with the rate given by the max resistance for this particular AB
                    # so if this max growth rate is less than 0.25, then we say it is susceptible
                    if max(l_icel_birthrates) < 0.25 : l_susceptibility[i] += 1 
            #At this point l_susceptibility can be turned into a histogram    
            hist,bins = np.histogram( l_susceptibility, bins = np.linspace(0, 1000, num=101) )
            mean_susc = np.mean(l_susceptibility)
            l_avr.append( np.mean(l_susceptibility) )
            print "Mean susceptibility: ", mean_susc
            #print l_susceptibility
            l_l_susceptibility.append(hist)

            # lav_F.append(np.mean(lF))
            # lstd_F.append(np.std(lF))
            # lav_A.append(np.mean(lA))
            # lstd_A.append(np.std(lA))
            # lav_B.append(np.mean(lB))
            # lstd_B.append(np.std(lB))

            lF=[];lA=[];lB=[];l_antib=[];lab_resistance=[];
            print "Next time step"
        try:
            antib_at_pos = line[ab_list_pos_in_file]
            if antib_at_pos != '0,':
                l_antib.extend( antib_at_pos.split(',')[:-1] ) #get ab deposited in the field
            # if genome, and if there is something after
            if line[genome_pos_in_file] != 'n' and len(line)==1+gene_ab_pos_in_file:
                ab_genome = list( set(line[gene_ab_pos_in_file].split(',')[:-1]))
                lab_resistance.append( ab_genome )

            # if l_antib and lab_resistance:
            #     print l_antib, lab_resistance
            #     #sys.exit(1)
            # genome=line[genome_pos_in_file]
            # lF.append(genome.count("F"))
            # lA.append(genome.count("A"))
            # lB.append(genome.count("B"))
        except Exception as e: 
            print(e)
            print "Faulty line number", i
            print line
            #sys.exit(1)
            pass
        i+=1

fig, (ax0, ax1,ax2) = plt.subplots(3,sharex=True)

#tot number AB
ax0.plot([int(x)+500 for x in ltime] , l_num_antib )
ax0.set_title("Total nr. of ABs in the field")

# Distribution of nr. of AB to which indiv are susc.
a = np.ma.masked_where(l_l_susceptibility < 0.000000001, l_l_susceptibility)
#ax0.pcolor(ltime, bins[:-1], [x[::-1] for x in zip(*a)])
ax1.pcolor(ltime, bins[:-1], zip(*a))
ax1.set_title("Distr of how many bact are susc to how many ABs")
# print l_num_antib
#plt.plot( [int(x)+500 for x in ltime] , l_num_antib)

# avrg. Fraction of AB susc.
ax2.plot( [int(x)+500 for x in ltime] , [ x/float(y) if y!=0 else 0 for x,y in zip(l_avr,l_num_antib) ] )
ax2.set_ylim(0,1)
ax2.set_title("Avrg. fraction of AB susceptibility")
ax2.set_xlabel("Time")
# plt.plot(ltime,lav_F,label = "avrg #F")
# plt.plot(ltime,lav_A,label = "avrg #A")
# plt.plot(ltime,lav_B,label = "avrg #B")

# plt.fill_between(ltime, [x-y for x,y  in zip(lav_F,lstd_F)], [x+y for x,y  in zip(lav_F,lstd_F)],alpha=0.5)
# plt.fill_between(ltime, [x-y for x,y  in zip(lav_A,lstd_A)], [x+y for x,y  in zip(lav_A,lstd_A)],alpha=0.5)
# plt.fill_between(ltime, [x-y for x,y  in zip(lav_B,lstd_B)], [x+y for x,y  in zip(lav_B,lstd_B)],alpha=0.5)
# plt.ylim(ymin=0)
# plt.legend()
plt.margins(0,0)

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
