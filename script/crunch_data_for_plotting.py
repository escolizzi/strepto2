#!/usr/bin/python3

'''
this version only outputs the time series as

T1 Fmin Fmed Fmax Amin Amed Amax Gmin Gmed Gmax
T2 Fmin Fmed Fmax Amin Amed Amax Gmin Gmed Gmax
...

'''


import sys,math,os,string
from subprocess import Popen, PIPE
import numpy as np
from scipy.stats.mstats import mquantiles
from pathlib import Path


lF=[]
lA=[]
lB=[]

l_data_to_save=[]

genome_pos_in_file = 4

filename = sys.argv[1]

lfilename = filename.split("/")
lfilename[-1] = "crunched_"+lfilename[-1]
fileoutputname = "/".join(lfilename)

# Check if output file already exists
print("file output is: ",fileoutputname)
my_file = Path(fileoutputname)
if my_file.is_file():
    start_time = int((Popen(["tail","-n1",my_file], stdout=PIPE).stdout.read()).split()[0])
    time_from_infile = int((Popen(["tail","-n1",filename], stdout=PIPE).stdout.read()).split()[0])
    if time_from_infile == start_time: 
        print("Crunched file is up to date. No crunching! :)")
        sys.exit(1)
    print("Crunched file exists - appending! Last time crunched file: ", start_time, ", last time data file: ",time_from_infile)
    how_to_open_fileout="a"
else:
    # time = (Popen(["head","-n1",filename], stdout=PIPE).stdout.read()).split()[0]
    start_time=-1
    how_to_open_fileout="w"

print("Initial time =",start_time)

maxtime =float("inf") 
try :
    maxtime = int(sys.argv[2])
except:
    pass

first_pass=True

with open(filename,"r") as fin:
    for line in fin:
        line =line.split()
        timenow=int(line[0])
        
        if timenow<=start_time: continue
        
        if first_pass == True:
            time=timenow
            first_pass = False
        
        if len(line)<4: continue
        if "n" in line: continue
        
        if timenow != time:
            if time>maxtime: break
            if time%100000==0: print("Time: ",time)
            
            Fq25,Fmed,Fq75=mquantiles(lF)
            Aq25,Amed,Aq75=mquantiles(lA)
            Bq25,Bmed,Bq75=mquantiles(lB)
            
            l_data_to_save.append([time,Fq25,Fmed,Fq75,Aq25,Amed,Aq75,Bq25,Bmed,Bq75])
            
            time = timenow
            lF=[];lA=[];lB=[]
            
        try:
            genome=line[genome_pos_in_file]
            lF.append(genome.count("F"))
            lA.append(genome.count("A"))
            lB.append(genome.count("B"))
        except:
            pass

with open(fileoutputname,how_to_open_fileout) as fout:
    for line in l_data_to_save:
        for item in line:
            fout.write(str(item)+" ")
        fout.write("\n")
print("Done")