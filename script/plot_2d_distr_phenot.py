#!/usr/bin/python3

import sys,math
import matplotlib.pyplot as plt
import numpy as np

dc={}
tot_lines=0
tot_linesnoF=0
if len(sys.argv)<2:
  print("mutrate_percolony.py, calculates and plots the fraction of AB mutants *PER COLONY* in the system, feed it only one or few time steps")
filename=sys.argv[1]

l_abpr=[]
l_growth=[]

#filename="data_strpt_season1000_maxpAB0.01_betatradeoff3_bitstrABlen12_frnonproducing0._nomix_btwn_seasns_muttypeC.txt_t1000000.txt"
with open(filename, "r") as fin:
  for line in fin.readlines():
    line=line.split()
    if len(line)<5: continue
    if line[4]=='n' or line[4]==',' : continue
    time=line[0]
    tot_lines+=1
    
    cln=line[3]
    #print "line:",line
    #print "seq:",seq
    seq=line[4]
    noF=0
    
    F=seq.count("F")
    A=seq.count("A")
    key=time+'_'+cln
    if key in dc: 
      dc[key][0].append(F)
      dc[key][1].append(A)
    else:
      dc[key]=[[F],[A]]

    # abpr=A/(A+3.)*(math.exp(-1.*F))
    # growth=0.1*F/(F+10.)
    # l_growth.append(F)
    # l_abpr.append(A)


import alphashape
from descartes import PolygonPatch
from shapely.geometry import MultiPoint, Polygon

fig, ax = plt.subplots()
# my_cmap = plt.cm.get_cmap('tab20', len(dc))

dc={i:dc[i] for i in dc if len(dc[i][0])>1000}
# print(dc)
my_cmap = plt.cm.get_cmap('jet', len(dc))
print("Hello: there are len(dc) = ", len(dc), " hulls to draw")
for i,key in enumerate(dc):
  # print("Got: dc[key] =",dc[key])
  l_growth = dc[key][0]
  l_abpr = dc[key][1]
  print(i)
  
  points=[x for x in set(zip(l_growth,l_abpr))] #set to remove extra points... for concave hull they are not needed
  
  alpha = 0.5 * alphashape.optimizealpha(points)
  hull = alphashape.alphashape(points, alpha)
  # hull_pts = hull.exterior.coords.xy
  # print("Hello3")
  # fig, ax = plt.subplots()
  ax.scatter(l_growth,l_abpr,color=my_cmap(i), alpha=0.05)
  # ax.scatter(hull_pts[0], hull_pts[1], color=[my_cmap(i)], alpha=0.4)
  try:
    ax.add_patch( PolygonPatch( hull, fill=False, color=my_cmap(i), alpha=0.2 ) )

  except:
    print("These points give errors: ",points)
    pass

  gstep=1
  astep=1
  gedges = np.arange(0., max(l_growth)+gstep, step=gstep)
  aedges = np.arange(0., max(l_abpr)+astep, step=astep)
  H, xedges, yedges = np.histogram2d(l_growth, l_abpr, bins=(gedges, aedges))
  X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
  cs = ax.contour(X, Y, H.T, colors=[my_cmap(i)] )


  if i==0: break
ax.set_xlim([0,20])  
ax.set_ylim([0,120])  
plt.show()
sys.exit(1)

while False:

  gstep=1
  astep=1
  gedges = np.arange(0., max(l_growth)+gstep, step=gstep)
  aedges = np.arange(0., max(l_abpr)+astep, step=astep)
  # print("maxes: ",max(l_growth),max(l_abpr) )
  H, xedges, yedges = np.histogram2d(l_growth, l_abpr, bins=(gedges, aedges))
  X, Y = np.meshgrid(xedges[:-1], yedges[:-1])
  
  H = np.ma.masked_where(H<=0,H)
  # pcmesh=ax.pcolormesh(X, Y, H.T)#,cmap="gist_earth_r")
  # print (len(X))
  # print(len(Y))
  # print("Hello 1,",i," my_cmap(i)=",my_cmap(i))
  cs = plt.contour(X, Y, H.T, colors=[my_cmap(i)] ) 


  # print("Hello 2,",i)
  # plt.colorbar(pcmesh)
  # break
plt.show()
sys.exit(1)

# gstep=0.005
# astep=0.05
# gedges=np.arange(0., .1+gstep, step=gstep)
# aedges=np.arange(0., 1.+astep, step=astep)
#H, xedges, yedges = np.histogram2d(l_growth, l_abpr, bins=(gedges, aedges))
# print("Log10 transforming Histogram H")
# H = np.where(H> 0, np.log10(H), 0)

gstep=1
astep=1
gedges=np.arange(0., max(l_growth)+gstep, step=gstep)
aedges=np.arange(0., max(l_abpr)+astep, step=astep)
H, xedges, yedges = np.histogram2d(l_growth, l_abpr, bins=(gedges, aedges))


X, Y = np.meshgrid(xedges, yedges)

pcmesh=ax.pcolormesh(X[:-1], Y[:-1], H.T,cmap="gist_earth_r")
plt.colorbar(pcmesh)
plt.show()
sys.exit(1)

    # if abpr>0.1 and growth<0.01:
    # #if seq.count("F")<=2:
    #   # print( abpr)
    #   noF=1
    #   tot_linesnoF+=1
    #   #print "noF true", seq
    # key=time+'_'+cln
    # if key in dc: 
    #   dc[key][0]+=1
    #   dc[key][1]+=noF
    # else:
    #   dc[key]=[1,noF]

# l_mutfr = [ x[1]/float(x[0]) for key, x in dc.items() ]
# # step = 0.0005
# step = 0.002
# maxhist=0.15
# hist,bins= np.histogram(l_mutfr, bins=np.arange(0., maxhist, step=step))
# plt.plot(bins[:-1],hist,lw=2,label="distribution of AB mutants\nper colony")

# # hist2,bins2= np.histogram(l_abpr, bins=np.arange(0., 1., step=step))
# # plt.plot(bins2[10:-1],hist2[10:],lw=2,label="distribution of AB prod total")

# print ("tot_lines",tot_lines,"tot_linesnoF",tot_linesnoF)
# averagemutrate=tot_linesnoF/float(tot_lines)

# plt.plot([averagemutrate,averagemutrate],[0,0.667*max(hist)],lw=2,linestyle='dashed',label="avrg. fraction of\nAB mutants")
# plt.xlabel("fraction AB producing mutant per colony")
# plt.ylabel("colonies")
# plt.legend()
#print "No mut?"
#for key,x in dc.items():
#  if x[1]/float(x[0])==0: 
#    print key, x[0]

# plt.show()
