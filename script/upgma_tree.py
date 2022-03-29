#!/usr/bin/python3.6

import sys,math,os,string,random
from scipy.cluster.hierarchy import dendrogram, linkage
import scipy.spatial as scsp
from matplotlib import pyplot as plt
import numpy as np

# # From https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
# def fancy_dendrogram(*args, **kwargs):
#     max_d = kwargs.pop('max_d', None)
#     if max_d and 'color_threshold' not in kwargs:
#         kwargs['color_threshold'] = max_d
#     annotate_above = kwargs.pop('annotate_above', 0)

#     ddata = dendrogram(*args, **kwargs)

#     if not kwargs.get('no_plot', False):
#         plt.title('Hierarchical Clustering Dendrogram (truncated)')
#         plt.xlabel('sample index or (cluster size)')
#         plt.ylabel('distance')
#         for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
#             x = 0.5 * sum(i[1:3])
#             y = d[1]
#             if y > annotate_above:
#                 plt.plot(x, y, 'o', c=c)
#                 plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
#                              textcoords='offset points',
#                              va='top', ha='center')
#         if max_d:
#             plt.axhline(y=max_d, c='k')
#     return ddata

ab_list_pos_in_file = 2
genome_pos_in_file = 4
gene_ab_pos_in_file = 5

filename = sys.argv[1]
l_antib=[]
with open(filename,"r") as fin:
    for line in fin.readlines():
        line =line.split()
        if len(line)>4: 
            # print line
            if line[genome_pos_in_file] != 'n':
                if 'A' in line[genome_pos_in_file]:
                    # print line
                    l_antib.extend([x for x in line[gene_ab_pos_in_file].split(',')[:-1] ])
# print l_antib
# Make a set from l_antib, you need antib diversity
l_antib = list(set(l_antib))
l_antib = [format(int(x), '024b') for x in l_antib]
print( "So many different antibiotics: ", len(l_antib))
# print l_antib[1]

# sys.exit(1)
# # X = [[i] for i in [2, 8, 0, 4, 1, 9, 9, 0]]
# X = [[0,0,0,0,0,0,1],[0,1,1,1,0,1,1],[1,0,0,0,0,0,0],[1,0,0,0,1,1,1],[0,1,0,0,0,1,0],[1,1,0,1,0,1,1]]

# print X
# l_antib_forlinkage = random.sample(l_antib, len(l_antib)/100)
# l_antib_forlinkage = l_antib[0:100]

# random.shuffle(l_antib)
l_antib_forlinkage=l_antib[:]
l_antib_forlabel = l_antib_forlinkage[:]
# print "So many used for AB tree: ", len(l_antib_forlinkage)1

l_antib_forlinkage=[ [int(y) for y in x] for x in l_antib_forlinkage ]
# print l_antib_forlinkage[1]
# sys.exit(1)

# Calulate distance matrix
distmat = np.zeros((len(l_antib_forlinkage), len(l_antib_forlinkage)), dtype = int)
for i,X in enumerate(l_antib_forlinkage):
    for j,Y in enumerate(l_antib_forlinkage):
        distmat[i,j] = sum([0 if x==y else 1 for x,y in zip(X,Y)])
# print distmat[:3,:3]
# sys.exit(1)


# Z = linkage(X, 'ward')
print(sys.getrecursionlimit())
sys.setrecursionlimit(1500)
# sys.exit(1)
h, w = distmat.shape
Z = linkage(distmat[np.triu_indices(h, 1)], method='average', metric ='hamming')

###   ---   Z = linkage(l_antib_forlinkage, metric ='hamming') <- THIS IS A SUPER BUG IN SCIPY !!!!

# dn = dendrogram(Z,p=10,truncate_mode='lastp', ax=ax)#,labels=l_antib_forlabel)
#fig, (ax1,ax2) = plt.subplots(2)
fig, (ax1,ax2) = plt.subplots(2)
# dn = dendrogram(Z,p=90,truncate_mode='lastp', ax=ax1)#,labels=l_antib_forlabel)

# FOR DENDROGRAM USE THIS
# dn = dendrogram(Z, ax=ax1,color_threshold=2.1)

# For number of cluster finding
last = Z[-500:, 2]
last_rev = last[::-1]
idxs = np.arange(1, len(last) + 1)
ax2.plot(idxs, last_rev, label="data")

# Y=[-36,-25,-16,-9,-4,-1,0,1,2,3,4]
# X=[0,1,2,3,4,5,6,7,8,9,10]
# ax2.plot(X,Y,label="data")
# plt.show()
# sys.exit(1)
# smoothen
from scipy.signal import savgol_filter
yhat = savgol_filter(last_rev, 11, 3)
# yhat = savgol_filter(Y, 3, 2)
ax2.plot(idxs, yhat, label="smooth")
# ax2.plot(X, yhat, label="smooth")
# plt.show()
# sys.exit(1)
#2nd derivative
acceleration = np.diff(yhat, 2)  # 2nd derivative of the distances
acceleration_rev = acceleration #[::-1]
ax2.plot(idxs[:-2] + 1, acceleration_rev, label="2nd der")
# ax2.plot(np.array(X[:-2]) + 1, acceleration_rev, label="2nd der")
k = acceleration_rev.argmax() + 2  # if idx 0 is the max of this we want 2 clusters
print ("clusters:", k)
ax2.legend()

# Z = linkage(X, 'ward')

# fig = plt.figure(figsize=(25, 10))

# dn = dendrogram(Z,ax=ax2)

plt.show()
