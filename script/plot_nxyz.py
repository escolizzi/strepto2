#!/usr/bin/python3

"""
Simple plotting tool, takes a file and assumes that it is plotted so that first column is x axis
and all other columns are y axis
"""

import sys
import matplotlib as mpl
mpl.use('qt5agg')
import matplotlib.pyplot as plt


ldata=[]
with open(sys.argv[1],"r") as fin:
    for line in fin:
        line=line.split()
        ldata.append([ float(x) for x in line ])

ldata = [*zip(*ldata)]
for data in ldata[1:]:
    plt.plot(ldata[0],data)
plt.show()