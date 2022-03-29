#!/usr/bin/python3
import matplotlib as mpl
mpl.use('qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
import numpy as np

lx = np.linspace(0,15,150)
l_growth = np.array([0.1*x/(10.+x) for x in lx])
l_ab = np.array([np.exp(-x) for x in lx])

fig, ax1 = plt.subplots()

color_growthrate= 'tab:blue'
color_ab= 'tab:orange'

ax1.plot(lx,l_growth,color=color_growthrate, lw=2)
ax1.set_xlabel('nr. growth genes', fontweight='bold')

ax1.set_ylabel('Growth rate', color=color_growthrate, fontweight='bold')
ax1.tick_params(axis='y', labelcolor=color_growthrate)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
ax2.plot(lx,l_ab,color=color_ab, lw=2)
ax2.set_ylabel('Max AB production', color=color_ab, fontweight='bold')  # we already handled the x-label with ax1
ax2.tick_params(axis='y', labelcolor=color_ab)


fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()