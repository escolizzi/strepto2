#!/usr/bin/python3

# import sys,math
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from PIL import Image

sns.diverging_palette(220, 20, as_cmap=True)

lF = [ i for i in range(10)]
lA = [ i for i in range(10)]

grate = np.empty([len(lF),len(lA)], dtype=np.float32)
abrate = np.empty([len(lF),len(lA)], dtype=np.float32)
color = np.empty([len(lF),len(lA),3], dtype=int)
# color = []

for i,nF in enumerate(lF):
  for j,nA in enumerate(lA):
    grate = nF/(10. + nF)
    abrate = nA/(3. + nA) * np.exp(-nF)
    color[i,j] = [int(255*grate),int(255*abrate),0]
  # color.append([])

# im = Image.fromarray(color)
# im.show()

print(color)
plt.xticks(range(0,50,5))
plt.yticks(range(0,50,5))
plt.xlabel('nA')
plt.ylabel('nF')
plt.plot_surface(lF, lA, color)
# plt.imshow(color,origin='lower')
plt.show()