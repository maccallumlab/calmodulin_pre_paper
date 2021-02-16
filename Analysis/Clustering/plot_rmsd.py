##Generates histograms based on clustering results of mdrmsd.py

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle
clusts=[]
yz = []
binz =[]
for i in range (0, 3):
    x = []
    with open("clusters_0"+str(i)+".dat", 'r') as fp:
    	for line in fp:
            line = line.split()
            x.append(float(line[0]))
    y,bins, _ = plt.hist(x, bins=200, label='cluster '+str(i))
    clusts.append(x)
    yz.append(y)
    binz.append(bins)

a = np.matrix(clusts[0])
plt.ylabel('Count')
plt.xlabel('RMSD')
plt.legend(loc='best')
plt.xlim(1.5,8)
plt.savefig('clusters_200h.png')

print("mean: ", np.mean(a))
print("min: ", np.min(a))
bin_max = np.where(yz[0]==yz[0].max())
print("mode: val=", binz[0][bin_max][0], " size = ",yz[0].max())
