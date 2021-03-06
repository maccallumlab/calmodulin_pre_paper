##Plots KDEs for Trials 1-3##
#uses dat files generated by get_rmsds.py
#bandwidth of 0.1 which determines smoothness

import numpy as np 
import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity

trials = ['Trial1', 'Trial2', 'Trial3']
for i in trials:
    x = []
    with open('rmsds_'+i+'.dat', 'r') as fp:
        for line in fp:
            l = line.strip()
            x.append(float(l))
    x = np.array(x)
    kde = KernelDensity(bandwidth=0.1, kernel = 'gaussian').fit(x.reshape(-1,1))
    t = np.linspace(1, 9, 1000)
    kde_dens = np.exp(kde.score_samples(t[:, None]))
    plt.plot(t, kde_dens, '-', label = i, linewidth=1.5)

plt.legend(loc='upper right')
plt.xlim([0.75, 9.0])
plt.ylim([0.0, 2.25])
plt.ylabel('$density$')
plt.xlabel('$RMSD(\AA)$')
plt.savefig('kde_01_1-3_1ms.png')

