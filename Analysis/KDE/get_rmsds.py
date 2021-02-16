##Calculates RMSDs to reference crystal structure of backbone (omitting flexible tails)
#generates rmsd files for plot_kde.py

import numpy as np 
import mdtraj as md
import matplotlib.pyplot as plt
from sklearn.neighbors import KernelDensity


native = md.load('cam_fill.pdb')
native_indices = native.top.select('backbone and (resid 4 to 146 or resid >=149)')

traj = md.load_dcd('CaM_Trial3.dcd', top='trajectory_1.pdb')[-10000:]
traj_indices = traj.top.select('backbone and (resid 4 to 146 or resid >=149)')

print(len(native_indices), len(traj_indices))

ref = native.atom_slice(atom_indices=native_indices)
cam = traj.atom_slice(atom_indices = traj_indices)
print(ref, cam)


rmsds_to_native = md.rmsd(cam, ref)*10

with open('rmsds_Trial3.dat', 'w') as fp:
    for i in rmsds_to_native:
        fp.write(str(i)+'\n')

