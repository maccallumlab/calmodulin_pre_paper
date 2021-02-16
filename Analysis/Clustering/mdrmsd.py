##Clustering script for CaM_Trials##
#clusters using HDBSCAN the last 1 microsecond of simulation
#uses rmsd to native of backbone (excluding flexible tails but including peptide) as distance metric

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import hdbscan

MIN_SAMPLES = 200  #determined from data/trial and error

#Calculate rmsd to native
def compute_rmsd_matrix(traj):
    distances = np.empty((traj.n_frames, traj.n_frames))

    for i in range(traj.n_frames):
        distances[i] = md.rmsd(traj, traj, i, atom_indices=traj.top.select('backbone'))
    return distances

#Determines the k-plot (helpful in determining MIN_SAMPLES)
def plot_k_dist(distances):
    print('plotting dists')
    s = np.sort(distances, axis=0)
    counts = s[:, MIN_SAMPLES]

    plt.plot(counts)
    plt.xlabel('distance')
    plt.ylabel('num_steps')
    plt.savefig('k-dist.png')

#Clusters data using HDBSCAN
def make_clusters(native, traj):
    distances = compute_rmsd_matrix(traj)
    plot_k_dist(distances)

    #clustering set up
    clusterer = hdbscan.HDBSCAN(min_cluster_size=MIN_SAMPLES)
    cluster_indices = clusterer.fit_predict(distances)

    min_index = 0
    max_index = np.max(cluster_indices) + 1

    #clustering
    clusters = [traj[np.where(cluster_indices == index)]
        for index in range(min_index, max_index)]
    clusters = sorted(clusters, key=lambda x: x.n_frames, reverse=True)
    
    #now add the unclustered frames to last cluster
    clusters.append(traj[np.where(cluster_indices == -1)])

    cluster_sizes = [c.n_frames for c in clusters]
    total_frames = traj.n_frames

    print('Found {} total clusters.'.format(len(clusters)))

    #calculates important values and outputs to files
    for i, c in enumerate(clusters):
        rmsds_to_native = md.rmsd(c, native)*10
        mean = np.mean(rmsds_to_native)
        median = np.median(rmsds_to_native)
        min_ = np.min(rmsds_to_native)
        max_ = np.max(rmsds_to_native)
        std_ = np.std(rmsds_to_native)
        np.savetxt("clusters_0"+str(i)+".dat", rmsds_to_native, fmt="%f")
        print('Cluster {:02d} has population {:.1f}; RMSD: {:.2f} {:.2f} {:.2f} {:.2f} {:.2f}'.format(i, 100 * cluster_sizes[i] / float(total_frames), mean, median, min_, max_, std_))
        c.save('cluster_{:02d}.pdb'.format(i))

#native struct
native = md.load('cam_fill.pdb')
native_indices = native.top.select('backbone and (resid 4 to 146 or resid>=149)')

#last 1 microsecond of simulation
traj = md.load_dcd('CaM_Trial3.dcd', top='trajectory_1.pdb')[-10000:]
traj_indices = traj.top.select('backbone and (resid 4 to 146 or resid >=149)')

#gets indices of subsection
ref = native.atom_slice(atom_indices=native_indices)
cam = traj.atom_slice(atom_indices = traj_indices)

make_clusters(ref, cam)
