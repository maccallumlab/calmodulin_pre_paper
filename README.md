# README

This is a brief overview of the simulations and the analysis. Note that there is some redundency.

## Simulation

**setup.py**

Script to run MELD simulation for calmodulin (Trial3). All necesary restraint files are included.

See github.com/maccallumlab/meld for more information on running MELD.

**job.sh**

Submits MELD simulation to cluster. Requires simultaneous access to 48 GPUs with CUDA. 

See github.com/maccallumlab/meld for more information on running MELD.

### rest_files

Contains restraint files for PRE-derived restraints (split by location and distance) 


## Analysis

### Clustering

**mdrmsd.py**

Clustering protocol. HDBSCAN with min_cluster_size=200

Input | Description
----- | -----------
cam_fill.pdb	    |	PDB of native
CaM_Trial{}.dcd	    |	DCD of trajectory generated using meld "extract_trajectory extract_traj_dcd Cam_Trial{}.dcd --start=0 --end=50000"
trajectory_1.pdb    |	PDB of first frame of trajectory.pdb generated manually (functions as topology file)

Output | Description
------ | -----------
k-dist.png	    |	Visualized k-dist to help determine min_cluster_size
clusters_0{}.dat    |	RMSDs of clusters
cluster_0{}.pdb	    |	PDBs of clusters



**plot_rmsd.py**

Generates histograms of clusters.

Input | Description
----- | -----------
clusters_0{}.dat |	RMSDs of clusters (from mdrmsd.py)

Output | Description
------ | -----------
clusters_200h.png |	Histogram image


### KDE

**get_rmsds.py**

Determines rmsds to native for plotting KDEs (partially redundant, easier to use when only plotting KDEs than mdrmsd.py).

Excludes flexible tails, includes peptide.

Input | Description
----- | -----------
cam_fill.pdb	    |	PDB of native
CaM_Trial{}.dcd	    |	DCD of trajectory generated using meld "extract_trajectory extract_traj_dcd Cam_Trial{}.dcd --start=0 --end=50000"
trajectory_1.pdb    |	PDB of first frame of trajectory.pdb generated manually (used as topology file)

Output | Description
------ | -----------
rmsds_Trial{}.dat |	RMSDs of native to each frame in trajectory (1 microsecond)



**plot_kde.py**

Plots KDEs using .dat files from get_rmsds.py. Bandwidth of 0.1 for smoothness.

Input | Description
----- | -----------
rmsds_Trial{}.dat |	RMSDs of native to each frame in trajectory (1 microsecond) (from get_rmsds.py)

Output | Description
------ | -----------
kde_01_{}.png	|	KDE plot
