#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from meld.remd import ladder, adaptor, master_runner
from meld import system
from meld import comm, vault
from meld import parse
from meld import util
from meld.system import patchers
import meld
import simtk.openmm as mm


N_REPLICAS = 48
N_STEPS = 50000
BLOCK_SIZE = 100

parse.allowed_residues.append('CA')

def gen_state(s, index):
    pos = s._coordinates
    pos = pos - np.mean(pos, axis=0)
    vel = np.zeros_like(pos)
    alpha = index / (N_REPLICAS - 1.0)
    energy = 0
    return system.SystemState(pos, vel, alpha, energy, [99999., 99999., 99999.])


def get_dist_restraints_pre(filename, s, scaler, ramp):
    groups = []
    for line in open(filename):
        cols = line.split()
        index_i, name_i = int(cols[0]), cols[1]
        index_j, name_j = int(cols[2]), cols[3]
        r_min, r_max = float(cols[4]), float(cols[5])
        k = float(cols[6])

        r_min_linear = max(0.0, r_min - 0.2)
        r_max_linear = r_max + 0.2

        rest = s.restraints.create_restraint('distance', scaler=scaler, ramp=ramp,
                                             atom_1_res_index=index_i, atom_1_name=name_i,
                                             atom_2_res_index=index_j, atom_2_name=name_j,
                                             r1=r_min_linear, r2=r_min,
                                             r3=r_max, r4=r_max_linear, k=k)
        group = s.restraints.create_restraint_group([rest], 1)
        groups.append(group)
    return groups


def get_dist_restraints(filename, s, scaler, ramp):
    groups = []
    for line in open(filename):
        cols = line.split()
        index_i, name_i = int(cols[0]), cols[1]
        index_j, name_j = int(cols[2]), cols[3]
        dist, k = float(cols[4]), float(cols[5])

        rest = s.restraints.create_restraint('distance', scaler=scaler, ramp=ramp,
                                             atom_1_res_index=index_i, atom_1_name=name_i,
                                             atom_2_res_index=index_j, atom_2_name=name_j,
                                             r1=0.0, r2=0.0, r3=dist, r4=dist+0.3, k=k)
        group = s.restraints.create_restraint_group([rest], 1)
        groups.append(group)
    return groups


def setup_system():
    # load the sequence
    protein_sequence = parse.get_sequence_from_AA1(filename='protein.dat')
    peptide_sequence = parse.get_sequence_from_AA1(filename='peptide.dat')
  
    n_res_protein = len(protein_sequence.split())
    n_res_peptide = len(peptide_sequence.split())

    # build the system
    protein = system.ProteinMoleculeFromSequence(protein_sequence)
    peptide = system.ProteinMoleculeFromSequence(peptide_sequence)
    protein.set_translation([100, 100, 150])
    peptide.set_translation([100, 150, 100])
    calcium1 = system.ProteinMoleculeFromSequence('CA') 
    calcium2 = system.ProteinMoleculeFromSequence('CA') 
    calcium3 = system.ProteinMoleculeFromSequence('CA') 
    calcium4 = system.ProteinMoleculeFromSequence('CA') 
    calcium1.set_translation([100, 105, 50])  
    calcium2.set_translation([100, 110, 50])
    calcium3.set_translation([100, 115, 50]) 
    calcium4.set_translation([100, 120, 50]) 

    rdc_patcher = patchers.RdcAlignmentPatcher(n_tensors=1)
    ond_patcher = patchers.VirtualSpinLabelPatcher(
           {17: 'OND',
            34: 'OND',
            42: 'OND',
            53: 'OND',
            86: 'OND',
            110: 'OND',
            117: 'OND',
            127: 'OND',
            143: 'OND',
            149: 'OND'})
 
    b = system.SystemBuilder()
    s = b.build_system_from_molecules([protein, calcium1, calcium2, calcium3, calcium4, peptide],
                                      leap_header_cmds="source leaprc.water.tip3p",
                                      patchers=[rdc_patcher, ond_patcher])

    s.temperature_scaler = system.GeometricTemperatureScaler(0, 0.3, 300., 550.)

    ramp = s.restraints.create_scaler('linear_ramp', start_time=1, end_time=200,
                                      start_weight=0, end_weight=1)

    #
    # Secondary Structure
    #
    ss_scaler = s.restraints.create_scaler('constant')
    protein_ss_rests = parse.get_secondary_structure_restraints(filename='protein_ss.dat',
                                                                system=s,
                                                                scaler=ss_scaler,
                                                                ramp=ramp,
                                                                torsion_force_constant=2.5,
                                                                distance_force_constant=2.5,
                                                                min_secondary_match=5)

    peptide_ss_rests = parse.get_secondary_structure_restraints(filename='peptide_ss.dat',
                                                                system=s,
                                                                scaler=ss_scaler,
                                                                ramp=ramp,
                                                                torsion_force_constant=2.5,
                                                                distance_force_constant=2.5,
                                                                first_residue=int(n_res_protein)+5)  # + 4 due to calciums
                                                                                                     # + 1 for 1-based indexing

    protein_ss_keep = int(len(protein_ss_rests) * 0.95) 
    peptide_ss_keep = int(len(peptide_ss_rests) * 0.95)
    s.restraints.add_selectively_active_collection(protein_ss_rests, protein_ss_keep)
    s.restraints.add_selectively_active_collection(peptide_ss_rests, peptide_ss_keep)


    #
    # Confinement Restraints
    #
    conf_scaler = s.restraints.create_scaler('constant')
    confinement_rests = []
    n_res = n_res_protein + n_res_peptide + 4
    for index in range(1, n_res + 1): 
        protein_rest = s.restraints.create_restraint('confine',
                                                     conf_scaler,
                                                     ramp=ramp,
                                                     res_index=index,
                                                     atom_name='CA',
                                                     radius=5,
                                                     force_const=250.0)
        confinement_rests.append(protein_rest)
      
    s.restraints.add_as_always_active_list(confinement_rests)

    #
    # Calcium restraints
    #
    scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.5, alpha_max=1.0, factor=4.0)
    calcium_rests = get_dist_restraints('calcium_restraints.dat', s, scaler, ramp)
    n_keep_calcium = len(calcium_rests)
    s.restraints.add_selectively_active_collection(calcium_rests, n_keep_calcium)

    #
    # PRE restraints
    #
    scaler_short = s.restraints.create_scaler('nonlinear', alpha_min=0.6, alpha_max=1.0, factor=4.0)
    scaler_medium = s.restraints.create_scaler('nonlinear', alpha_min=0.5, alpha_max=0.6, factor=4.0)
    scaler_long = s.restraints.create_scaler('nonlinear', alpha_min=0.4, alpha_max=0.5, factor=4.0)
    scalers = [scaler_short, scaler_medium, scaler_long]

    OND_list = [17, 34, 42, 53, 86, 110, 117, 127, 143, 149]
    for ond in OND_list:
        for length , i in zip(['short','medium','long'], range(3)):
            scaler = scalers[int(i)]
            pre_restraints = get_dist_restraints_pre('rest_files/'+str(ond)+'-pre-'+length+'.dat', s, scaler, ramp)
            n_keep_pre = int(len(pre_restraints) * 0.90)
            s.restraints.add_selectively_active_collection(pre_restraints, n_keep_pre)

    #
    # RDC Restraints
    # 
    rdc_scaler = s.restraints.create_scaler('nonlinear', alpha_min=0.3, alpha_max=0.4, factor=4.0, strength_at_alpha_max=1.0e-2)
    rdc_rests = parse.get_rdc_restraints(system=s,
                                         patcher=rdc_patcher, 
                                         scaler=rdc_scaler,
                                         ramp=ramp,
                                         quadratic_cut=1.0,
                                         scale_factor=1.0e4,
                                         filename='rdc.dat') 
    s.restraints.add_as_always_active_list(rdc_rests)


    # create the options
    options = system.RunOptions()
    options.implicit_solvent_model = 'obc'
    options.use_big_timestep = False
    options.cutoff = 1.8
    options.remove_com = True

    options.use_amap = True
    options.amap_beta_bias = 1.0
    options.timesteps = 25000
    options.minimize_steps = 5000

    # create a store
    store = vault.DataStore(s.n_atoms, N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner
    l = ladder.NearestNeighborLadder(n_trials=48 * 48)
    policy_1 = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy_1)

    remd_runner = master_runner.MasterReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS,
                                                            ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS)
    store.save_communicator(c)

    # create and save the initial states
    states = [gen_state(s, i) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()
    return s.n_atoms

setup_system()
