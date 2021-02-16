#!/bin/bash
#SBATCH --gres=gpu:4
#SBATCH --nodes=12
#SBATCH --ntasks=48
#SBATCH --time=00-50:00
#SBATCH --tasks-per-node=8
#SBATCH -o CaM_Trial6_2.out
#SBATCH -e CaM_Trial6_2.err
#SBATCH --job-name="CaM6final2"

module load cuda

source activate meld-0.4.14

srun -l hostname
srun -l echo '$CUDA_VISIBLE_DEVICES'

mpirun -np 48 launch_remd --debug
