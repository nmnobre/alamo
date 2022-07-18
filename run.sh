#!/bin/bash
#SBATCH -p compute
#SBATCH -J disconnections             # job name
#SBATCH -o /mmfs1/scratch/output_%j_stdout          # output file name (%j expands)
#SBATCH -N 1                 # total number of nodes 
#SBATCH -n 128
#SBATCH -t 24:00:00          # run time (hh:mm:ss)

module swap openmpi4 mvapich2/2.3.4

mpirun ~/alamo/bin/alamo-2d-debug-g++ input plot_file=/mmfs1/scratch/output_${SLURM_JOB_ID} 
