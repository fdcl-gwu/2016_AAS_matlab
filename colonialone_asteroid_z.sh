#!/bin/bash

# set output and error output filenames, %j will be replaced by Slurm with the jobid
#SBATCH -o stage3_z_%j.out
#SBATCH -e stage3_z_%j.err 

#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulumani@gwu.edu

# single node in the "short" partition
#SBATCH -N 1
#SBATCH -p short

# set the correct directory - cloned via git
#SBATCH -D /home/skulumani/5con_all

#SBATCH -J ast_z_2
#SBATCH --export=NONE

# half hour timelimit
#SBATCH -t 2-00:00:00

module load matlab

# run asteroid shooting script
# matlab -nodesktop < asteroid_shooting.m
matlab -nodesktop < hpc_test_z.m


