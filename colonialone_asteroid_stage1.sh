#!/bin/bash

# set output and error output filenames, %j will be replaced by Slurm with the jobid
#SBATCH -o stage1_hightol_%j.out
#SBATCH -e stage1_hightol_%j.err 

#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulumani@gwu.edu

# single node in the "short" partition
#SBATCH -N 1
#SBATCH -p short

# set the correct directory - cloned via git
#SBATCH -D /home/skulumani/asteroid

#SBATCH -J stage1
#SBATCH --export=NONE

# half hour timelimit
#SBATCH -t 2-00:00:00

module load matlab

# run asteroid shooting script
# matlab -nodesktop < asteroid_shooting.m
matlab -nodesktop < hpc_test_stage1.m


