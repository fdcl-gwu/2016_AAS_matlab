Running Asteroid Shooting README

##ColonialOne login
ssh sklumani@login.colonialone.gwu.edu

skulumani and GWU lastpass password

## Scheduler commands
sbatch name_of_script.sh
squeue -l 
sinfo

## Copy a file from colonialone to the local compute
pscp skulumani@login.colonialone.gwu.edu:/home/skulumani/asteroid/fine_asteroid_1.mat C:\Users\skulumani\Desktop\fine_asteroid_1.mat

## Send a file from SEAS computer to colonial one
pscp C:\Users\skulumani\Desktop\fine_asteroid_1.mat skulumani@login.colonialone.gwu.edu:/home/skulumani/asteroid/fine_asteroid_1.mat 

## Running the Code on local computer
run hpc_test.m

The number of angles is defined in asteroid_shooting.m - this defines the density of the reachability set (discretization)

## Run on HPC
Make sure you have the compiled version of polyhedron_potential_mex_1024
Modify file names in colonialone_asteroid.sh
Modify file name in asteroid_shooting.m for final mat file

## MEX Function
Need to compile polyhedron_potential.m using MEX
Use codegen to generate MEX function for polyhedron_potential.m
Test function is given by polyhedron_potential_test.m 
Modify the number of faces using 'false' - full number of faces or 'true' - 1024 faces

## Transfer examples

hpc_test.m/colonialone_asteroid.sh - 4 constraint transfer
hpc_test_4convert.m/colonialone_asteroid_4convert.sh - Extra vertical weighting 
hpc_test_z.m/colonialone_asteroid_z.sh - 5 constraint with the added z=0 constraint

asteroid_transfer_analysis.m - Plots the transfers using the saved mat files in ./results

## Final transfer

The results in ./results/5constraint_transfer shows the results used in the paper for AAS