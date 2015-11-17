#!/bin/bash

#SBATCH --job-name chaffin-populate-parameter-space
#SBATCH --time 10:00
#SBATCH --nodes 1
#SBATCH --ntasks 1
#SBATCH --output janus_test_run.out

cd /projects/nick/mike/3d_thomas_model/

make simulate_coronal_scan

lb test_simulation_command
