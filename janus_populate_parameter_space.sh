#!/bin/bash

#SBATCH --job-name chaffin-populate-parameter-space
#SBATCH --time 240:00
#SBATCH --nodes 64
#SBATCH --ntasks 1024
#SBATCH --output janus_populate_parameter_space.out

cd /projects/nick/mike/3d_thomas_model/

make simulate_coronal_scan

lb paramater_space_midpoint_simulation_commands
