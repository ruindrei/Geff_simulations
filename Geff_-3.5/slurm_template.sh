#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name= 
#SBATCH --time= 
#SBATCH --nodes= 
#SBATCH --ntasks-per-node=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbird@ucr.edu
#SBATCH -A AST21005
export OMP_NUM_THREADS=28
ibrun MP-GenIC _genic_params.ini
ibrun MP-Gadget mpgadget.param
