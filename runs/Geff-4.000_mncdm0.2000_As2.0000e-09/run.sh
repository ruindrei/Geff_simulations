#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=Geff-4.000_mncdm0.2000_As2.0000e-09
#SBATCH --time= 
#SBATCH --nodes= 
#SBATCH --ntasks-per-node=2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sbird@ucr.edu
#SBATCH -A AST21005
export OMP_NUM_THREADS=28
ibrun MP-GenIC paramfile.genic
ibrun MP-Gadget paramfile.gadget
