This directory tests if we can recreate Adam He's Figure 1 from his paper: https://arxiv.org/abs/2503.15592v2.
The .ini files containing Geff must be run with a modified Github class code: 
https://github.com/ash2223/class-interacting-neutrinos-PT-lya

template_geff.ini:
Has the best-fit parameters of Table IV from Adam's paper. 
This has 1 massive neutrino (m_nu = 0.303) and 2 massless ones.
Log10_Geff_nu = -4.477
N_ur is scaled down to match 2 neutrino species.
NOTE: Omega_Lambda may be scaled down from the best-fit value in order to ensure flatness. 

template_lcdm_planck.ini:
Has the parameters of the 2nd to last column of Table 2 from https://arxiv.org/pdf/1807.06209.
This has 1 massive neutrino (m_nu = 0.06) and 2 massless ones. 
Log10_Geff_nu = 0
NOTE: Omega_Lambda may also be scaled down from the best-fit value in order to ensure flatness.

Once the files are run with the modified CLASS code, they should output pk and tk.dat files that are included. 
The z4_pk.dat files are at redshift 2.0.

Plot_Pk.py is a python code that can be used for plotting a Pk ratio from terminal
It takes two arguments, --base and --Pk
--base is the string filepath to the lcdm.dat file that we want to divide by.
--Pk is the string filepath to the numerator geff.dat file that contains our self-interacting neutrinos. 

Example:
python Plot_pk.py --base "lcdm_planck_z4_pk.dat" --Pk "lcdm_geff_z4_pk.dat"

This will output and save a plot that shows the Pk ratio of the Pk shown in lcdm_geff_z4_pk.dat divided by lcdm_planck_z4_pk.dat.

The plot for these example files is shown in Pk_ratio_geff_to_planck.png. It matches the correct shape from Figure 1. 
