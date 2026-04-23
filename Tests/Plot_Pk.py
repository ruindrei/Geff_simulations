## Function for plotting Pk ratios of self-interacting neutrinos / LCDM


## The idea is to input a LCDM Tk base file and then any number of other Tk files to plot.
## Format of command line:
## python Plot_Pk.py --base "lcdm.dat" --Pk "Pk_one.dat" "Pk_two.dat" ...

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import argparse

## Create the argument parser. 
parser = argparse.ArgumentParser(description='Plot power spectra from CLASS output files.')
parser.add_argument('--base', type=str, required=True, help='Path to the base CLASS output file containing the power spectrum.')
parser.add_argument('--Pk', type=str, nargs='+', required=True, help='Paths to the CLASS output files containing the power spectra to compare.')
args = parser.parse_args()

## Load the base power spectrum data.
base_data = np.loadtxt(args.base, skiprows=4)
k_base = base_data[:,0]
Pk_base = base_data[:,1]

## Now we can loop over the other Pk files and plot them.
if args.Pk is not None:
    for Pk_file in args.Pk:
        Pk_data = np.loadtxt(Pk_file, skiprows=4)
        k_geff = Pk_data[:,0]
        Pk_geff = Pk_data[:,1]
        

        ## We need to interpolate the base Pk to the same k values as the new file to compute the ratio.
        ## They are both loglog so we can interpolate in log space to be more accurate.
        interp_Pk_base = interp1d(k_base, Pk_base, kind='cubic', bounds_error=False, fill_value='extrapolate')
        Pk_base_interp = interp_Pk_base(k_geff)

        ## Now plot the ratio in loglog space.
        plt.plot(k_geff, Pk_geff / Pk_base_interp, label=Pk_file + ' / ' + args.base)
        plt.xlabel('k [h/Mpc]')
        plt.xscale('log')
        plt.ylabel('P(k) / P_LCDM(k)')
        plt.title('Power Spectrum Ratio')
        plt.ylim(0.7, 1.2)
        plt.legend()
        plt.show()
