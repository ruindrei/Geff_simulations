## Code that computes the correct rescaling for A_s to mimic a massive neutrino cosmology.

## First let's load in our MPks for both cosmologies.
## THe idea is that the amplitudes of Pk need to match at k=1 h/Mpc and z=99 which is what our output redshift already is. 

import numpy as np
import matplotlib
#import pandas as pd ## Probably pandas will be most useful for .txt files
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d


A_s = 1.987e-9

## Both simulations have the same h. 

massless = np.loadtxt("Geff_find_rescaling_tests/lcdm_base_class_00_pk.dat", skiprows=4)
k_massless = massless[:,0]
Pk_massless = massless[:,1]

massive = np.loadtxt("Geff_find_rescaling_tests/lcdm_nu_base_class_00_pk.dat", skiprows=4)
k_massive = massive[:,0]
Pk_massive = massive[:,1]

## Plot 
plt.loglog(k_massless, Pk_massless, label = "Massless Neutrinos")
plt.loglog(k_massive, Pk_massive, label = "Massive Neutrinos")
plt.xlabel("k [h/Mpc]")
plt.ylabel("P(k) [(Mpc/h)^3]")
plt.legend()
plt.show()

## Now we need to rescale A_s to match the amplitudes at k=1 h/Mpc.
## We can do this by interpolating the Pk values at k=1 h/Mpc
k_target = 1.0  # h/Mpc
# Interpolate Pk for massless neutrinos at k_target
log_k_massless = np.log10(k_massless)
log_Pk_massless = np.log10(Pk_massless)
interp_massless = interp1d(log_k_massless, log_Pk_massless, kind='cubic', bounds_error=False, fill_value='extrapolate')
log_Pk_massless_target = interp_massless(np.log10(k_target))
Pk_massless_target = 10**log_Pk_massless_target
# Interpolate Pk for massive neutrinos at k_target
log_k_massive = np.log10(k_massive)
log_Pk_massive = np.log10(Pk_massive)
interp_massive = interp1d(log_k_massive, log_Pk_massive, kind='cubic', bounds_error=False, fill_value='extrapolate')
log_Pk_massive_target = interp_massive(np.log10(k_target))
Pk_massive_target = 10**log_Pk_massive_target

## Now compute the rescaling factor for A_s
rescaling_factor = Pk_massive_target / Pk_massless_target
print(f"Rescaling factor for A_s to match P(k) at k={k_target} h/Mpc: {rescaling_factor:.4f}")

## The new A_s for the massless neutrino cosmology to match the amplitude at k=1 h/Mpc is:
A_s_rescaled = A_s * rescaling_factor
print(f"Rescaled A_s: {A_s_rescaled:.4e}")

### Now we shall run the rescaled CLASS cosmology and see that it worked. 
rescaled = np.loadtxt("Geff_find_rescaling_tests/lcdm_rescaled_base_class_00_pk.dat", skiprows=4)
k_rescaled = rescaled[:,0]
Pk_rescaled = rescaled[:,1]

## Plot all three together to check
plt.loglog(k_massless, Pk_massless, label = "Massless Neutrinos")
plt.loglog(k_massive, Pk_massive, label = "Massive Neutrinos")
plt.loglog(k_rescaled, Pk_rescaled, label = "Rescaled Massless Neutrinos", linestyle='--')
plt.xlabel("k [h/Mpc]")
plt.ylabel("P(k) [(Mpc/h)^3]")
plt.legend()
plt.show()

## Check by interpolation as well:
log_k_rescaled = np.log10(k_rescaled)
log_Pk_rescaled = np.log10(Pk_rescaled)
interp_rescaled = interp1d(log_k_rescaled, log_Pk_rescaled, kind='cubic', bounds_error=False, fill_value='extrapolate')
log_Pk_rescaled_target = interp_rescaled(np.log10(k_target))
Pk_rescaled_target = 10**log_Pk_rescaled_target
print(f"Pk at k={k_target} h/Mpc for rescaled cosmology: {Pk_rescaled_target:.4e}")
print(f"Pk at k={k_target} h/Mpc for massive neutrino cosmology: {Pk_massive_target:.4e}")


## PLot the ratio of rescaled to massive case
## First interpolate both to a common k array (we can just use the rescaled k array since it should be the same as the others)
log_k_common = np.log10(k_rescaled)
log_Pk_massive_interp = interp_massive(log_k_common)
log_Pk_rescaled_interp = interp_rescaled(log_k_common)
Pk_massive_interp = 10**log_Pk_massive_interp
Pk_rescaled_interp = 10**log_Pk_rescaled_interp
plt.plot(k_rescaled, Pk_rescaled_interp / Pk_massive_interp, label="Rescaled / Massive", color='purple')
plt.xscale('log')
plt.xlabel("k [h/Mpc]")
plt.ylabel("P(k) Ratio")
plt.axhline(1, color='gray', linestyle=':')
plt.legend()
plt.show()
