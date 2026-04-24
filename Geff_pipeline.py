## Pipeline for generating .ini files, running CLASS to generate the Pk files, computing MP-Gadget inputs and filling out templates.
"""
Script inputs:
Can be a text file? Or from terminal? 

--class          ## str filepath to desired CLASS code.
--A_s            ## A_s value will be a some kind of float value "2.109e-9"
--n_s            ## n_s spectral tilt value
--h              ## h cosmological value
--m_ncdm         ## neutrino mass in eV (if not rescaling)
--log10_G_eff_nu ## Geff self-interacting neutrino strength. 
--z_reio	 ## For later usage in creating the HI patchy reionization files. 


---> Create .ini file with default unchanging and the passed parameters.
---> Run CLASS code for .ini file and output pk file and make sure you know which file has been outputted.
---> Save plots of Pk files. IDK if we want to renormalize them to LCDM+massnu or what. 

We now need inputs for the MP-Gadget paramfiles for paramfile.gadget and paramfile.genic.
These inputs are derived from the initial inputs and should require no additional input from the user.

paramfile.genic:

--Omega0	# The total matter density (needs to be calculated from massive nu + cdm + baryons)
--OmegaLambda	# Omega_Lambda needs to be calculated to ensure flatness
--OmegaBaryon	# Should be constant.
--Omega_ur	# May need to calculate this if we are using massive neutrinos (no rescaling). Ask Simeon about this equation.
--HubbleParam	# Derived from h of course.

--FileWithInputSpectrum	# Taken from the _pk.dat output file from earlier
--FileWithTransferFunction	# Taken from the _tk.dat output file from earlier

--PrimordialIndex	# n_s value from earlier
--PrimordialAmp		# A_s value from earlier

--Seed		# Random number for seed, IDK if we want to change this. 


paramfile.gadget:

--OutputList	# List of output scale factors (I have a script to compute this)
--UVFluctuationFile	# NEED TO GENERATE THIS FILE FROM THE paramfile.genic FILE. Once this is output: write the name of the file. 

slurm.sh file:

--job_name	# What we want to call the job. This can be a mix of all the parameters and values in this case Geff (since that is what we are changing)

VALIDATION:

Use MP-Gadget/tools to run the genic file and see what .ini file is generated. Compare to the original .ini file to make sure all parameters are consistent. 
"""



## First let's start writing the .ini file parameters that DO NOT CHANGE.
import numpy as np
import os
import subprocess
import argparse
from textwrap import dedent

def _base_ini_content(A_s, n_s, h, m_ncdm, log10_G_eff_nu=None):
    content = dedent(f"""
        # CLASS .ini file for Geff pipeline

        # Precision Parameters (DO NOT CHANGE)

        k_per_decade_for_pk = 50
        k_bao_width = 8
        k_per_decade_for_bao = 200
        neglect_CMB_sources_below_visibility = 1e-30
        transfer_neglect_late_source = 3000.0
        l_max_g = 50
        l_max_ur = 150
        gauge = synchronous
        non_linear = none # default

        # Cosmological parameters
        A_s = {A_s}
        n_s = {n_s}
        h = {h}
        m_ncdm = {m_ncdm}
    """)

    if log10_G_eff_nu is not None:
        content += f"log10_G_eff_nu = {log10_G_eff_nu}\n\n"

    root_value = f"{'geff_' + str(log10_G_eff_nu) + '_' if log10_G_eff_nu is not None else 'lcdm_'}"

    content += dedent(f"""
        #
        omega_cdm = 0.1215
        omega_b = 0.02218
        Omega_k = 0
        Omega_fld = 0.0

        alpha_s = 0.0
        T_cmb = 2.7255
        N_ur = 2.030796442906864
        N_ncdm = 1.

        T_ncdm = 0.71611
        tau_reio = 0.0534

        ## Output parameters
        output = dTk vTk mPk
        root = {root_value}
        P_k_max_h/Mpc = 214.46605848506323
        z_max_pk = 100.0
        z_pk = 99.0, 4.0, 3.0, 2.0
        extra metric transfer functions = y
    """)
    return content


def write_ini_file(A_s, n_s, h, m_ncdm, log10_G_eff_nu):
    print("Writing .ini file with given parameters...")
    ini_content = _base_ini_content(A_s, n_s, h, m_ncdm, log10_G_eff_nu=log10_G_eff_nu)
    ini_filename = f"class_geff_{log10_G_eff_nu}.ini"
    with open(ini_filename, "w") as ini_file:
        ini_file.write(ini_content)
    print(f".ini file '{ini_filename}' has been written with the specified parameters.")

    print("Writing LCDM .ini file with given parameters...")
    ini_content_lcdm = _base_ini_content(A_s, n_s, h, m_ncdm, log10_G_eff_nu=None)
    ini_filename_lcdm = "class_lcdm.ini"
    with open(ini_filename_lcdm, "w") as ini_file:
        ini_file.write(ini_content_lcdm)
    print(f"LCDM .ini file '{ini_filename_lcdm}' has been written with the specified parameters.")

    return ini_filename, ini_filename_lcdm


## Now make a function to run the CLASS code with the generated .ini file and output the pk file.
def run_class(class_path, ini_filename,log10_G_eff_nu):
    print("Running CLASS from {class_path} with the generated .ini file...")
    # Construct the command to run CLASS
    ## The command will look something like /path/to/class/class /path/to/ini/file.ini

    command = f"{class_path} {ini_filename}"
    # Run the command using subprocess
    try:
        result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print("CLASS ran successfully.")
        print("Output:", result.stdout.decode())
    except subprocess.CalledProcessError as e:
        print("Error running CLASS:")
        print(e.stderr.decode())
        raise e
    
    ## This should have run CLASS and generated our output .dat files. 
    ## Return the expected pk files based on the .ini filename for later usage.
    ## The way CLASS does it is: lcdm_geff_+str(log10_G_eff_nu) + "_" + "z1_pk.dat" for the pk files where z() changes based on the redshift. For this purpose we want z=99 so z1 is okay. 
    pk_filename = f"geff_{log10_G_eff_nu}_00_z1_pk.dat"
    tk_filename = f"geff_{log10_G_eff_nu}_00_z1_tk.dat"
    return pk_filename, tk_filename ## Saved for later in MP-Gadget. 


## Now we need to make plots of the Pk files compared to their LCDM couunterparts. 

def plot_pk_files_matching(geff_pk_file, lcdm_pk_file, log10_G_eff_nu, scale=True):
    ## Can be scaled or base lcdm depending. 
    ## scale = True means you are plotting Pk with geff over Pk without geff.
    ## scale = False means you plot Pk with geff over planck 2018 data. 
    import matplotlib.pyplot as plt

    # Load the Pk data from the files
    geff_data = np.loadtxt(geff_pk_file, skiprows = 4)
    lcdm_data = np.loadtxt(lcdm_pk_file, skiprows = 4)

    k_geff = geff_data[:, 0]
    pk_geff = geff_data[:, 1]

    k_lcdm = lcdm_data[:, 0]
    pk_lcdm = lcdm_data[:, 1]

    ## Interpolate the k_lcdm and pk_lcdm to the k_geff values for a more direct comparison.
    from scipy.interpolate import interp1d
    interp_lcdm = interp1d(k_lcdm, pk_lcdm, kind='cubic', bounds_error=False, fill_value="extrapolate")
    pk_lcdm_interp = interp_lcdm(k_geff)


    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(k_geff, pk_geff/pk_lcdm_interp, label='Geff Pk', color='blue') ## Plot the ratio
    plt.xscale('log')
    plt.xlabel('k [h/Mpc]')
    plt.ylabel('P(k) / P_LCDM(k)')
    plt.title('')
    plt.legend()
    if scale == True:
        plt.savefig(f"pk_comparison_{log10_G_eff_nu}.png")
        print(f"Pk comparison plot saved as 'pk_comparison_{log10_G_eff_nu}.png'.")

    else:
        plt.savefig(f"pk_comparison_{log10_G_eff_nu}_planck.png")
        print(f"Pk comparison plot saved as 'pk_comparison_{log10_G_eff_nu}_planck.png'.")
    plt.show()

### Now we need functions to create the paramfile.genic file for MP-Gadget. 
def create_paramfile_gadget(pk_filename, tk_filename, A_s, n_s, h, m_ncdm, log10_G_eff_nu):
    print("Creating paramfile.genic for MP-Gadget...")

    ## Calculate the derived parameters for the genic file. 
    T_CMB = 2.7255 ## K
    Omega_fld = 0.0
    Omega_g = 4.480075654158969e-07 * T_CMB**4 / h**2
    Omega_b = 0.02218 / (h**2)
    Omega_cdm = 0.1215 / (h**2)
    Omega_0  = Omega_b + Omega_cdm + m_ncdm / (93.14 * h**2) ## Assuming one massive neutrino species for simplicity.

    """ISSUE: BE CAREFUL WITH THE CALCULATION FOR 1 NEUTRINO ASK SIMEON!!!!!!!!!!!"""
    Omega_ur = Omega_g * 0.22710731766023898 * (3.046 - 1.013198221453432 * 2) ## For 3 massive neutrinos? Fix it for 1 massive neutrino. 
    N_ur = Omega_ur / Omega_g / 0.22710731766023898 + 1.013198221453432* (3 - N_ncdm) 

    print(f"Calculated Omega_ur: {Omega_ur}, N_ur: {N_ur}")
    Omega_Lambda = 1 - Omega_fld - Omega_0 - Omega_g - Omega_ur
    MNue = m_ncdm ## Assuming one massive neutrino species for simplicity.

    paramfile_genic_content = dedent(f"""
        OutputDir = output # Directory for outputi
        FileBase = IC              # Base-filename of output files

        Ngrid = 128 # Size of cubic grid on which to create particles.

        BoxSize = 15000   # Periodic box size of simulation

        Omega0 = {Omega_0}      # Total matter density  (at z=0)
        OmegaLambda = {Omega_Lambda}      # Cosmological constant (at z=0)
        OmegaBaryon = {Omega_b}     # Baryon density        (at z=0)
        Omega_ur = {Omega_ur} # You do need to set this to change N_ur for the CLASS output. CHECK THIS NUMBER!!!

        ProduceGas = 1         # 1 = Produce gas  0 = no gas, just DM.
        RadiationOn = 1         #Enable background radiation
        HubbleParam = {h}      # Hubble parameter (may be used for power spec parameterization)

        Redshift = 99        # Starting redshift

        DifferentTransferFunctions = 1

        # filename of tabulated input spectrum in Mpc/h units
        FileWithInputSpectrum = {pk_filename}
        FileWithTransferFunction = {tk_filename}

        PrimordialIndex = {n_s}
        PrimordialAmp = {A_s}
        Seed = 181170    #  seed for IC-generator


        # Massive neutrinos
        MNue = {m_ncdm} # eV, mass of neutrino species.
 
        UnitLength_in_cm = 3.085678e21   # defines length unit of output (in cm/h) 
        UnitMass_in_g = 1.989e43      # defines mass unit of output (in g/cm)
        UnitVelocity_in_cm_per_s = 1e5 # defines velocity unit of output (in cm/sec)
    """)

    ## Save the paramfile.genic content to a file.
    genic_filename = "paramfile.genic"
    with open(genic_filename, "w") as f:
        f.write(paramfile_genic_content)
    print(f"paramfile.genic has been created and saved as '{genic_filename}' with the specified parameters.")

    ## Now write the paramfile.gadget file with the appropriate parameters.

    print("Creating paramfile.gadget for MP-Gadget...")
    paramfile_gadget_content = dedent(f"""
        InitCondFile = output/IC
        OutputDir = output
        TreeCoolFile = ../TREECOOL_fg20_eff_default.dat
        OutputList = "0.16666666666666666,0.18181818181818182,0.2222222222222222,0.25,0.3333333333333333"

        # CPU time -limit

        TimeLimitCPU = 43200 #= 12 hrs

        #  Characteristics of run
        TimeMax = 0.33333

        Omega0 = {Omega_0}   # Total matter density  (at z=0)
        OmegaLambda = {Omega_Lambda}    # Cosmological constant (at z=0)
        OmegaBaryon = {Omega_b}     # Baryon density        (at z=0)
        Omega_ur = {Omega_ur} # You do need to set this to change N_ur for the CLASS output. CHECK THIS NUMBER!!!

        HubbleParam = {h}      # Hubble paramater (may be used for power spec parameterization)

        CoolingOn = 1
        StarformationOn = 1
        RadiationOn = 1
        HydroOn = 1
        DensityIndependentSphOn = 1 # Should be on?
        MetalReturnOn = 0 ## Set this to 0 for PRIYA.
        BlackHoleOn = 1 ## also on. Defaults okay
        WindOn = 1 ## Should be on! Defaults okay
        StarformationCriterion = density #,h2?
        MassiveNuLinRespOn = 1 # For neutrinos.
        UVFluctuationFile=HI_table_massive
        HIReionTemp=15000 # Can go to 15000 as well

        MetalCoolFile = ../cooling_metal_UVB

        SnapshotWithFOF = 1 # Don't know if necessary
        FOFHaloLinkingLength = 0.2
        FOFHaloMinLength = 32

        MinGasTemp = 5.0


        #----------------Neutrinos---------------------
        MNue = {m_ncdm}


        #  Further parameters of SPH
        #  #Only kernel supported by fake_spectra
        DensityKernelType = cubic
        InitGasTemp = 270.
        MinGasTemp = 100

        # Remove He bubbles:
        QSOLightupOn=0 # Turn off if box size < 60 cMpc/h

        #----------------------BH Stuff-------------------------
        BlackHoleKineticOn = 1 # switch to kinetic feedback mode when the BH accretion rate is low

        BlackHoleFeedbackFactor = 0.05
        BlackHoleFeedbackMethod = spline | mass
        SeedBlackHoleMass = 5.0e-5
        BlackHoleAccretionFactor = 100.0
        BlackHoleNgbFactor = 2.0
        BlackHoleEddingtonFactor = 3.0

        MinFoFMassForNewSeed = 1
        TimeBetweenSeedingSearch = 1.03
        WriteBlackHoleDetails = 1

        #---------------------Wind-Stuff---------------
        WindModel = ofjt10 #,isotropic
        WindEfficiency = 2.0 # Default
        WindEnergyFraction = 1.0 # Default
        WindSigma0 = 353.0 #km/s Default
        WindSpeedFactor = 3.7 # Default

        WindFreeTravelLength = 20 # Default 
        WindFreeTravelDensFac = 0.1 # Not Default!
    """)

    gadget_filename = "paramfile.gadget"
    with open(gadget_filename, "w") as f:
        f.write(paramfile_gadget_content)
    print(f"paramfile.gadget has been created and saved as '{gadget_filename}' with the specified parameters.")



## Now run the pipeline with the desired parameters from terminal.
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pipeline for generating .ini files, running CLASS, and creating MP-Gadget paramfiles.")
    parser.add_argument("--class_path", type=str, required=True, help="Filepath to the CLASS binary.")
    parser.add_argument("--A_s", type=float, required=True, help="Amplitude of the primordial power spectrum (e.g., 2.1e-9).")
    parser.add_argument("--n_s", type=float, required=True, help="Spectral index of the primordial power spectrum (e.g., 0.965).")
    parser.add_argument("--h", type=float, required=True, help="Dimensionless Hubble parameter (e.g., 0.678).")
    parser.add_argument("--m_ncdm", type=float, required=True, help="Mass of the neutrino in eV (if not rescaling).")
    parser.add_argument("--log10_G_eff_nu", type=float, required=False, help="Log10 of the effective self-interaction strength of neutrinos (optional).")

    args = parser.parse_args()

    # Step 1: Write .ini file
    ini_filename, ini_filename_lcdm = write_ini_file(args.A_s, args.n_s, args.h, args.m_ncdm, args.log10_G_eff_nu)

    # Step 2: Run CLASS with the generated .ini file
    geff_pk_file, geff_tk_file = run_class(args.class_path, ini_filename,args.log10_G_eff_nu)
    lcdm_pk_file, lcdm_tk_file = run_class(args.class_path, ini_filename_lcdm,args.log10_G_eff_nu)
    lcdm_planck_pk_file  = "planck_2018_pk.dat" ## This is the file with the Planck 2018 data for P(k).

    # Step 3: Plot Pk files compared to LCDM
    plot_pk_files_matching(geff_pk_file, lcdm_pk_file, args.log10_G_eff_nu)
    plot_pk_files_matching(geff_pk_file, lcdm_planck_pk_file, args.log10_G_eff_nu, scale=False)

    # Step 4: Create paramfiles for MP-Gadget
    create_paramfile_gadget(geff_pk_file, geff_tk_file, args.A_s, args.n_s, args.h, args.m_ncdm,args.log10_G_eff_nu)

