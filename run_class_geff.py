#!/usr/bin/env python3
"""
Run CLASS with Geff cosmology, generate simulation parameter files, and plot P(k) ratios.

Usage:
    python run_class_geff.py --class-path /path/to/class [--h 0.67] [--log10Geff -4.5] \
        [--m-ncdm 0.3] [--ns 0.96] [--As 2.1e-9] [--adjust-lambda]
"""

import argparse
import os
import re
import subprocess
import sys
import numpy as np
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TEMPLATE_GEFF = os.path.join(SCRIPT_DIR, "Tests", "template_geff.ini")
TEMPLATE_CDM = os.path.join(SCRIPT_DIR, "Tests", "template_cdm.ini")
TEMPLATE_GENIC = os.path.join(SCRIPT_DIR, "Paramfile_templates", "paramfile.genic")
TEMPLATE_GADGET = os.path.join(SCRIPT_DIR, "Paramfile_templates", "paramfile.gadget")
TEMPLATE_SLURM = os.path.join(SCRIPT_DIR, "Paramfile_templates", "slurm_template.sh")
CDM_REF_DIR = os.path.join(SCRIPT_DIR, "Tests", "cdm_reference")

REDSHIFTS = [99.0, 4.0, 3.0, 2.0]
# CLASS output index (1-based) corresponding to each redshift in z_pk
REDSHIFT_LABELS = {1: "z=99", 2: "z=4", 3: "z=3", 4: "z=2"}


# ---------------------------------------------------------------------------
# INI file helpers
# ---------------------------------------------------------------------------

def parse_ini(filepath):
    """Parse a CLASS ini file into an ordered list of (key, value) tuples
    and a dict for quick look-ups.  Commented-out lines are ignored."""
    params = {}
    with open(filepath) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            if "=" not in stripped:
                continue
            key_val = stripped.split("#")[0]  # strip inline comments
            key, val = key_val.split("=", 1)
            params[key.strip()] = val.strip()
    return params


def write_ini_from_template(template_path, output_path, overrides):
    """Read *template_path* line-by-line; for any key present in *overrides*
    replace its value.  Write result to *output_path*.

    Keys not in *overrides* are left exactly as in the template.
    """
    override_keys_used = set()
    lines_out = []
    with open(template_path) as f:
        for line in f:
            stripped = line.strip()
            if stripped and not stripped.startswith("#") and "=" in stripped:
                key_val = stripped.split("#")[0]
                key, _ = key_val.split("=", 1)
                key_clean = key.strip()
                if key_clean in overrides:
                    # preserve any inline comment from original line
                    comment = ""
                    if "#" in stripped:
                        comment = "  #" + stripped.split("#", 1)[1]
                    lines_out.append(f"{key_clean} = {overrides[key_clean]}{comment}\n")
                    override_keys_used.add(key_clean)
                    continue
            lines_out.append(line)

    # Append any override keys that were not found in the template
    extra = set(overrides) - override_keys_used
    if extra:
        lines_out.append("\n# Additional parameters\n")
        for k in sorted(extra):
            lines_out.append(f"{k} = {overrides[k]}\n")

    with open(output_path, "w") as f:
        f.writelines(lines_out)


# ---------------------------------------------------------------------------
# Physics helpers
# ---------------------------------------------------------------------------

def omega_ncdm_from_mass(m_ncdm_eV, T_ncdm=0.71611):
    """Approximate omega_ncdm = Omega_ncdm * h^2 for a single species.

    For the standard temperature ratio T_ncdm = (4/11)^{1/3} this reduces
    to m_ncdm / 93.14.  For non-standard T_ncdm the number density scales
    as (T_ncdm / T_standard)^3.
    """
    T_standard = (4.0 / 11.0) ** (1.0 / 3.0)
    return m_ncdm_eV / 93.14 * (T_ncdm / T_standard) ** 3


def compute_omega_radiation(T_cmb, N_ur):
    """Return omega_gamma and omega_ur (both = Omega*h^2)."""
    omega_gamma = 2.4728e-5 * (T_cmb / 2.7255) ** 4
    omega_ur = N_ur * (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0) * omega_gamma
    return omega_gamma, omega_ur


def compute_big_omegas(omega_b, omega_cdm, omega_ncdm, omega_gamma, omega_ur, h):
    """Convert little-omega -> big-Omega and derive OmegaLambda for flat universe."""
    h2 = h ** 2
    Omega_b = omega_b / h2
    Omega_cdm = omega_cdm / h2
    Omega_ncdm = omega_ncdm / h2
    Omega_gamma = omega_gamma / h2
    Omega_ur = omega_ur / h2
    Omega_m = Omega_b + Omega_cdm + Omega_ncdm
    Omega_Lambda = 1.0 - Omega_m - Omega_gamma - Omega_ur
    return Omega_m, Omega_b, Omega_cdm, Omega_ncdm, Omega_Lambda


# ---------------------------------------------------------------------------
# Parameter-file generation for MP-GenIC / MP-Gadget
# ---------------------------------------------------------------------------

def write_genic_paramfile(template_path, output_path, cosmology, tk_path, pk_path):
    """Fill cosmology values into the MP-GenIC template."""
    lines_out = []
    with open(template_path) as f:
        for line in f:
            stripped = line.strip()
            if "=" in stripped and not stripped.startswith("#"):
                key = stripped.split("=", 1)[0].strip()
                comment = ""
                if "#" in stripped:
                    comment = "  #" + stripped.split("#", 1)[1]

                if key == "Omega0":
                    lines_out.append(f"Omega0 = {cosmology['Omega_m']:.10f}{comment}\n")
                elif key == "OmegaLambda":
                    lines_out.append(f"OmegaLambda = {cosmology['Omega_Lambda']:.10f}{comment}\n")
                elif key == "OmegaBaryon":
                    lines_out.append(f"OmegaBaryon = {cosmology['Omega_b']:.10f}{comment}\n")
                elif key == "HubbleParam":
                    lines_out.append(f"HubbleParam = {cosmology['h']}{comment}\n")
                elif key == "PrimordialIndex":
                    lines_out.append(f"PrimordialIndex = {cosmology['n_s']}{comment}\n")
                elif key == "PrimordialAmp":
                    lines_out.append(f"PrimordialAmp = {cosmology['A_s']}{comment}\n")
                elif key == "FileWithTransferFunction":
                    lines_out.append(f"FileWithTransferFunction = {tk_path}{comment}\n")
                elif key == "FileWithInputSpectrum":
                    lines_out.append(f"FileWithInputSpectrum = {pk_path}{comment}\n")
                else:
                    lines_out.append(line)
            else:
                lines_out.append(line)

    with open(output_path, "w") as f:
        f.writelines(lines_out)


def write_gadget_paramfile(template_path, output_path, cosmology):
    """Fill cosmology values into the MP-Gadget template."""
    lines_out = []
    with open(template_path) as f:
        for line in f:
            stripped = line.strip()
            if "=" in stripped and not stripped.startswith("#"):
                key = stripped.split("=", 1)[0].strip()
                comment = ""
                if "#" in stripped:
                    comment = "  #" + stripped.split("#", 1)[1]

                if key == "Omega0":
                    lines_out.append(f"Omega0 = {cosmology['Omega_m']:.10f}{comment}\n")
                elif key == "OmegaLambda":
                    lines_out.append(f"OmegaLambda = {cosmology['Omega_Lambda']:.10f}{comment}\n")
                elif key == "OmegaBaryon":
                    lines_out.append(f"OmegaBaryon = {cosmology['Omega_b']:.10f}{comment}\n")
                elif key == "HubbleParam":
                    lines_out.append(f"HubbleParam = {cosmology['h']}{comment}\n")
                else:
                    lines_out.append(line)
            else:
                lines_out.append(line)

    with open(output_path, "w") as f:
        f.writelines(lines_out)


def write_slurm_script(template_path, output_path, run_name, genic_filename, gadget_filename):
    """Fill job name and file references into SLURM template."""
    lines_out = []
    with open(template_path) as f:
        for line in f:
            if line.startswith("#SBATCH --job-name="):
                lines_out.append(f"#SBATCH --job-name={run_name}\n")
            elif "MP-GenIC" in line:
                lines_out.append(line.replace("_genic_params.ini", genic_filename))
            elif "MP-Gadget" in line:
                lines_out.append(line.replace("mpgadget.param", gadget_filename))
            else:
                lines_out.append(line)

    with open(output_path, "w") as f:
        f.writelines(lines_out)


# ---------------------------------------------------------------------------
# CLASS runner
# ---------------------------------------------------------------------------

def run_class(class_exe, ini_path):
    """Execute CLASS and return stdout+stderr."""
    print(f"  Running CLASS: {class_exe} {ini_path}")
    result = subprocess.run(
        [class_exe, ini_path],
        capture_output=True, text=True,
    )
    if result.returncode != 0:
        print("CLASS STDERR:\n", result.stderr, file=sys.stderr)
        print("CLASS STDOUT:\n", result.stdout, file=sys.stderr)
        sys.exit(f"CLASS failed with exit code {result.returncode}")
    return result.stdout + result.stderr


def find_class_outputs(root_prefix, run_dir):
    """Discover the CLASS output files matching root_prefix in run_dir.

    CLASS appends a two-digit counter (00, 01, …) to avoid overwriting.
    We find the latest counter for each (z-index, type) combination.
    Returns dict: {z_index: {"pk": path, "tk": path, "pk_cb": path}}
    """
    outputs = {}
    pattern = re.compile(
        re.escape(os.path.basename(root_prefix))
        + r"(\d{2})_z(\d+)_(pk|pk_cb|tk)\.dat$"
    )
    for fname in sorted(os.listdir(run_dir)):
        m = pattern.match(fname)
        if m:
            counter, zi, ftype = m.group(1), int(m.group(2)), m.group(3)
            full = os.path.join(run_dir, fname)
            if zi not in outputs:
                outputs[zi] = {}
            outputs[zi][ftype] = full
    return outputs


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_pk_ratios(geff_outputs, cdm_outputs, run_dir, h_value):
    """Plot P(k)_geff / P(k)_cdm at each redshift and save to run_dir."""
    fig, ax = plt.subplots(figsize=(9, 6))

    colors = {1: "C0", 2: "C1", 3: "C2", 4: "C3"}
    for zi in sorted(geff_outputs.keys()):
        if "pk" not in geff_outputs[zi] or zi not in cdm_outputs or "pk" not in cdm_outputs[zi]:
            continue

        geff_data = np.loadtxt(geff_outputs[zi]["pk"], skiprows=4)
        cdm_data = np.loadtxt(cdm_outputs[zi]["pk"], skiprows=4)

        k_geff, pk_geff = geff_data[:, 0], geff_data[:, 1]
        k_cdm, pk_cdm = cdm_data[:, 0], cdm_data[:, 1]

        interp_cdm = interp1d(
            np.log10(k_cdm), np.log10(pk_cdm),
            kind="cubic", bounds_error=False, fill_value="extrapolate",
        )
        pk_cdm_interp = 10.0 ** interp_cdm(np.log10(k_geff))

        label = REDSHIFT_LABELS.get(zi, f"z-index {zi}")
        ax.plot(k_geff, pk_geff / pk_cdm_interp, label=label, color=colors.get(zi))

    ax.set_xscale("log")
    ax.set_xlabel(r"$k$ [$h$/Mpc]")
    ax.set_ylabel(r"$P(k)_{\rm Geff} / P(k)_{\rm CDM}$")
    ax.set_title("Matter Power Spectrum Ratio (Geff / CDM)")
    ax.axvline(x=1.0 / h_value, color="gray", ls="--", alpha=0.5, label=r"$k = 1/h$")
    ax.set_xlim(1e-3, 1e2)
    ax.set_ylim(0.7, 1.2)
    ax.legend()
    ax.grid(True, alpha=0.3)

    outpath = os.path.join(run_dir, "Pk_ratio_geff_to_cdm.png")
    fig.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Plot saved: {outpath}")


# ---------------------------------------------------------------------------
# CDM reference generation
# ---------------------------------------------------------------------------

def ensure_cdm_reference(class_exe):
    """Generate CDM reference P(k) and T(k) if they don't already exist."""
    if os.path.isdir(CDM_REF_DIR) and any(f.endswith("_pk.dat") for f in os.listdir(CDM_REF_DIR)):
        print("  CDM reference files already exist, skipping.")
        return find_class_outputs("cdm_ref_", CDM_REF_DIR)

    print("  Generating CDM reference from template_cdm.ini ...")
    os.makedirs(CDM_REF_DIR, exist_ok=True)

    root_prefix = os.path.join(CDM_REF_DIR, "cdm_ref_")
    overrides = {"root": root_prefix}
    ini_path = os.path.join(CDM_REF_DIR, "cdm_reference.ini")
    write_ini_from_template(TEMPLATE_CDM, ini_path, overrides)

    run_class(class_exe, ini_path)
    outputs = find_class_outputs("cdm_ref_", CDM_REF_DIR)
    if not outputs:
        sys.exit("ERROR: CLASS ran but produced no output files for CDM reference.")
    print(f"  CDM reference files generated in {CDM_REF_DIR}")
    return outputs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_run_name(log10Geff, m_ncdm, As):
    """Construct a directory name from key parameter values."""
    geff_str = f"Geff{log10Geff:.3f}"
    mncdm_str = f"_mncdm{m_ncdm:.4f}"
    as_str = f"_As{As:.4e}"
    return (geff_str + mncdm_str + as_str).replace("+", "")


def main():
    parser = argparse.ArgumentParser(
        description="Run CLASS with Geff cosmology and generate MP-Gadget parameter files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--class-path", type=str,
        default="/Users/verag/repositories/snu-class-update/class_snu_uptodate/class",
        help="Path to CLASS executable (default: snu-class-update)",
    )
    parser.add_argument("--h", type=float, default=None, dest="h_param",
                        help="Hubble parameter h")
    parser.add_argument("--log10Geff", type=float, default=None,
                        help="log10(Geff / MeV^-2)")
    parser.add_argument("--m-ncdm", type=float, default=None,
                        help="Mass of ncdm species [eV]")
    parser.add_argument("--ns", type=float, default=None,
                        help="Scalar spectral index n_s")
    parser.add_argument("--As", type=float, default=None,
                        help="Scalar amplitude A_s")
    parser.add_argument("--adjust-lambda", action="store_true",
                        help="Compensate flatness via Omega_Lambda instead of omega_cdm")
    args = parser.parse_args()

    class_exe = os.path.abspath(args.class_path)
    if not os.path.isfile(class_exe):
        sys.exit(f"CLASS executable not found: {class_exe}")

    # ------------------------------------------------------------------
    # 1. Read template defaults
    # ------------------------------------------------------------------
    tmpl = parse_ini(TEMPLATE_GEFF)

    h = args.h_param if args.h_param is not None else float(tmpl["h"])
    log10Geff = args.log10Geff if args.log10Geff is not None else float(tmpl["log10_G_eff_nu"])
    m_ncdm = args.m_ncdm if args.m_ncdm is not None else float(tmpl["m_ncdm"])
    n_s = args.ns if args.ns is not None else float(tmpl["n_s"])
    A_s = args.As if args.As is not None else float(tmpl["A_s"])
    omega_cdm_tmpl = float(tmpl["omega_cdm"])
    omega_b = float(tmpl["omega_b"])
    T_ncdm = float(tmpl["T_ncdm"])
    T_cmb = float(tmpl["T_cmb"])
    N_ur = float(tmpl["N_ur"])
    m_ncdm_tmpl = float(tmpl["m_ncdm"])

    # ------------------------------------------------------------------
    # 2. Flatness: adjust omega_cdm or let Lambda float
    # ------------------------------------------------------------------
    omega_ncdm_tmpl = omega_ncdm_from_mass(m_ncdm_tmpl, T_ncdm)
    omega_ncdm_new = omega_ncdm_from_mass(m_ncdm, T_ncdm)
    delta_omega_ncdm = omega_ncdm_new - omega_ncdm_tmpl

    if args.adjust_lambda:
        omega_cdm = omega_cdm_tmpl
        print(f"  --adjust-lambda: keeping omega_cdm = {omega_cdm:.6f}, "
              f"CLASS will adjust Omega_Lambda for flatness.")
    else:
        omega_cdm = omega_cdm_tmpl - delta_omega_ncdm
        print(f"  Flatness adjustment: omega_ncdm changed by {delta_omega_ncdm:+.6f}")
        print(f"  omega_cdm: {omega_cdm_tmpl:.6f} -> {omega_cdm:.6f}")

    # ------------------------------------------------------------------
    # 3. Derived big-Omega quantities (for param files)
    # ------------------------------------------------------------------
    omega_gamma, omega_ur = compute_omega_radiation(T_cmb, N_ur)
    Omega_m, Omega_b, Omega_cdm, Omega_ncdm, Omega_Lambda = compute_big_omegas(
        omega_b, omega_cdm, omega_ncdm_new, omega_gamma, omega_ur, h,
    )
    print(f"  Omega_m = {Omega_m:.6f}  Omega_Lambda = {Omega_Lambda:.6f}  "
          f"Omega_b = {Omega_b:.6f}  h = {h}")

    # ------------------------------------------------------------------
    # 4. Create run directory
    # ------------------------------------------------------------------
    run_name = build_run_name(log10Geff, m_ncdm, A_s)
    run_dir = os.path.join(SCRIPT_DIR, "runs", run_name)
    os.makedirs(run_dir, exist_ok=True)
    print(f"  Run directory: {run_dir}")

    # ------------------------------------------------------------------
    # 5. Write CLASS ini file (only override what user specified + omega_cdm + root)
    # ------------------------------------------------------------------
    root_prefix = os.path.join(run_dir, "geff_")
    overrides = {"root": root_prefix}

    if args.h_param is not None:
        overrides["h"] = str(h)
    if args.log10Geff is not None:
        overrides["log10_G_eff_nu"] = str(log10Geff)
    if args.m_ncdm is not None:
        overrides["m_ncdm"] = str(m_ncdm)
    if args.ns is not None:
        overrides["n_s"] = str(n_s)
    if args.As is not None:
        overrides["A_s"] = str(A_s)

    if not args.adjust_lambda:
        overrides["omega_cdm"] = f"{omega_cdm:.10f}"

    ini_path = os.path.join(run_dir, "class_geff.ini")
    write_ini_from_template(TEMPLATE_GEFF, ini_path, overrides)
    print(f"  CLASS ini written: {ini_path}")

    # ------------------------------------------------------------------
    # 6. Run CLASS
    # ------------------------------------------------------------------
    run_class(class_exe, ini_path)
    geff_outputs = find_class_outputs("geff_", run_dir)
    if not geff_outputs:
        sys.exit("ERROR: CLASS ran but produced no output files.")
    print(f"  CLASS outputs: {len(geff_outputs)} redshift(s)")

    # ------------------------------------------------------------------
    # 7. Generate CDM reference (once)
    # ------------------------------------------------------------------
    cdm_outputs = ensure_cdm_reference(class_exe)

    # ------------------------------------------------------------------
    # 8. Generate MP-GenIC, MP-Gadget, SLURM param files
    # ------------------------------------------------------------------
    cosmology = {
        "h": h, "n_s": n_s, "A_s": A_s,
        "Omega_m": Omega_m, "Omega_b": Omega_b, "Omega_Lambda": Omega_Lambda,
    }

    # z=99 files for ICs
    z99_tk = geff_outputs.get(1, {}).get("tk", "")
    z99_pk = geff_outputs.get(1, {}).get("pk", "")

    genic_out = os.path.join(run_dir, "paramfile.genic")
    write_genic_paramfile(TEMPLATE_GENIC, genic_out, cosmology, z99_tk, z99_pk)
    print(f"  Genic params: {genic_out}")

    gadget_out = os.path.join(run_dir, "paramfile.gadget")
    write_gadget_paramfile(TEMPLATE_GADGET, gadget_out, cosmology)
    print(f"  Gadget params: {gadget_out}")

    slurm_out = os.path.join(run_dir, "run.sh")
    write_slurm_script(TEMPLATE_SLURM, slurm_out, run_name, "paramfile.genic", "paramfile.gadget")
    print(f"  SLURM script: {slurm_out}")

    # ------------------------------------------------------------------
    # 9. Plot P(k) ratios
    # ------------------------------------------------------------------
    plot_pk_ratios(geff_outputs, cdm_outputs, run_dir, h)

    print("\nDone. All outputs in:", run_dir)


if __name__ == "__main__":
    main()
