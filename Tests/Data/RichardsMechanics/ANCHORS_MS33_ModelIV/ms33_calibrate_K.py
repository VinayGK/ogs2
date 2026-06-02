#!/usr/bin/env python3
"""
ms33_calibrate_K.py — Pi-path augmentation amplitude K calibration for ANCHORS MS33 Model I.

Calibrates vdw_augmentation_prefactor (K, J/kg) for each dry density case in
Model I by running OGS and comparing the simulated final swelling pressure at
t=200 d to the Villar empirical target:

    p_sw = exp(6.77 * rho_d[g/cm³] - 9.07)  MPa

Prerequisites (must be done before running this script):
  1. The early return bug in computeReferenceMicroPorositySwellingStressIncrement
     must be fixed (Action 3 in AGENTS.md). With slope=0 and the bug present,
     sigma_sw = 0 always and this script will never converge.
  2. OGS must be rebuilt: ninja RichardsMechanics ogs
  3. All three Model I PRJ files must have micro_water_content_swelling_slope = 0.

Usage:
    python ms33_calibrate_K.py               # calibrate all three densities
    python ms33_calibrate_K.py --verify      # verify existing K, do not change PRJ
    python ms33_calibrate_K.py --density 1400  # calibrate only one density

Outputs:
    - Updated PRJ files with calibrated K values
    - Console report: K, simulated p_sw, Villar target, % error for each density
    - ms33_calibrate_K_results.txt: machine-readable summary

Provenance: ANCHORS_MS33_ModelIV/AGENTS.md, Action 5 (Pi-path K calibration).
"""

import argparse
import os
import re
import glob
import subprocess
import sys

import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
OGS = "/Users/vinaykumar/git/build/release-omp-mfront/bin/ogs"
MODEL_I_DIR = os.path.join(
    os.path.dirname(__file__), "..", "ANCHORS_MS33_ModelI"
)
MODEL_I_DIR = os.path.realpath(MODEL_I_DIR)

CASES = {
    1400: {
        "prj": os.path.join(MODEL_I_DIR, "ms33_modelI_dd1400.prj"),
        "prefix": "ms33_modelI_dd1400",
    },
    1600: {
        "prj": os.path.join(MODEL_I_DIR, "ms33_modelI_dd1600.prj"),
        "prefix": "ms33_modelI_dd1600",
    },
    1800: {
        "prj": os.path.join(MODEL_I_DIR, "ms33_modelI_dd1800.prj"),
        "prefix": "ms33_modelI_dd1800",
    },
}

# Villar targets: p_sw = exp(6.77 * rho_d[g/cm³] - 9.07) MPa
def villar_target_pa(rho_d_kg_m3: float) -> float:
    """Villar swelling pressure target in Pa for given dry density [kg/m³]."""
    return np.exp(6.77 * (rho_d_kg_m3 / 1000.0) - 9.07) * 1e6


# ---------------------------------------------------------------------------
# PRJ editing
# ---------------------------------------------------------------------------
PRJ_ENCODING = "iso-8859-1"
K_TAG_RE = re.compile(
    r"<vdw_augmentation_prefactor>[^<]+</vdw_augmentation_prefactor>"
)


def get_K_from_prj(prj_path: str) -> float:
    """Read current K value from PRJ file."""
    with open(prj_path, encoding=PRJ_ENCODING) as f:
        content = f.read()
    m = K_TAG_RE.search(content)
    if not m:
        raise RuntimeError(f"vdw_augmentation_prefactor not found in {prj_path}")
    val_str = m.group(0).split(">")[1].split("<")[0].strip()
    return float(val_str)


def set_K_in_prj(prj_path: str, K: float) -> None:
    """Write K into PRJ file (in-place)."""
    with open(prj_path, encoding=PRJ_ENCODING) as f:
        content = f.read()
    n_matches = len(K_TAG_RE.findall(content))
    if n_matches != 1:
        raise RuntimeError(
            f"Expected exactly 1 vdw_augmentation_prefactor tag in {prj_path}, "
            f"found {n_matches}"
        )
    new_content = K_TAG_RE.sub(
        f"<vdw_augmentation_prefactor>{K:.8g}</vdw_augmentation_prefactor>",
        content,
    )
    with open(prj_path, "w", encoding=PRJ_ENCODING) as f:
        f.write(new_content)
    print(f"  [PRJ] Set K={K:.4f} J/kg in {os.path.basename(prj_path)}")


# ---------------------------------------------------------------------------
# OGS run
# ---------------------------------------------------------------------------
def run_ogs(prj_path: str, outdir: str, logfile: str, timeout_s: int = 900) -> bool:
    """Run OGS. Return True on success (returncode == 0)."""
    print(f"  [OGS] Running {os.path.basename(prj_path)} ...", flush=True)
    with open(logfile, "w") as f:
        result = subprocess.run(
            [OGS, "-o", outdir, "-l", "warn", prj_path],
            stdout=f,
            stderr=subprocess.STDOUT,
            timeout=timeout_s,
        )
    if result.returncode != 0:
        print(f"  [OGS] FAILED (returncode={result.returncode}), see {logfile}")
        return False
    print(f"  [OGS] Done → log: {os.path.basename(logfile)}")
    return True


# ---------------------------------------------------------------------------
# VTU extraction
# ---------------------------------------------------------------------------
def find_final_vtu(outdir: str, prefix: str) -> str:
    """Find VTU at t=17280000 s (= 200 days)."""
    pattern = os.path.join(outdir, f"{prefix}_ts_*_t_17280000.000000.vtu")
    files = glob.glob(pattern)
    if not files:
        raise FileNotFoundError(
            f"No VTU at t=17280000 s in {outdir} matching prefix {prefix!r}. "
            "Check that the simulation completed and the output file exists."
        )
    return sorted(files)[-1]


def extract_mean_swelling_pressure(vtu_path: str) -> float:
    """Return the spatially averaged swelling pressure [Pa] from VTU.

    In the constant-volume isotropic Model I test at t=200 d, the stress state
    is approximately isotropic.  Swelling pressure = -mean_stress:
        p_sw = -(sigma_xx + sigma_yy + sigma_zz) / 3   [Pa, positive = compressive]

    The sigma field is a Kelvin vector (2D axisymmetric):
        [sigma_xx, sigma_yy, sigma_zz, sigma_xy * sqrt(2)]
    Components 0,1,2 are the diagonal entries.

    Falls back to the swelling_stress field if sigma is not available.
    """
    try:
        import vtk
    except ImportError:
        raise ImportError(
            "vtk Python package not found. Install with: pip install vtk"
        )

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtu_path)
    reader.Update()
    mesh = reader.GetOutput()
    pd = mesh.GetPointData()

    # Prefer 'sigma' (total effective stress); fall back to 'swelling_stress'
    arr = pd.GetArray("sigma")
    if arr is None:
        arr = pd.GetArray("swelling_stress")
    if arr is None:
        raise RuntimeError(
            f"Neither 'sigma' nor 'swelling_stress' found in {vtu_path}"
        )

    n_pts = arr.GetNumberOfTuples()
    n_comp = arr.GetNumberOfComponents()
    if n_comp < 3:
        raise RuntimeError(
            f"Expected ≥3 Kelvin components in sigma, got {n_comp} in {vtu_path}"
        )

    sig = np.array([[arr.GetComponent(i, c) for c in range(n_comp)] for i in range(n_pts)])
    # Kelvin diagonal: components 0, 1, 2
    mean_stress = (sig[:, 0] + sig[:, 1] + sig[:, 2]) / 3.0
    p_sw = -float(np.mean(mean_stress))  # compressive = positive
    print(f"  [VTU] p_sw = {p_sw / 1e6:.5f} MPa  ({os.path.basename(vtu_path)})")
    return p_sw


# ---------------------------------------------------------------------------
# Calibration loop
# ---------------------------------------------------------------------------
def simulate_pressure(rho_d: int, K: float, outdir: str, n_call: list) -> float:
    """Set K, run OGS for given dry density, return simulated swelling pressure [Pa]."""
    n_call[0] += 1
    case = CASES[rho_d]
    prj_path = case["prj"]
    prefix = case["prefix"]
    logfile = os.path.join(outdir, f"calib_dd{rho_d}_call{n_call[0]}_K{K:.0f}.log")

    set_K_in_prj(prj_path, K)
    success = run_ogs(prj_path, outdir, logfile)
    if not success:
        raise RuntimeError(f"OGS failed for rho_d={rho_d}, K={K:.1f}. See {logfile}")

    vtu = find_final_vtu(outdir, prefix)
    return extract_mean_swelling_pressure(vtu)


def calibrate_one_density(rho_d: int, K_lo: float, K_hi: float,
                           tol_rel: float = 0.002, max_iter: int = 20,
                           verbose: bool = True) -> float:
    """Bisect K to match Villar target for one dry density.

    Uses scipy.optimize.brentq.  The target is the Villar swelling pressure.
    Returns the calibrated K value.
    """
    target = villar_target_pa(rho_d)
    n_call = [0]
    outdir = MODEL_I_DIR

    print(f"\n{'='*60}")
    print(f"Calibrating rho_d = {rho_d} kg/m³")
    print(f"Villar target: {target/1e6:.5f} MPa")
    print(f"Search bracket: K ∈ [{K_lo:.1f}, {K_hi:.1f}] J/kg")
    print(f"{'='*60}")

    def residual(K: float) -> float:
        p_sw = simulate_pressure(rho_d, K, outdir, n_call)
        res = p_sw - target
        if verbose:
            err_pct = res / target * 100
            print(f"    residual: {res/1e6:+.5f} MPa  ({err_pct:+.3f}%)")
        return res

    # Verify bracket: residual must change sign
    r_lo = residual(K_lo)
    r_hi = residual(K_hi)
    if r_lo * r_hi > 0:
        raise RuntimeError(
            f"Bracket [{K_lo}, {K_hi}] does not straddle zero: "
            f"residual({K_lo}) = {r_lo/1e6:.4f}, residual({K_hi}) = {r_hi/1e6:.4f}. "
            "Widen the bracket."
        )

    K_opt = brentq(residual, K_lo, K_hi, xtol=0.1, rtol=tol_rel, maxiter=max_iter)

    # Final run with optimal K
    set_K_in_prj(CASES[rho_d]["prj"], K_opt)
    p_sw_final = simulate_pressure(rho_d, K_opt, outdir, n_call)
    err_pct = abs(p_sw_final - target) / target * 100

    print(f"\n  RESULT rho_d={rho_d} kg/m³:")
    print(f"    K_opt    = {K_opt:.4f} J/kg")
    print(f"    p_sw     = {p_sw_final/1e6:.5f} MPa")
    print(f"    target   = {target/1e6:.5f} MPa")
    print(f"    error    = {err_pct:.4f}%")
    print(f"    OGS calls = {n_call[0]}")
    return K_opt


def verify_existing(densities: list) -> None:
    """Run each density with the existing K in the PRJ and report errors."""
    print("\n=== VERIFY mode: running with existing K values ===\n")
    outdir = MODEL_I_DIR
    n_call = [0]
    rows = []
    for rho_d in densities:
        case = CASES[rho_d]
        K_existing = get_K_from_prj(case["prj"])
        print(f"rho_d={rho_d}: existing K = {K_existing:.4f} J/kg")
        p_sw = simulate_pressure(rho_d, K_existing, outdir, n_call)
        target = villar_target_pa(rho_d)
        err_pct = (p_sw - target) / target * 100
        rows.append((rho_d, K_existing, p_sw, target, err_pct))

    print("\n=== Verification summary ===")
    print(f"{'rho_d':>8}  {'K':>12}  {'p_sw [MPa]':>12}  {'target [MPa]':>12}  {'err [%]':>8}")
    for rho_d, K, p_sw, tgt, err in rows:
        flag = "OK" if abs(err) < 1.0 else "FAIL"
        print(f"{rho_d:>8}  {K:>12.2f}  {p_sw/1e6:>12.5f}  {tgt/1e6:>12.5f}  {err:>+8.4f}  {flag}")

    mean_abs_err = np.mean([abs(r[4]) for r in rows])
    print(f"\nMean absolute relative error: {mean_abs_err:.4f}%")
    if mean_abs_err < 1.0:
        print("GATE PASSED: all densities < 1% MAE")
    else:
        print("GATE FAILED: recalibration required → run without --verify")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--verify", action="store_true",
                        help="Verify existing K values without changing PRJ files")
    parser.add_argument("--density", type=int, choices=[1400, 1600, 1800],
                        help="Calibrate only one dry density")
    parser.add_argument("--K-lo", type=float, default=None,
                        help="Lower bracket for K bisection (default: auto)")
    parser.add_argument("--K-hi", type=float, default=None,
                        help="Upper bracket for K bisection (default: auto)")
    args = parser.parse_args()

    densities = [args.density] if args.density else [1400, 1600, 1800]

    # Default brackets: wide enough to encompass physically realistic range
    default_brackets = {
        1400: (500.0, 5e4),
        1600: (2000.0, 2e5),
        1800: (1e4, 1e6),
    }

    if args.verify:
        verify_existing(densities)
        return

    results = {}
    for rho_d in densities:
        K_lo = args.K_lo or default_brackets[rho_d][0]
        K_hi = args.K_hi or default_brackets[rho_d][1]
        K_opt = calibrate_one_density(rho_d, K_lo, K_hi)
        results[rho_d] = K_opt

    print("\n=== Calibration complete ===")
    print(f"{'rho_d':>8}  {'K_opt [J/kg]':>16}")
    for rho_d, K_opt in results.items():
        print(f"{rho_d:>8}  {K_opt:>16.4f}")

    # Write results file
    result_path = os.path.join(
        os.path.dirname(__file__), "ms33_calibrate_K_results.txt"
    )
    with open(result_path, "w") as f:
        f.write("# ms33_calibrate_K results\n")
        f.write(f"# Villar formula: p_sw = exp(6.77*rho_d[g/cm3] - 9.07) MPa\n")
        f.write(f"# {'rho_d[kg/m3]':>15}  {'K_opt[J/kg]':>15}  {'target[MPa]':>12}\n")
        for rho_d, K_opt in results.items():
            tgt_mpa = villar_target_pa(rho_d) / 1e6
            f.write(f"  {rho_d:>15}  {K_opt:>15.6f}  {tgt_mpa:>12.6f}\n")
    print(f"\nResults written to {result_path}")

    print("\nNext steps:")
    print("  1. Propagate K_opt(dd1600) to ms33_modelI_dd1600.prj,")
    print("     ms33_modelIII_gap2mm.prj, ms33_modelIV_pellets.prj,")
    print("     ms33_modelVII_freeswelling.prj (all use rho_d=1600 reference).")
    print("  2. Rerun all ANCHORS models (see AGENTS.md Action 5).")
    print("  3. Check transport_porosity >= 0 gate for Model IV.")


if __name__ == "__main__":
    main()
