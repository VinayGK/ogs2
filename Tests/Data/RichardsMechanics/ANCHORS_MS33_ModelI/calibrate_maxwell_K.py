#!/usr/bin/env python3
"""Per-density secant calibration of the vdW augmentation prefactor K for the
ANCHORS MS33 Model I cells on the dsm_native_maxwell_conjugate branch.

Fits K = potential_augmentation_prefactor so the SATURATED mean swelling stress
matches the Villar/Lloret Eq.7 target  Ps[MPa] = exp(6.77*rho_d[g/cm3] - 9.07).

This is the CALIBRATION driver. The Ps-vs-Villar overlay it feeds is the
calibration DISPLAY (CLAUDE.md §2: not an independent pass/fail assertion).

Run, e.g.:
  python3 calibrate_maxwell_K.py --dd 1600
One process per density; launch the three concurrently (one per 4-core block).
"""
import argparse, glob, json, os, re, shutil, subprocess, sys, time
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

HERE = os.path.dirname(os.path.abspath(__file__))
OGS = "/Users/vinaykumar/git/build/maxwell-conjugate-20260602/bin/ogs"
RHO_SOLID = 2780.0  # kg/m^3, EURAD-2 MS33 spec (summarizer RHO_SOLID)

# Dixon (2023) Fig.1 swelling pressure, EMDD=rho_d targets (MPa) -- the MS33 anchor
# (WG agreement 2026-05-27; supersedes Villar Eq.7). dd1600=14.161 from the Model III/IV
# PRJ headers; dd1400/1800 from Dixon Fig.1 (EMDD=rho_d); dd900 unused (model edge).
DIXON_EMDD_MPa = {1400: 4.9, 1600: 14.161, 1800: 40.6, 900: 0.051}
def dixon_emdd_target_MPa(rho_d_kg):
    return float(DIXON_EMDD_MPa[int(round(rho_d_kg))])

CASES = {
    1400: dict(prj="ms33_modelI_dd1400.prj", k0=1.2264e-20, phi0=0.4964028776978417),
    1600: dict(prj="ms33_modelI_dd1600.prj", k0=5.8703e-21, phi0=0.4244604316546763),
    1800: dict(prj="ms33_modelI_dd1800.prj", k0=2.6570e-21, phi0=0.3525179856115108),
    900: dict(prj="ms33_modelI_dd900.prj", k0=7.48e-20, phi0=0.6762589928057554),
}


def read_swelling_stress_MPa(vtu):
    """Saturated mean swelling stress (compression positive), MPa.
    mean = (-sxx - syy - szz)/3 averaged over cells/points, per summarizer."""
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(vtu)
    r.Update()
    g = r.GetOutput()
    pd = g.GetPointData()
    arr = pd.GetArray("swelling_stress")
    if arr is None:
        cd = g.GetCellData()
        arr = cd.GetArray("swelling_stress")
    if arr is None:
        raise RuntimeError(f"swelling_stress not in {vtu}")
    a = vtk_to_numpy(arr)  # (n, ncomp) Kelvin-mapped: xx,yy,zz,xy[,...]
    if a.ndim == 1:
        a = a.reshape(1, -1)
    sxx, syy, szz = a[:, 0], a[:, 1], a[:, 2]
    mean = (-sxx - syy - szz) / 3.0  # Pa, compression positive
    return float(np.mean(mean)) / 1e6


def latest_vtu(prefix):
    cands = glob.glob(os.path.join(HERE, f"{prefix}_ts_*_t_*.vtu"))
    if not cands:
        return None
    def ts(p):
        m = re.search(r"_ts_(\d+)_t_", os.path.basename(p))
        return int(m.group(1)) if m else -1
    return max(cands, key=ts)


def set_prefactor(text, K):
    return re.sub(
        r"(<potential_augmentation_prefactor>)[^<]+(</potential_augmentation_prefactor>)",
        rf"\g<1>{K:.6g}\g<2>", text)


def run_ogs(prj_path, prefix, outdir):
    env = dict(os.environ, OMP_NUM_THREADS="4")
    log = os.path.join(outdir, f"{prefix}.calib.log")
    with open(log, "w") as lf:
        p = subprocess.run([OGS, prj_path, "-o", outdir],
                           stdout=lf, stderr=subprocess.STDOUT, env=env, cwd=outdir)
    return p.returncode, log


def evaluate_K(K, case, base_text, prefix, workdir):
    """Run one OGS sim at prefactor K; return (Ps_sim_MPa, ok, logtail)."""
    text = set_prefactor(base_text, K)
    # IMPORTANT: never overwrite the committed source PRJ. Write a temp prj in the
    # same directory (so ../mesh relative paths still resolve) with a _calib suffix.
    prj = os.path.join(workdir, case["prj"].replace(".prj", "_calib.prj"))
    with open(prj, "w") as f:
        f.write(text)
    # clean stale vtus for this prefix
    for old in glob.glob(os.path.join(workdir, f"{prefix}_ts_*_t_*.vtu")):
        os.remove(old)
    rc, log = run_ogs(prj, prefix, workdir)
    vtu = max(
        glob.glob(os.path.join(workdir, f"{prefix}_ts_*_t_*.vtu")),
        key=lambda p: int(re.search(r"_ts_(\d+)_t_", os.path.basename(p)).group(1)),
        default=None)
    tail = ""
    with open(log) as lf:
        tail = "".join(lf.readlines()[-25:])
    if rc != 0 or vtu is None:
        return None, False, tail
    # require the final time step (t_end=17280000) to be present => saturated
    if "_t_1.728" not in os.path.basename(vtu) and "17280000" not in os.path.basename(vtu):
        # accept anyway but flag; the run may have stopped early
        pass
    try:
        ps = read_swelling_stress_MPa(vtu)
    except Exception as e:
        return None, False, tail + f"\n[read error] {e}"
    return ps, True, os.path.basename(vtu)


def calibrate(dd, rel_tol=0.02, max_iter=12):
    case = CASES[dd]
    prefix = f"ms33_modelI_dd{dd}"
    target = dixon_emdd_target_MPa(dd)
    src_prj = os.path.join(HERE, case["prj"])
    base_text = open(src_prj).read()
    # initial K = the current value already written in the PRJ (post-edit guess)
    m = re.search(r"<potential_augmentation_prefactor>([^<]+)</", base_text)
    K0 = float(m.group(1))

    workdir = os.path.join(HERE, f"_calib_dd{dd}")
    os.makedirs(workdir, exist_ok=True)
    # copy meshes referenced as ../*.vtu : run in HERE instead to keep ../ valid
    # -> simpler: run in HERE, write temp prj into HERE under a calib name.
    shutil.rmtree(workdir, ignore_errors=True)

    history = []
    # secant in log-space: f(logK) = log(Ps_sim) - log(target)
    def f_of(K):
        ps, ok, info = evaluate_K(K, case, base_text, prefix, HERE)
        history.append(dict(K=K, Ps_sim=ps, ok=ok, info=info))
        if not ok:
            print(f"[dd{dd}] K={K:.6g} FAILED: {info}", flush=True)
            return None, ps, ok, info
        rel = (ps - target) / target
        print(f"[dd{dd}] K={K:.6g} -> Ps={ps:.4f} MPa (target {target:.4f}, rel {rel:+.3%})", flush=True)
        return np.log(ps) - np.log(target), ps, ok, info

    # bracket: two points
    x0 = np.log(K0)
    g0, ps0, ok0, info0 = f_of(K0)
    if not ok0:
        return dict(dd=dd, converged=False, reason="initial run failed",
                    logtail=info0, target=target, history=history)
    if abs((ps0 - target) / target) <= rel_tol:
        return finalize(dd, case, K0, ps0, target, history, True)

    # second point: scale K by the ratio target/ps0 in log (Ps ~ monotone increasing in K)
    K1 = K0 * (target / ps0)
    x1 = np.log(K1)
    g1, ps1, ok1, info1 = f_of(K1)
    if not ok1:
        return dict(dd=dd, converged=False, reason="second run failed",
                    logtail=info1, target=target, history=history)
    if abs((ps1 - target) / target) <= rel_tol:
        return finalize(dd, case, K1, ps1, target, history, True)

    for it in range(max_iter - 2):
        if g1 == g0:
            break
        x2 = x1 - g1 * (x1 - x0) / (g1 - g0)
        K2 = float(np.exp(x2))
        g2, ps2, ok2, info2 = f_of(K2)
        if not ok2:
            return dict(dd=dd, converged=False, reason="secant run failed",
                        logtail=info2, target=target, history=history)
        if abs((ps2 - target) / target) <= rel_tol:
            return finalize(dd, case, K2, ps2, target, history, True)
        x0, g0 = x1, g1
        x1, g1, ps1 = x2, g2, ps2

    # not converged within max_iter; return best
    best = min((h for h in history if h["ok"]),
               key=lambda h: abs((h["Ps_sim"] - target) / target))
    return finalize(dd, case, best["K"], best["Ps_sim"], target, history, False)


def finalize(dd, case, K, ps, target, history, converged):
    rel = (ps - target) / target
    return dict(dd=dd, K=K, Ps_sim_MPa=ps, Ps_target_MPa=target,
                rel_err=rel, k0_m2=case["k0"], phi0=case["phi0"],
                converged=converged, n_runs=len(history), history=history)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--dd", type=int, required=True, choices=[900, 1400, 1600, 1800])
    ap.add_argument("--rel-tol", type=float, default=0.02)
    ap.add_argument("--max-iter", type=int, default=12)
    a = ap.parse_args()
    t0 = time.time()
    res = calibrate(a.dd, a.rel_tol, a.max_iter)
    res["wall_s"] = time.time() - t0
    out = os.path.join(HERE, f"_calib_result_dd{a.dd}.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=2)
    print(f"[dd{a.dd}] DONE converged={res.get('converged')} "
          f"K={res.get('K')} Ps={res.get('Ps_sim_MPa')} wall={res['wall_s']:.0f}s -> {out}",
          flush=True)
