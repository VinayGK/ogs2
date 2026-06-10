#!/usr/bin/env python3
# Model VII free-swelling + load/unload: film_strain_coupling off vs kinematic
# vs equilibrium. Extracts per output VTU: mean porosity -> void ratio
# e = phi/(1-phi), mean micro water content n_l, top axial displacement, and
# the applied axial stress from the PRJ loading schedule (see PRJ header
# comment; traction is the curve value, compression negative).
# Outputs: filmstrain_compare.csv, e_vs_time.png, e_vs_sigma.png, summary.txt
import glob
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meshio
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
VARIANTS = ["off", "kinematic", "equilibrium"]
COLORS = {"off": "0.45", "kinematic": "tab:blue", "equilibrium": "tab:red"}

# Axial stress schedule [Pa] from the PRJ header (loading/unloading path);
# compression positive here for plotting.
SCHEDULE = [  # (t_start, t_end, sigma_start, sigma_end)
    (0.0, 17280000.0, 0.2e6, 0.2e6),
    (17280000.0, 17712000.0, 0.2e6, 0.2e6),
    (17712000.0, 18144000.0, 0.2e6, 0.4e6),
    (18144000.0, 18576000.0, 0.4e6, 1.0e6),
    (18576000.0, 19008000.0, 1.0e6, 2.5e6),
    (19008000.0, 19440000.0, 2.5e6, 5.0e6),
    (19440000.0, 19872000.0, 5.0e6, 2.5e6),
    (19872000.0, 20304000.0, 2.5e6, 1.0e6),
    (20304000.0, 20736000.0, 1.0e6, 0.4e6),
]


def sigma_axial(t):
    for t0, t1, s0, s1 in SCHEDULE:
        if t <= t1:
            return s0 + (s1 - s0) * (t - t0) / max(t1 - t0, 1.0)
    return SCHEDULE[-1][3]


def extract(variant):
    rows = []
    for f in glob.glob(os.path.join(HERE, f"out_{variant}", "*.vtu")):
        m = re.search(r"_ts_(\d+)_t_([0-9.]+)\.vtu$", f)
        if not m:
            continue
        ts, t = int(m.group(1)), float(m.group(2))
        mesh = meshio.read(f)
        cd = {k: v[0] for k, v in mesh.cell_data.items()}
        pd = mesh.point_data
        phi = float(np.mean(cd["porosity_avg"])) if "porosity_avg" in cd else (
            float(np.mean(pd["porosity"])) if "porosity" in pd else np.nan)
        nl_keys = [k for k in list(cd) + list(pd)
                   if "micro_water_content" in k]
        if nl_keys:
            k = nl_keys[0]
            nl = float(np.mean(cd[k] if k in cd else pd[k]))
        else:
            nl = np.nan
        uz_top = float(np.max(mesh.points[:, 1] * 0 +
                              pd["displacement"][:, 1])) \
            if "displacement" in pd else np.nan
        rows.append((ts, t, phi, nl, uz_top))
    rows.sort()
    return np.array(rows)


def main():
    data = {}
    for v in VARIANTS:
        d = extract(v)
        if d.size == 0:
            print(f"WARNING: no outputs for {v}")
            continue
        data[v] = d

    with open(os.path.join(HERE, "filmstrain_compare.csv"), "w") as fh:
        fh.write("variant,ts,t_s,t_d,phi_mean,e_void,n_l_mean,uz_top_m,"
                 "sigma_axial_Pa\n")
        for v, d in data.items():
            for ts, t, phi, nl, uz in d:
                e = phi / (1.0 - phi)
                fh.write(f"{v},{int(ts)},{t},{t/86400.0},{phi},{e},{nl},"
                         f"{uz},{sigma_axial(t)}\n")

    # e(t)
    fig, ax = plt.subplots(figsize=(8, 5))
    for v, d in data.items():
        e = d[:, 2] / (1.0 - d[:, 2])
        ax.plot(d[:, 1] / 86400.0, e, "-", color=COLORS[v], label=v, lw=1.8)
    ax.set_xlabel("time [d]")
    ax.set_ylabel("void ratio e = phi/(1-phi) [-]")
    ax.axvline(200.0, color="k", lw=0.6, ls=":")
    ax.text(200.5, ax.get_ylim()[0], "load/unload", fontsize=8, rotation=90,
            va="bottom")
    ax.legend(title="film_strain_coupling")
    ax.set_title("MS33 Model VII free swelling + load/unload "
                 "(dsm_native_h_of_eps @7ff8861847)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "e_vs_time.png"), dpi=160)

    # e vs sigma (the inter-team overlay shape for VII)
    fig, ax = plt.subplots(figsize=(8, 5))
    for v, d in data.items():
        e = d[:, 2] / (1.0 - d[:, 2])
        sig = np.array([sigma_axial(t) for t in d[:, 1]]) / 1e6
        mask = d[:, 1] >= 17280000.0  # the mechanical path
        ax.plot(sig[mask], e[mask], "o-", color=COLORS[v], label=v, ms=3,
                lw=1.4)
    ax.set_xscale("log")
    ax.set_xlabel("axial stress [MPa]")
    ax.set_ylabel("void ratio e [-]")
    ax.legend(title="film_strain_coupling")
    ax.set_title("Model VII load/unload path (t >= 200 d)")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "e_vs_sigma.png"), dpi=160)

    with open(os.path.join(HERE, "summary.txt"), "w") as fh:
        fh.write("variant   e@200d(end swelling)  e@5MPa(peak load)  "
                 "e@end(unload 0.4MPa)  n_l@end\n")
        for v, d in data.items():
            e = d[:, 2] / (1.0 - d[:, 2])
            t = d[:, 1]

            def at(tt):
                i = int(np.argmin(np.abs(t - tt)))
                return e[i]

            fh.write(f"{v:12s}  {at(17280000.0):8.4f}  "
                     f"{at(19440000.0):8.4f}  {at(20736000.0):8.4f}  "
                     f"{d[-1, 3]:.6g}\n")
    print(open(os.path.join(HERE, "summary.txt")).read())


if __name__ == "__main__":
    main()
