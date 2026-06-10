#!/usr/bin/env python3
"""
FEBEX native-DSM K calibration (Vinay approved 2026-06-01).
Fit vdw_augmentation_prefactor K so the single-element constant-volume swelling
pressure hits the FEBEX density-compliant target (Villar/ENRESA, b=6.77):
  block  rho_d=1.60 -> Ps = 5.82 MPa
  pellet rho_d=1.30 -> Ps = 0.76 MPa
sigma_swell is K-governed (BishopsSaturationCutoff=1, chi=0) -> macro WRC inert;
specific_surface=725 (FEBEX) DOES enter the micro film thickness, so K is NOT the
MX-80 value. We FIT it; the swelling pressure is the calibration target (NOT an
independent assertion) -> CLAUDE.md §2-clean.

FEBEX params (spec/FEBEX_PARAMETER_PROPOSAL.md, all verified):
  rho_s=2700 (Gs 2.70, D5.2.1:4264), specific_surface=725 m2/g (D5.2.1:4269),
  A_Hamaker=2.2e-20 (Israelachvili-Adams 1978), decay 7.5e-7 (prior-commit),
  micro EOS a=1e-16/b=1/rho_l0=100 (MS33 carry, Vinay 2026-06-01).
Calibration runs use the native SaturationTuller macro (inert here); the MGR
column run will swap in the FEBEX vG macro WRC.
"""
import os, re, subprocess, glob, sys
import numpy as np, meshio

ROOT = "/Users/vinaykumar/git/ogs/beacon_wp5_mgr_repro_2026-06-01"
MODEL = ROOT + "/model"
OGS = "/Users/vinaykumar/git/build/dsm-native-omp-release/bin/ogs"

MAT = {
    "block":  dict(rho_d=1600, phi0=0.4074074074, n_s=0.5925925926,
                   k0=5.0e-21, suction=1.2e8, target=5.82, K0=25000.0),
    "pellet": dict(rho_d=1300, phi0=0.5185185185, n_s=0.4814814815,
                   k0=1.3e-19, suction=1.2e8, target=0.76, K0=3300.0),
}

PRJ = """<?xml version="1.0" encoding="ISO-8859-1"?>
<!-- FEBEX K-calibration ({mat}, rho_d={rho_d}) — single-element constant volume.
     K fitted to FEBEX Ps target {target} MPa (Villar/ENRESA b=6.77). Autonomy folder. -->
<OpenGeoSysProject>
    <meshes>
        <mesh axially_symmetric="true">square_1e-2_quad_1e0.vtu</mesh>
        <mesh axially_symmetric="true">square_1e-2_quad_1e0_left.vtu</mesh>
        <mesh axially_symmetric="true">square_1e-2_quad_1e0_right.vtu</mesh>
        <mesh axially_symmetric="true">square_1e-2_quad_1e0_top.vtu</mesh>
        <mesh axially_symmetric="true">square_1e-2_quad_1e0_bottom.vtu</mesh>
    </meshes>
    <processes>
        <process>
            <name>RM</name><type>RICHARDS_MECHANICS</type><integration_order>2</integration_order>
            <micro_porosity>
                <mass_exchange_coefficient>1e-13</mass_exchange_coefficient>
                <nonlinear_solver><maximum_iterations>100</maximum_iterations>
                    <residuum_tolerance>1e-8</residuum_tolerance><increment_tolerance>1e-20</increment_tolerance></nonlinear_solver>
            </micro_porosity>
            <potential_exchange>
                <enabled>true</enabled><pressure_tolerance>1e-12</pressure_tolerance>
                <hamaker_constant>2.2e-20</hamaker_constant>
                <specific_surface>725.0</specific_surface>
                <micro_solid_density_reference>2700.0</micro_solid_density_reference>
                <micro_solid_volume_fraction_reference>{n_s}</micro_solid_volume_fraction_reference>
                <micro_liquid_density_reference>100.0</micro_liquid_density_reference>
                <micro_liquid_density_a>1e-16</micro_liquid_density_a>
                <micro_liquid_density_b>1.0</micro_liquid_density_b>
                <initial_micro_water_content>1.1699e-3</initial_micro_water_content>
                <local_nonlinear_solve_mode>scalar_micro_macro_mass_storage_mode</local_nonlinear_solve_mode>
                <micro_potential_convention>negative_attractive</micro_potential_convention>
                <use_micro_liquid_density_for_micro_pressure>true</use_micro_liquid_density_for_micro_pressure>
                <accumulate_swelling_contributions>true</accumulate_swelling_contributions>
                <vdw_augmentation_prefactor>{K}</vdw_augmentation_prefactor>
                <vdw_augmentation_decay_length>7.5e-7</vdw_augmentation_decay_length>
            </potential_exchange>
            <constitutive_relation><type>LinearElasticIsotropic</type>
                <youngs_modulus>E</youngs_modulus><poissons_ratio>nu</poissons_ratio></constitutive_relation>
            <process_variables><pressure>pressure</pressure><displacement>displacement</displacement></process_variables>
            <secondary_variables>
                <secondary_variable name="sigma"/><secondary_variable name="swelling_stress"/>
                <secondary_variable name="saturation"/><secondary_variable name="porosity"/>
                <secondary_variable name="micro_pressure"/><secondary_variable name="micro_water_content"/>
            </secondary_variables>
            <specific_body_force>0 0</specific_body_force>
            <initial_stress>sigma0</initial_stress><mass_lumping>true</mass_lumping>
        </process>
    </processes>
    <media><medium id="0"><phases>
        <phase><type>AqueousLiquid</type><properties>
            <property><name>viscosity</name><type>Constant</type><value>1e-3</value></property>
            <property><name>density</name><type>Constant</type><value>1000.0</value></property></properties></phase>
        <phase><type>Solid</type><properties>
            <property><name>density</name><type>Constant</type><value>2700.0</value></property></properties></phase>
    </phases><properties>
        <property><name>biot_coefficient</name><type>Constant</type><value>1.0</value></property>
        <property><name>permeability</name><type>KozenyCarman</type>
            <initial_permeability>k0</initial_permeability><initial_porosity>phi0</initial_porosity></property>
        <property><name>porosity</name><type>PorosityFromMassBalance</type>
            <initial_porosity>phi0</initial_porosity><minimal_porosity>0</minimal_porosity><maximal_porosity>1</maximal_porosity></property>
        <property><name>reference_temperature</name><type>Constant</type><value>293.15</value></property>
        <property><name>relative_permeability</name><type>RelativePermeabilityGeneralizedPower</type>
            <residual_liquid_saturation>0.0</residual_liquid_saturation><residual_gas_saturation>0.0</residual_gas_saturation>
            <min_relative_permeability>1e-2</min_relative_permeability><enhancement_factor>1.0</enhancement_factor><exponent>3</exponent></property>
        <property><name>saturation</name><type>SaturationTuller</type>
            <residual_liquid_saturation>0.0</residual_liquid_saturation><residual_gas_saturation>0.0</residual_gas_saturation>
            <area_factor_tuller>1.0</area_factor_tuller><pore_area_shapefactor_tuller>0.8584073464102069</pore_area_shapefactor_tuller>
            <characteristic_pore_size>1e-5</characteristic_pore_size><surface_tension>0.0715</surface_tension>
            <cavitation_pressure>1.4e8</cavitation_pressure></property>
        <property><name>bishops_effective_stress</name><type>BishopsSaturationCutoff</type><cutoff_value>1</cutoff_value></property>
    </properties></medium></media>
    <time_loop><processes><process ref="RM">
        <nonlinear_solver>basic_newton</nonlinear_solver>
        <convergence_criterion><type>PerComponentDeltaX</type><norm_type>NORM2</norm_type>
            <abstols>1 1e-12 1e-12</abstols><reltols>1e-6 1e-6 1e-6</reltols></convergence_criterion>
        <time_discretization><type>BackwardEuler</type></time_discretization>
        <time_stepping><type>IterationNumberBasedTimeStepping</type>
            <t_initial>0.0</t_initial><t_end>3456000</t_end>
            <initial_dt>10</initial_dt><minimum_dt>0.1</minimum_dt><maximum_dt>86400</maximum_dt>
            <number_iterations>1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25</number_iterations>
            <multiplier>1.225 1.2 1.175 1.15 1.125 1.1 1.075 1.05 1.025 0.975 0.95 0.925 0.9 0.875 0.85 0.825 0.8 0.775 0.75 0.725 0.7 0.675 0.65 0.625</multiplier>
        </time_stepping></process></processes>
        <output><type>VTK</type><prefix>{prefix}</prefix><suffix>_ts_{{:timestep}}_t_{{:time}}</suffix>
            <fixed_output_times>1728000 3456000</fixed_output_times>
            <variables><variable>sigma</variable><variable>swelling_stress</variable><variable>saturation</variable>
                <variable>porosity</variable><variable>micro_pressure</variable><variable>micro_water_content</variable></variables>
        </output>
    </time_loop>
    <parameters>
        <parameter><name>sigma0</name><type>Function</type>
            <expression>-1.0e4</expression><expression>-1.0e4</expression><expression>-1.0e4</expression><expression>0</expression></parameter>
        <parameter><name>E</name><type>Constant</type><value>52e6</value></parameter>
        <parameter><name>nu</name><type>Constant</type><value>0.3</value></parameter>
        <parameter><name>phi0</name><type>Constant</type><value>{phi0}</value></parameter>
        <parameter><name>k0</name><type>Constant</type><value>{k0}</value></parameter>
        <parameter><name>displacement0</name><type>Constant</type><values>0 0</values></parameter>
        <parameter><name>pressure_ic</name><type>Constant</type><value>-{suction}</value></parameter>
        <parameter><name>zero</name><type>Constant</type><value>0.0</value></parameter>
        <parameter><name>p_bc_scale</name><type>Constant</type><value>1</value></parameter>
        <parameter><name>pressure_bc</name><type>CurveScaled</type><curve>suction_ramp</curve><parameter>p_bc_scale</parameter></parameter>
    </parameters>
    <process_variables>
        <process_variable><name>displacement</name><components>2</components><order>1</order>
            <initial_condition>displacement0</initial_condition><boundary_conditions>
            <boundary_condition><mesh>square_1e-2_quad_1e0_left</mesh><type>Dirichlet</type><component>0</component><parameter>zero</parameter></boundary_condition>
            <boundary_condition><mesh>square_1e-2_quad_1e0_right</mesh><type>Dirichlet</type><component>0</component><parameter>zero</parameter></boundary_condition>
            <boundary_condition><mesh>square_1e-2_quad_1e0_bottom</mesh><type>Dirichlet</type><component>1</component><parameter>zero</parameter></boundary_condition>
            <boundary_condition><mesh>square_1e-2_quad_1e0_top</mesh><type>Dirichlet</type><component>1</component><parameter>zero</parameter></boundary_condition>
        </boundary_conditions></process_variable>
        <process_variable><name>pressure</name><components>1</components><order>1</order>
            <initial_condition>pressure_ic</initial_condition><boundary_conditions>
            <boundary_condition><mesh>square_1e-2_quad_1e0_bottom</mesh><type>Dirichlet</type><parameter>pressure_bc</parameter></boundary_condition>
            <boundary_condition><mesh>square_1e-2_quad_1e0_top</mesh><type>Dirichlet</type><parameter>pressure_bc</parameter></boundary_condition>
        </boundary_conditions></process_variable>
    </process_variables>
    <nonlinear_solvers><nonlinear_solver><name>basic_newton</name><type>Newton</type><max_iter>60</max_iter><linear_solver>ls</linear_solver></nonlinear_solver></nonlinear_solvers>
    <linear_solvers><linear_solver><name>ls</name><eigen><solver_type>SparseLU</solver_type><scaling>true</scaling></eigen></linear_solver></linear_solvers>
    <curves><curve><name>suction_ramp</name><coords>0 1728000 3456000</coords><values>-{suction} 0 0</values></curve></curves>
</OpenGeoSysProject>
"""

def run_K(mat, K):
    m = MAT[mat]; prefix = f"cal_{mat}"
    outdir = f"{ROOT}/results/cal_{mat}"
    os.makedirs(outdir, exist_ok=True)
    prj = PRJ.format(mat=mat, rho_d=m["rho_d"], phi0=m["phi0"], n_s=m["n_s"],
                     k0=m["k0"], suction=m["suction"], target=m["target"], K=K, prefix=prefix)
    pf = f"{MODEL}/{prefix}.prj"
    open(pf, "w").write(prj)
    r = subprocess.run([OGS, f"{prefix}.prj", "-o", outdir], cwd=MODEL,
                       capture_output=True, text=True)
    fs = glob.glob(f"{outdir}/{prefix}_*.vtu")
    if not fs:
        return None, "no vtu (run failed)"
    fs.sort(key=lambda f: float(re.search(r"_t_([0-9.]+)\.vtu", f).group(1)))
    mesh = meshio.read(fs[-1])
    g = lambda n: float(np.asarray(mesh.point_data[n]).mean())
    sig = np.asarray(mesh.point_data["sigma"]).mean(axis=0)
    sig_ax = -sig[1] / 1e6  # MPa, compression positive
    S = g("saturation"); Pi = g("micro_pressure") / 1e6
    return sig_ax, f"S={S:.3f} Pi={Pi:.2f}MPa"

def calibrate(mat, tol=0.02, itmax=5):
    m = MAT[mat]; target = m["target"]; K = m["K0"]
    log = [f"\n===== CALIBRATE {mat} (target Ps={target} MPa, rho_d={m['rho_d']}) ====="]
    for it in range(itmax):
        sig, info = run_K(mat, K)
        if sig is None:
            log.append(f"  it{it}: K={K:.1f} -> FAILED ({info})"); break
        rel = (sig - target) / target
        log.append(f"  it{it}: K={K:9.1f} -> sigma_ax={sig:7.4f} MPa  ({rel*100:+5.1f}%)  [{info}]")
        if abs(rel) <= tol:
            log.append(f"  CONVERGED: K={K:.1f} J/kg -> Ps={sig:.4f} MPa (target {target}, within {tol*100:.0f}%)")
            return K, sig, log
        K = K * target / sig  # sigma ~ K (approx) -> Newton-ish scaling
    log.append(f"  STOPPED at K={K:.1f}, last sigma={sig}")
    return K, sig, log

if __name__ == "__main__":
    out = []
    for mat in (sys.argv[1:] or ["block", "pellet"]):
        Kf, sf, log = calibrate(mat)
        out += log
        out.append(f"  >>> {mat}: FITTED K = {Kf:.1f} J/kg  (Ps = {sf:.4f} MPa)")
    txt = "\n".join(out)
    open(f"{ROOT}/results/CALIBRATION_RESULT.txt", "w").write(txt + "\n")
    print(txt)
