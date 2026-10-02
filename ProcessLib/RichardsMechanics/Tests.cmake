if (NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE RichardsMechanics/gravity.prj)
    OgsTest(PROJECTFILE RichardsMechanics/mechanics_linear.prj)
    OgsTest(PROJECTFILE RichardsMechanics/confined_compression_fully_saturated.prj RUNTIME 7)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_linear.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_anisotropic.prj)
    OgsTest(PROJECTFILE RichardsMechanics/flow_fully_saturated_coordinate_system.prj)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_small.prj RUNTIME 9)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_small_masslumping.prj RUNTIME 10)
    OgsTest(PROJECTFILE RichardsMechanics/RichardsFlow_2d_quasinewton.prj RUNTIME 80)
    OgsTest(PROJECTFILE RichardsMechanics/double_porosity_swelling.prj RUNTIME 20)
    OgsTest(PROJECTFILE RichardsMechanics/deformation_dependent_porosity.prj RUNTIME 8)
    OgsTest(PROJECTFILE RichardsMechanics/deformation_dependent_porosity_swelling.prj RUNTIME 11)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_power_law_permeability_xyz.prj RUNTIME 80)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_swelling_xyz.prj)
    OgsTest(PROJECTFILE RichardsMechanics/orthotropic_swelling_xy.prj)
    OgsTest(PROJECTFILE RichardsMechanics/bishops_effective_stress_power_law.prj)
    OgsTest(PROJECTFILE RichardsMechanics/bishops_effective_stress_saturation_cutoff.prj)
    OgsTest(PROJECTFILE RichardsMechanics/alternative_mass_balance_anzInterval_10.prj)
    if(NOT ENABLE_ASAN)
        OgsTest(PROJECTFILE RichardsMechanics/rotated_consolidation.prj RUNTIME 2)
    endif()
    OgsTest(PROJECTFILE RichardsMechanics/LiakopoulosHM/liakopoulos.prj RUNTIME 17)
    OgsTest(PROJECTFILE RichardsMechanics/LiakopoulosHM/liakopoulos_restart.xml RUNTIME 17)
    OgsTest(PROJECTFILE RichardsMechanics/LiakopoulosHM/liakopoulos_QN.prj RUNTIME 50)
    OgsTest(PROJECTFILE RichardsMechanics/A2.prj RUNTIME 20)
    OgsTest(PROJECTFILE RichardsMechanics/restart_w_backfill.prj RUNTIME 20)

    # ANCHORS EURAD-2 MS33 theoretical benchmarking — DSM native hierarchical runs
    # Re-timed 2026-08-31 on the dd900 K-knot adopt runs, OGS's own timer:
    # Model III 10.1542 s / 889 steps (III/run.log, OMP_NUM_THREADS=4) and
    # Model VII 207.022 s / 920 steps (VII/run_stdout.log, OMP_NUM_THREADS
    # unset -> 18-thread fallback), both under
    # /Users/vinaykumar/ogs-models/dd900_adopt_run_2026-08-31/, ogs
    # archive/dsm_native_Pi_fofnlev_branchtip_2026-08-11-49-gbed3e395. Both
    # stay inside their present RUNTIME, which is left unchanged; Model VII's
    # margin to RUNTIME 300 was measured at 18 threads, not at one.
    # Model I was NOT re-timed. Its figures — 0.111132 / 0.118636 / 0.149283 s
    # for dd1400 / dd1600 / dd1800, verbatim from the tracked
    # rerun_ms33_modelI_dd*.log — are 2026-05-22 runs of the now stale binary
    # vdw-baseline-2026-05-08-41-g9a1b956c, predating both commit 1bb414ac05
    # (log-linear live K, 2026-08-26) and the dd900 knot; those three decks
    # carry no <prefactors> table, so they were never part of the cascade. The
    # numbers are kept as the only timings on record for Model I, not as
    # current measurements, and RUNTIME 120 is left unchanged.
    # Lowering Model I/III RUNTIME below large_runtime = 60 would rename the
    # ctests (-LARGE suffix, scripts/cmake/test/OgsTest.cmake). Only Model IV
    # below was stale.
    # IC-FIX V0, 2026-09-30 (branch gen5_icfix_conformant_2026-09-30; Vinay 2026-09-30 "elevate the consistent IC to
    # submission", K_fix accepted): the decks below carry the K table refitted with the binary that realises the declared
    # initial_micro_water_content (commit 686fcd6ef8), and the references of Reference, III, IV and both VII are the fx outputs
    # of the readiness suite (Shilpa's MBP, bin/ogs md5 23640465bff1ef6ce5c83d6c5dcd48d3); the replaced frames are in
    # superseded_references_2026-09-30_icfix/ next to each. Final step counts: Reference 855 -> 895, III 1116 -> 970,
    # IV 16527 -> 16494, VII 240 d 1023 -> 1006, VII ladder550 1304 -> 1287 (0 rejected in all five; III had 5 rejected before).
    # Model I dd1400/1600/1800: references UNCHANGED (shipped-K frames); the K_fix re-run outputs pass all 11 checks of each
    # deck against them (vtkdiff replay, MEASURED on the mini), step counts 308/311/308 unchanged. RUNTIME values are unchanged;
    # measured wall times of the readiness runs (Shilpa's M1 Max, loaded, single-threaded assembly): Reference 20.3 s, III 349 s,
    # VII 369 s / 466 s, IV 16112.7 s (above RUNTIME 10736; the MBP run took 8957.6 s). NOT run through ctest (the builds used
    # for the fx runs have OGS_BUILD_TESTING=OFF): the ctest run on the MBP is the first open item of the V0 manifest.
    # MASS-CONSERVING DSM, CANDIDATE 2a, 2026-10-02 (branch gen5_masscons_1a_conformant_2026-10-02, cut from the v5 tip
    # 3101058538; Vinay 2026-10-02 "do both": candidate 2 cut twice, 2a = 1a, 2b = 1b; PLAN P2): the MS33 test tree (decks,
    # meshes, references, K_fix) is the one of gen5_icfix_conformant_2026-09-30 (d64c47cb30), and the nine committed decks
    # carry, inside <potential_exchange>, the KKT micro ceiling, F3, the ruled tangents and latched gate, the T_m drop,
    # variant 1a (darcy_relative_permeability_mobility = kirchhoff_element_mean) and macro_storage_exact_time_levels = true
    # (no every-step output, no diagnostic _ip variables). The other candidate-2 branch differs only in that variant line.
    # References: still the IC-fix frames listed above; whether they pass under this configuration is NOT tested yet.
    # <P7a: re-registered frames, step counts old -> new; replaced frames in superseded_references_2026-10-02_masscons_1a/>.
    # Model I: <P7a: replay result (D22)>. RUNTIME: <P7a/D21: measured walls>. ctest: <P7c: host, result>.
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj RUNTIME 120)
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj RUNTIME 120)
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1800.prj RUNTIME 120)
    # Model III ships the GAP-SWITCH deck (Vinay 2026-08-17). The outer radial
    # boundary swells free until u_r reaches the 2 mm technological gap, then
    # switches to a rigid Dirichlet wall — true container contact.
    # 2026-09-08 (spec-conformance fix 2, Vinay "approve all five, run them as one
    # campaign"): the deck now runs the clay r = 23 mm mesh ms33_clay_r23_h70
    # (920 quads) with the wall latching at r = 25 mm; measured 211.7 s / 943 steps
    # (OMP 2) and 2xx s at OMP 1 on 2026-09-08 -> RUNTIME 120 -> 300.
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIII/ms33_modelIII_gapswitch.prj RUNTIME 300)
    # Reference Configuration (confined dd1600 column, no gap), TRACKED and registered
    # 2026-09-08 (spec-conformance fix 4): 160x KC-base permeability + the live K(rho_d)
    # table of III/IV/VII, so the reference curve is on the models' hydraulics.
    # Measured 9.6 s / 771 steps (OMP 1) on 2026-09-08.
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_Reference/ms33_reference_dd1600.prj RUNTIME 60)
    # DEPRECATED 2026-08-17 — soft 2-medium gap annulus surrogate: no contact
    # mechanics, over-closes to ~67% with a residual aperture. Superseded by the
    # gap-switch deck above. Deck and reference retained (CLAUDE.md §6.2/§6.3);
    # registration commented out, not deleted.
    # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj RUNTIME 240)
    # RUNTIME covers the SLOWEST of three measured runs of this exact deck,
    # each 21500 accepted / 0 rejected steps to t = 200 d, OGS's own timer:
    #   10720    s  2026-08-28 09:02:27+0200, OMP_NUM_THREADS=4,
    #               ogs 6.5.8-565-gbea47887  (cascade_refit/IV_rep1/run.log)
    #   10735.4  s  2026-08-28 09:02:27+0200, OMP_NUM_THREADS=4,
    #               ogs 6.5.8-565-gbea47887  (cascade_refit/IV_rep2/run.log)
    #    8666.15 s  2026-08-31 09:23:34+0200, OMP_NUM_THREADS=6, ogs
    #            archive/dsm_native_Pi_fofnlev_branchtip_2026-08-11-49-gbed3e395
    #            (dd900_adopt_run_2026-08-31/IV/out/run.log, the run quoted in
    #            DSM/AGENTS.md entry 30)
    # The spread tracks the thread count, not the deck: the 6-thread run is the
    # fastest configuration on record, not the representative one. RUNTIME is
    # therefore taken from the slowest, 10736 = ceil(10735.4 s), which keeps the
    # repo's 2x convention against the worst case on record — RUNTIME > 750
    # makes OgsTest emit an explicit TIMEOUT = 2*RUNTIME = 21472 s, where
    # RUNTIME 8667 would give 17334 s, only 1.61x the slowest run. At the old
    # RUNTIME 240 no TIMEOUT is emitted at all, so ctest's default 1500 s would
    # apply and would kill a healthy run of this deck. No test is renamed by
    # the change: 240 and 10736 are both above large_runtime = 60, so this case
    # was already -LARGE. Effect on the OGS_CTEST_MAX_RUNTIME gate
    # (OgsTest.cmake) is band-dependent: a cap below 240 dropped the deck
    # before the change and still does — the only in-tree consumer,
    # scripts/ci/jobs/build-linux.yml, sets 60 — while a cap in [240, 10735]
    # kept the deck before and now drops it. Timings are
    # machine-local (this workstation, non-MPI, OpenMP thread counts as
    # listed); CI headroom is expected, not verified.
    # 2026-09-08 (spec-conformance fix 1, Vinay "approve all five, run them as one
    # campaign"): the deck now runs the CONCENTRIC mesh ms33_pellets_concentric_r25_h70
    # (clay core r <= 15 mm inside the pellet annulus) with the gen-4 900 knot retained
    # (the gen-5 knot stalls on this mesh at 27.28 d, see the deck); measured 7946 s /
    # 16008 steps at OMP 6 on 2026-09-08 -> RUNTIME 10736 kept (worst case on record).
    # 2026-09-08 campaign v2 (spec-audit 3x, LANE2 finding 3): output grid = the data-collection template's
    # 5-day rows (41 fixed output times); physics/BCs/knot unchanged; reference re-registered ts_16008 -> ts_16527
    # (old one in superseded_references_2026-09-08_v2/). Measured 7957.37 s / 16527 steps at OMP 6 on 2026-09-08
    # (campaign v2 runs/IV/run.log) -> RUNTIME 10736 kept.
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj RUNTIME 10736)
    # K(rho_d) equivalence pair (each material's k0 x20 spec, for speed): the
    # table-K variant resolves K = K(dry_density) at parse time and must
    # reproduce, bit-for-bit, the per-material scalar-K reference. Verified
    # 2026-06-08 (abs max diff = 0 on all 14 output fields at t=200 d). Both
    # registered run-only here.
    # DE-REGISTERED 2026-08-12 (Vinay): cannot pass as registered (no <test_definition>; OGS hard-fails at parse under the ctest wrapper). The two ModelIV variants additionally DIVERGE on the merged code (die ts #825 FD / #2333 analytic). Decks kept per never-delete; re-register only with ratified references.
    # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets_kref20x.prj RUNTIME 240)
    # DE-REGISTERED 2026-08-12 (Vinay): cannot pass as registered (no <test_definition>; OGS hard-fails at parse under the ctest wrapper). The two ModelIV variants additionally DIVERGE on the merged code (die ts #825 FD / #2333 analytic). Decks kept per never-delete; re-register only with ratified references.
    # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets_kofdd.prj RUNTIME 240)
    # 2026-09-08 (spec-conformance fix 3): step-and-hold traction ladder + 7 half-interval
    # output frames; measured 217.8 s / 973 steps (OMP 2) on 2026-09-08 -> RUNTIME 300 kept.
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj RUNTIME 300)
    # 2026-09-08 campaign v2 (spec-audit 3x, LANE1 F1 / LANE2 #1): LADDER-550 variant of Model VII on
    # the CURRENT spec ladder (theoretical_benchmarking.tex Table Loading_unloading_path: 50-d stages
    # 0.2/0.4/1/2.5/5/2.5/1/0.4 MPa, 10-d ramps, t_end 550 d; see the deck header for the citation).
    # SUPPLEMENTS the 240-d Q&A-schedule deck above (CLAUDE.md §3: supplement, never replace); which
    # ladder is the headline is pending Vinay's ruling. Own output prefix
    # ms33_modelVII_freeswelling_ladder550, so the two VII ctests never share an output file name.
    # Same two-tier tolerances as the parent deck. Measured 278.75 s / 1304 steps (OMP 2) on
    # 2026-09-08 (campaign v2 runs/VII_550/run.log) -> RUNTIME 300, as the parent.
    OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling_ladder550.prj RUNTIME 300)
    # K(rho_d) feature on a 2nd model (single-material Model VII -> table resolves
    # to the rho_d=1600 node, a physical no-op; k0 x50 spec for speed). Run to
    # t_end 2026-06-08. Exercises the table-resolution path on the free-swelling cell.
    # DE-REGISTERED 2026-08-12 (Vinay): cannot pass as registered (no <test_definition>; OGS hard-fails at parse under the ctest wrapper). The two ModelIV variants additionally DIVERGE on the merged code (die ts #825 FD / #2333 analytic). Decks kept per never-delete; re-register only with ratified references.
    # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling_kofdd.prj RUNTIME 300)
    # HISTORICAL (MS33 registrations of the cut point 3101058538 = v4 tip ca3c9faf00, superseded on this branch by the
    # gen-5 block above, 2026-10-02; kept, never deleted):
    #| # ANCHORS EURAD-2 MS33 theoretical benchmarking — DSM native hierarchical runs
    #| OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1400.prj RUNTIME 120)
    #| OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1600.prj RUNTIME 120)
    #| OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelI/ms33_modelI_dd1800.prj RUNTIME 120)
    #| # Model III ships the GAP-SWITCH deck (Vinay 2026-08-17). The outer radial
    #| # boundary swells free until u_r reaches the 2 mm technological gap, then
    #| # switches to a rigid Dirichlet wall — true container contact.
    #| OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIII/ms33_modelIII_gapswitch.prj RUNTIME 120)
    #| # DEPRECATED 2026-08-17 — soft 2-medium gap annulus surrogate: no contact
    #| # mechanics, over-closes to ~67% with a residual aperture. Superseded by the
    #| # gap-switch deck above. Deck and reference retained (CLAUDE.md §6.2/§6.3);
    #| # registration commented out, not deleted.
    #| # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIII/ms33_modelIII_gap2mm.prj RUNTIME 240)
    #| OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets.prj RUNTIME 240)
    #| # K(rho_d) equivalence pair (each material's k0 x20 spec, for speed): the
    #| # table-K variant resolves K = K(dry_density) at parse time and must
    #| # reproduce, bit-for-bit, the per-material scalar-K reference. Verified
    #| # 2026-06-08 (abs max diff = 0 on all 14 output fields at t=200 d). Both
    #| # registered run-only here.
    #| # DE-REGISTERED 2026-08-12 (Vinay): cannot pass as registered (no <test_definition>; OGS hard-fails at parse under the ctest wrapper). The two ModelIV variants additionally DIVERGE on the merged code (die ts #825 FD / #2333 analytic). Decks kept per never-delete; re-register only with ratified references.
    #| # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets_kref20x.prj RUNTIME 240)
    #| # DE-REGISTERED 2026-08-12 (Vinay): cannot pass as registered (no <test_definition>; OGS hard-fails at parse under the ctest wrapper). The two ModelIV variants additionally DIVERGE on the merged code (die ts #825 FD / #2333 analytic). Decks kept per never-delete; re-register only with ratified references.
    #| # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelIV/ms33_modelIV_pellets_kofdd.prj RUNTIME 240)
    #| OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling.prj RUNTIME 300)
    #| # K(rho_d) feature on a 2nd model (single-material Model VII -> table resolves
    #| # to the rho_d=1600 node, a physical no-op; k0 x50 spec for speed). Run to
    #| # t_end 2026-06-08. Exercises the table-resolution path on the free-swelling cell.
    #| # DE-REGISTERED 2026-08-12 (Vinay): cannot pass as registered (no <test_definition>; OGS hard-fails at parse under the ctest wrapper). The two ModelIV variants additionally DIVERGE on the merged code (die ts #825 FD / #2333 analytic). Decks kept per never-delete; re-register only with ratified references.
    #| # OgsTest(PROJECTFILE RichardsMechanics/ANCHORS_MS33_ModelVII/ms33_modelVII_freeswelling_kofdd.prj RUNTIME 300)
endif()

if (NOT OGS_USE_MPI AND OGS_USE_MFRONT)
    OgsTest(PROJECTFILE RichardsMechanics/mfront_restart_part1.prj RUNTIME 1)
    OgsTest(PROJECTFILE RichardsMechanics/mfront_restart_part2.xml RUNTIME 1)
    OgsTest(PROJECTFILE RichardsMechanics/DoubleStructureBenchmark/double_porosity_swelling_RM.prj RUNTIME 1)
endif()

AddTest(
    NAME RichardsMechanics_double_porosity_swelling_dsm_micromacro_constbc_reference
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS double_porosity_swelling_dsm_micromacro_constbc.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu pressure pressure 1e-16 1e-12
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu saturation saturation 4e-15 0
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu porosity porosity 1e-16 0
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu transport_porosity transport_porosity 1e-16 0
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu micro_pressure micro_pressure 1e-16 1e-12
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu micro_saturation micro_saturation 4e-15 0
    dsm_micromacro_constbc_reference_t_1000.000000.vtu double_porosity_swelling_dsm_micromacro_constbc_t_1000.000000.vtu swelling_stress swelling_stress 5e-14 0
)

AddTest(
    NAME RichardsMechanics_beacon_1a01_dsm_micromacro_smoke
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1a01_dsm_micromacro_smoke.prj
    WRAPPER time
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME RichardsMechanics_beacon_1a01_dsm_micromacro_stressprobe
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1a01_dsm_micromacro_stressprobe.prj
    WRAPPER time
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME RichardsMechanics_beacon_1a01_dsm_micromacro_inflow
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1a01_dsm_micromacro_inflow.prj
    WRAPPER time
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME RichardsMechanics_beacon_1a01_dsm_micromacro_reference
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1a01_dsm_micromacro_smoke.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu displacement displacement 1e-12 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu pressure pressure 1e-12 1e-12
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu saturation saturation 1e-12 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu micro_pressure micro_pressure 1e-12 1e-12
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu micro_saturation micro_saturation 1e-12 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu swelling_stress swelling_stress 1e-12 0
    beacon_1a01_reference_t_1000.000000.vtu beacon_1a01_dsm_micromacro_smoke_t_1000.000000.vtu sigma sigma 1e-12 1e-10
)

AddTest(
    NAME RichardsMechanics_beacon_1a01_dsm_micromacro_inflow_reference
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1a01_dsm_micromacro_inflow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu displacement displacement 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu pressure pressure 1e-12 1e-12
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu saturation saturation 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_pressure micro_pressure 1e-12 1e-12
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_saturation micro_saturation 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_water_content micro_water_content 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_porosity micro_porosity 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu micro_exchange_source micro_exchange_source 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu swelling_stress swelling_stress 1e-12 0
    beacon_1a01_dsm_micromacro_inflow_reference_t_100000.000000.vtu beacon_1a01_dsm_micromacro_inflow_t_100000.000000.vtu sigma sigma 1e-12 1e-10
)

AddTest(
    NAME RichardsMechanics_beacon_1b_dsm_micromacro_smoke
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1b_dsm_micromacro_smoke.prj
    WRAPPER time
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME RichardsMechanics_beacon_1b_dsm_micromacro_reference
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1b_dsm_micromacro_smoke.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu displacement displacement 1e-12 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu pressure pressure 1e-12 1e-12
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu saturation saturation 1e-12 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu micro_pressure micro_pressure 1e-12 1e-12
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu micro_saturation micro_saturation 1e-12 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu swelling_stress swelling_stress 1e-12 0
    beacon_1b_reference_t_1000.000000.vtu beacon_1b_dsm_micromacro_smoke_t_1000.000000.vtu sigma sigma 1e-12 1e-10
)

AddTest(
    NAME RichardsMechanics_beacon_1c_dsm_micromacro_smoke
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1c_dsm_micromacro_smoke.prj
    WRAPPER time
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME RichardsMechanics_beacon_1c_dsm_micromacro_reference
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 20
    EXECUTABLE_ARGS beacon_1c_dsm_micromacro_smoke.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu displacement displacement 1e-12 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu pressure pressure 1e-12 1e-12
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu saturation saturation 1e-12 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu porosity porosity 1e-12 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu transport_porosity transport_porosity 1e-12 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu micro_pressure micro_pressure 1e-12 1e-12
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu micro_saturation micro_saturation 1e-12 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu swelling_stress swelling_stress 1e-12 0
    beacon_1c_reference_t_1000.000000.vtu beacon_1c_dsm_micromacro_smoke_t_1000.000000.vtu sigma sigma 1e-12 1e-10
)

AddTest(
    NAME RichardsMechanics_square_1e2_confined_compression_restart
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 8
    EXECUTABLE_ARGS confined_compression_fully_saturated_restart.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    # Does not exist?
    # PROPERTIES DEPENDS ogs-RichardsMechanics_square_1e2_confined_compression-time-vtkdiff
    DIFF_DATA
    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu displacement displacement 1e-16 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu displacement displacement 1e-16 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu displacement displacement 1e-16 0

    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu pressure pressure 1e-16 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu pressure pressure 1e-16 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu pressure pressure 1e-16 0

    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu sigma sigma 5e-14 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu sigma sigma 5e-14 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu sigma sigma 5e-14 0

    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu epsilon epsilon 5e-14 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu epsilon epsilon 5e-14 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu epsilon epsilon 5e-14 0

    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu saturation saturation 4e-15 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu saturation saturation 4e-15 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu saturation saturation 4e-15 0

    confined_compression_fully_saturated_ts_20_t_100.000000.vtu confined_compression_fully_saturated_restart_ts_0_t_100.000000.vtu velocity velocity 1e-16 0
    confined_compression_fully_saturated_ts_120_t_1000.000000.vtu confined_compression_fully_saturated_restart_ts_100_t_1000.000000.vtu velocity velocity 1e-16 0
    confined_compression_fully_saturated_ts_420_t_4000.000000.vtu confined_compression_fully_saturated_restart_ts_400_t_4000.000000.vtu velocity velocity 1e-16 0
)

AddTest(
    NAME RichardsMechanics_A2_total_initial_stress
    PATH RichardsMechanics
    EXECUTABLE ogs
    RUNTIME 15
    EXECUTABLE_ARGS A2_total_stress0.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    A2_ts_3_t_4320.000000.vtu A2_total_stess0_test_ts_3_t_4320.000000.vtu displacement displacement 1e-16 0
    A2_ts_42_t_20736.000000.vtu A2_total_stess0_test_ts_42_t_20736.000000.vtu displacement displacement 1e-16 0
    A2_ts_76_t_2764800.000000.vtu A2_total_stess0_test_ts_76_t_2764800.000000.vtu displacement displacement 1e-16 0

    A2_ts_3_t_4320.000000.vtu A2_total_stess0_test_ts_3_t_4320.000000.vtu pressure pressure 1e-16 1e-12
    A2_ts_42_t_20736.000000.vtu A2_total_stess0_test_ts_42_t_20736.000000.vtu pressure pressure 1e-16 1e-12
    A2_ts_76_t_2764800.000000.vtu A2_total_stess0_test_ts_76_t_2764800.000000.vtu pressure pressure 1e-16 1e-12

    A2_ts_3_t_4320.000000.vtu A2_total_stess0_test_ts_3_t_4320.000000.vtu sigma sigma 5e-8 0
    A2_ts_42_t_20736.000000.vtu A2_total_stess0_test_ts_42_t_20736.000000.vtu sigma sigma 5e-8 0
    A2_ts_76_t_2764800.000000.vtu A2_total_stess0_test_ts_76_t_2764800.000000.vtu sigma sigma 5e-8 0

    A2_ts_3_t_4320.000000.vtu A2_total_stess0_test_ts_3_t_4320.000000.vtu epsilon epsilon 5e-14 0
    A2_ts_42_t_20736.000000.vtu A2_total_stess0_test_ts_42_t_20736.000000.vtu epsilon epsilon 5e-14 0
    A2_ts_76_t_2764800.000000.vtu A2_total_stess0_test_ts_76_t_2764800.000000.vtu epsilon epsilon 5e-14 0

    A2_ts_3_t_4320.000000.vtu A2_total_stess0_test_ts_3_t_4320.000000.vtu saturation saturation 4e-15 0
    A2_ts_42_t_20736.000000.vtu A2_total_stess0_test_ts_42_t_20736.000000.vtu saturation saturation 4e-15 0
    A2_ts_76_t_2764800.000000.vtu A2_total_stess0_test_ts_76_t_2764800.000000.vtu saturation saturation 4e-15 0

    A2_ts_3_t_4320.000000.vtu A2_total_stess0_test_ts_3_t_4320.000000.vtu velocity velocity 1e-16 0
    A2_ts_42_t_20736.000000.vtu A2_total_stess0_test_ts_42_t_20736.000000.vtu velocity velocity 1e-16 0
    A2_ts_76_t_2764800.000000.vtu A2_total_stess0_test_ts_76_t_2764800.000000.vtu velocity velocity 1e-16 0
)
