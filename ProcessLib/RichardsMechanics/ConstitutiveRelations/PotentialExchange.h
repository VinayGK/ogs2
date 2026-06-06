// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cmath>

#include "BaseLib/Error.h"

namespace ProcessLib::RichardsMechanics
{
struct YoungLaplaceMacroPotentialData
{
    double mu_LR = 0.0;
    double dmu_LR_dpLR = 0.0;
    double dmu_LR_drho_LR = 0.0;
    bool saturated_branch = true;
};

// DSM dsm_micromacro Phase-2 macro potential helper (Young-Laplace side):
// mu_LR = 0            for pLR > -ptol
// mu_LR = pLR / rho_LR otherwise
inline YoungLaplaceMacroPotentialData computeYoungLaplaceMacroPotential(
    double const p_LR, double const rho_LR, double const pressure_tolerance = 0.0)
{
    if (!(rho_LR > 0.0))
    {
        OGS_FATAL(
            "computeYoungLaplaceMacroPotential requires rho_LR > 0, got {:g}.",
            rho_LR);
    }
    if (!(pressure_tolerance >= 0.0))
    {
        OGS_FATAL(
            "computeYoungLaplaceMacroPotential requires pressure_tolerance >= "
            "0, got {:g}.",
            pressure_tolerance);
    }

    YoungLaplaceMacroPotentialData out;
    out.saturated_branch = (p_LR > -pressure_tolerance);
    if (out.saturated_branch)
    {
        return out;
    }

    out.mu_LR = p_LR / rho_LR;
    out.dmu_LR_dpLR = 1.0 / rho_LR;
    out.dmu_LR_drho_LR = -p_LR / (rho_LR * rho_LR);
    return out;
}

struct VanDerWaalsMicroPotentialData
{
    double omega_l = 0.0;
    double mu_lR = 0.0;

    double domega_l_dnl = 0.0;
    double domega_l_drho_lR = 0.0;
    double domega_l_dnS = 0.0;
    double domega_l_drho_SR = 0.0;

    double dmu_lR_dnl = 0.0;
    double dmu_lR_drho_lR = 0.0;   // NON-zero in the current form: computed as
                                   // -mu_lR/rho_lR (see line ~181, "/rho_lR fix"),
                                   // because mu_lR carries an explicit 1/rho_lR.
                                   // (Was exactly zero in the earlier reduced form.)
    double dmu_lR_dnS = 0.0;
    double dmu_lR_drho_SR = 0.0;
};

// DSM dsm_micromacro microscale vdW potential helper:
// omega_l = n_l * rho_lR / (nS * rho_SR)
// mu_lR_vdW = (A * Sa^3 / (6*pi)) * (nS^3 * rho_SR^3) / (n_l^3 * rho_lR)   [J/kg]
//
//   A      = hamaker_constant   [J]    Hamaker constant for clay-water-clay
//                                      Literature: montmorillonite-water-montmorillonite
//                                      A = 2.2e-20 J (Israelachvili & Adams 1978, SFA mica proxy)
//                                      Range: 1–5e-20 J (smectite, DLVO literature)
//                                      DO NOT calibrate A — it is a material constant.
//
// Dimensional derivation:
//   Film thickness h = n_l / (nS * rho_SR * Sa)                         [m]
//   Surface energy/area: E = -A/(12*pi*h^2)                             [J/m^2]
//   Specific free energy: mu_lR_vdW = E * (nS*rho_SR*Sa) / (nS*rho_lR)  [J/kg]
//   = -A*Sa^3*nS^3*rho_SR^3 / (12*pi*n_l^3*rho_lR) (adsorption sign)
//   Disjoining pressure Pi=-dE/dh gives factor 2: A*Sa^3*nS^3*rho_SR^3 / (6*pi*n_l^3*rho_lR)
//   Consistent with p_L_m = -rho_lR * mu_lR  [Pa]  (impl.h lines 276, 1044)
//
// Optional lumped exponential augmentation (activated when
// vdw_augmentation_prefactor > 0):
// h = n_l / (nS * rho_SR * Sa)               mean water film thickness [m]
// mu_lR_aug = K * exp(-h / lambda)
//   K      = vdw_augmentation_prefactor    [J/kg]  augmentation amplitude (calibrate to Villar)
//   lambda = vdw_augmentation_decay_length [m]     characteristic film thickness (calibrate)
// Total: mu_lR = sign * (mu_lR_vdW + mu_lR_aug)   — ADDITIVE, augmentation never replaces vdW
// Setting K = 0 (default) reduces exactly to the pure vdW form.
inline VanDerWaalsMicroPotentialData computeVanDerWaalsMicroPotential(
    double const n_l, double const rho_lR, double const nS, double const rho_SR,
    double const hamaker_constant, double const specific_surface,
    double const potential_sign_factor = 1.0,
    double const vdw_augmentation_prefactor = 0.0,
    double const vdw_augmentation_decay_length = 0.0)
{
    if (!(n_l > 0.0))
    {
        OGS_FATAL("computeVanDerWaalsMicroPotential requires n_l > 0, got {:g}.",
                  n_l);
    }
    if (!(rho_lR > 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires rho_lR > 0, got {:g}.",
            rho_lR);
    }
    if (!(nS > 0.0))
    {
        OGS_FATAL("computeVanDerWaalsMicroPotential requires nS > 0, got {:g}.",
                  nS);
    }
    if (!(rho_SR > 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires rho_SR > 0, got {:g}.",
            rho_SR);
    }
    if (!(hamaker_constant > 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires hamaker_constant > 0, "
            "got {:g}.",
            hamaker_constant);
    }
    if (!(specific_surface > 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires specific_surface > 0, "
            "got {:g}.",
            specific_surface);
    }
    if (!(vdw_augmentation_prefactor >= 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires "
            "vdw_augmentation_prefactor >= 0, got {:g}.",
            vdw_augmentation_prefactor);
    }
    if (vdw_augmentation_prefactor > 0.0 &&
        !(vdw_augmentation_decay_length > 0.0))
    {
        OGS_FATAL(
            "computeVanDerWaalsMicroPotential requires "
            "vdw_augmentation_decay_length > 0 when "
            "vdw_augmentation_prefactor > 0, got {:g}.",
            vdw_augmentation_decay_length);
    }

    constexpr double pi = 3.141592653589793238462643383279502884;

    VanDerWaalsMicroPotentialData out;

    out.omega_l = n_l * rho_lR / (nS * rho_SR);

    out.domega_l_dnl = rho_lR / (nS * rho_SR);
    out.domega_l_drho_lR = n_l / (nS * rho_SR);
    out.domega_l_dnS = -out.omega_l / nS;
    out.domega_l_drho_SR = -out.omega_l / rho_SR;

    double const prefactor = hamaker_constant * specific_surface *
                                 specific_surface * specific_surface /
                             (6.0 * pi);
    // Units (n_l, nS are dimensionless volume fractions, so [1]):
    //   A · Sa^3 · nS^3 · rho_SR^3 / (n_l^3 · rho_lR)
    //   [J] · [m^2/kg]^3 · [1] · [kg/m^3]^3 / ([1] · [kg/m^3])
    //   = [J · m^6/kg^3 · kg^3/m^9] / [kg/m^3]
    //   = [J/m^3] / [kg/m^3]
    //   = J/kg  ✓
    // The leading /rho_lR converts J/m^3 (Pa) to J/kg — dimensionally
    // required by the exchange equation rho_l_hat = alpha * (mu_LR - mu_lR).
    out.mu_lR = potential_sign_factor * prefactor * (nS * nS * nS) *
                (rho_SR * rho_SR * rho_SR) / (n_l * n_l * n_l * rho_lR);

    out.dmu_lR_dnl = -3.0 * out.mu_lR / n_l;
    out.dmu_lR_drho_lR = -out.mu_lR / rho_lR;  // non-zero after /rho_lR fix
    out.dmu_lR_dnS = 3.0 * out.mu_lR / nS;
    out.dmu_lR_drho_SR = 3.0 * out.mu_lR / rho_SR;

    // Lumped exponential force augmentation:
    // h = n_l / (nS * rho_SR * Sa)  [mean water film thickness, m]
    // mu_lR_aug = sign * K * exp(-h / lambda)
    if (vdw_augmentation_prefactor > 0.0)
    {
        double const xi = n_l / (vdw_augmentation_decay_length * nS *
                                 rho_SR * specific_surface);
        double const mu_aug =
            potential_sign_factor * vdw_augmentation_prefactor *
            std::exp(-xi);

        out.mu_lR += mu_aug;
        out.dmu_lR_dnl += -mu_aug * xi / n_l;
        // xi = n_l / (lambda * nS * rho_SR * Sa) — independent of rho_lR,
        // so mu_aug has no rho_lR dependence; dmu_lR_drho_lR from vdW term only
        out.dmu_lR_dnS += mu_aug * xi / nS;
        out.dmu_lR_drho_SR += mu_aug * xi / rho_SR;
    }

    return out;
}

// ── DSM Maxwell-conjugate term (B1) ─────────────────────────────────────────
// Restores the mean-effective-stress dependence of mu_lR — the Maxwell partner
// of the swelling eigenstress sigma_sw = -phi_m * Pi — so that (sigma, mu_lR)
// derive from ONE free energy Psi. Design + derivation:
//   ProcessLib/RichardsMechanics/DSM/MAXWELL_CONJUGATE_IMPLEMENTATION.md
//
//   mu_lR_mech = (1/rho_lR) * S1 * eps_v        for p_conf >= Pi   (sharp gate)
//              = 0                              otherwise
//   S1 = d sigma_sw,m / d n_l  (frozen phi, B1) [Pa] — supplied by the caller
//        from the eigenstress site (S1 = -n_S*(Pi + n_l*dPi_dnl)).
//
// Decisions (Vinay 2026-06-02): load EXPELS micro water; OGS effective stress,
// tension-positive; gate on confining pressure p_conf = -tr(sigma')/3 >= Pi;
// gate SHARP (the snap-drain at p_conf = Pi is the expected physical
// repercussion, not a numerical artifact); freeze phi when forming S1 (B1;
// B1.5 keep-phi-live deferred — see design doc ss.6.1).
struct MaxwellConjugateMicroPotentialData
{
    double mu_lR_mech = 0.0;           // J/kg
    double dmu_lR_mech_deps_v = 0.0;   // J/kg            per unit volumetric strain
    double dmu_lR_mech_dnl = 0.0;      // J/kg            per unit n_l
    double dmu_lR_mech_drho_lR = 0.0;  // (J/kg)/(kg/m^3)
    bool gate_open = false;
};

inline MaxwellConjugateMicroPotentialData computeMaxwellConjugateMicroPotential(
    double const S1, double const dS1_dnl, double const eps_v,
    double const p_conf, double const Pi, double const rho_lR,
    double const n_S)
{
    if (!(rho_lR > 0.0))
    {
        OGS_FATAL(
            "computeMaxwellConjugateMicroPotential requires rho_lR > 0, got "
            "{:g}.",
            rho_lR);
    }
    // n_S = 1 - phi_M = REV macro-solid (aggregate) fraction. mu_lR is a
    // SPECIFIC potential [J/kg]; its thermodynamic conjugate is the micro liquid
    // mass per REV volume m_l = rho_lR * phi_m = rho_lR * (1 - phi_M) * n_l, with
    // d m_l / d n_l = rho_lR * n_S. The eigenstress half is already REV-referenced
    // (sigma_sw = -(1 - phi_M)*n_l*Pi), so the conjugate half must divide by the
    // REV micro-liquid mass density rho_lR * n_S, NOT rho_lR alone. Omitting n_S
    // under-references the term by (1 - phi_M): the two cross-partials
    // d sigma_sw,m/d m_l and d mu_lR_mech/d eps_v then differ by (1 - phi_M),
    // the integrability condition d sigma/d m_l = d mu_lR/d eps fails, and the
    // (sigma, mu_lR) response is non-conservative (loop integral != 0). See
    // DSM/MAXWELL_CONJUGATE_REV_REFERENCING.md.
    if (!(n_S > 0.0))
    {
        OGS_FATAL(
            "computeMaxwellConjugateMicroPotential requires n_S = 1 - phi_M > 0, "
            "got {:g}.",
            n_S);
    }
    MaxwellConjugateMicroPotentialData out;
    // Sharp gate (B1): the term switches on only once the confining pressure
    // reaches the disjoining pressure. Below it the films out-push the wall, so
    // stress cannot expel water and the term is EXACTLY zero (the model stays
    // identical to the pre-Maxwell code there). This branch is a constitutive
    // choice (sharp Pi), not a numerical guard.
    out.gate_open = (p_conf >= Pi);
    if (!out.gate_open)
    {
        return out;
    }
    // REV-referenced measure: divide by rho_lR * n_S (per-REV micro-liquid mass
    // density), not rho_lR. This scales the conjugate UP by 1/(1 - phi_M), the
    // reciprocal of the (1 - phi_M) the eigenstress carries -> the pair derives
    // from one Psi (energy-conserving). phi_M is frozen at the GP (B1), so
    // d n_S / d rho_lR = 0 and the rho_lR-derivative keeps its form.
    double const rho_rev = rho_lR * n_S;                 // kg/m^3, per-REV measure
    out.mu_lR_mech = S1 * eps_v / rho_rev;               // J/kg
    out.dmu_lR_mech_deps_v = S1 / rho_rev;               // J/kg per unit strain
    out.dmu_lR_mech_dnl = dS1_dnl * eps_v / rho_rev;     // J/kg per unit n_l
    out.dmu_lR_mech_drho_lR = -out.mu_lR_mech / rho_lR;  // (J/kg)/(kg/m^3)
    return out;
}

// ── Film-pressure micro potential (maxwell beamer sec.5, consolidated form) ──
// SUPERSEDES the strain-view computeMaxwellConjugateMicroPotential above when
// the film_pressure_coupling flag is ON. The two halves (disjoining + the
// mechanical Maxwell partner) fuse into ONE potential of the *film pressure*
// p_film = Pi - b*p_conf:
//
//   mu_lR = -(Pi - b*p_conf) / rho_lR = -p_film / rho_lR,
//
// with Pi = -rho_lR*mu_lR_vdw > 0 (disjoining) and p_conf = -tr(sigma_eff)/3 > 0
// in compression (confining; OGS sigma_eff tension-positive). The vdW helper
// already returns mu_lR_vdw = -Pi/rho_lR, so this routine returns ONLY the FILM
// DELTA to ADD on top of mu_lR_vdw:
//
//   delta = mu_lR_film - mu_lR_vdw = +g * b * p_conf / rho_lR,
//
// where g is a C1 activation (smooth gate) of x = (p_conf - Pi_gate),
// Pi_gate = phi_m*Pi = n_S*n_l*Pi (option-1 REV-consistent threshold, decision
// #4 2026-06-02), of width w = film_pressure_gate_width. As p_conf rises past
// Pi_gate, mu_lR rises -> rho_l_hat = alpha*(mu_LR - mu_lR) turns negative ->
// the micro DRAINS (the smooth replacement of the sharp Macaulay snap).
//
//   w == 0  -> g = Heaviside(x)  (sharp; reproduces the old gate SHAPE exactly).
//   w  > 0  -> g = smoothstep cubic t^2*(3-2t), t = x/w  (C1: g' = 0 at both
//              ends), so the drain emerges continuously.
//
// Derivatives (the film delta D = g(x)*b*p_conf/rho_lR; p_conf and rho_lR are
// independent of n_l in this partial, x depends on n_l only through Pi_gate):
//   dD/dp_conf = g'(x)*b*p_conf/rho_lR + g*b/rho_lR
//   dD/dn_l    = g'(x)*(-dPi_gate/dn_l)*b*p_conf/rho_lR,
//                dPi_gate/dn_l = n_S*(Pi + n_l*dPi_dnl)
//   dD/drho_lR = -D/rho_lR
//
// Equipresence note: the (1 - phi_M) contact-area factor that the eigenstress
// carries CANCELS the per-REV-mass referencing of mu_lR (sec.5 slide "Contact
// area and density"), so the SPECIFIC potential mu_lR = -p_film/rho_lR carries
// NO explicit (1-phi_M) and divides by the intrinsic rho_lR — exactly mirroring
// the disjoining term mu_lR_vdw = -Pi/rho_lR it augments. n_S enters ONLY
// through the gate threshold Pi_gate = phi_m*Pi. (Contrast the strain-view
// helper above, whose conjugate divides by rho_lR*n_S because its driver S1 is
// the REV-scaled eigenstress slope.)
struct FilmPressureMicroPotentialData
{
    double mu_lR_film_delta = 0.0;        // J/kg            (add to mu_lR_vdw)
    double dmu_lR_film_dnl = 0.0;         // J/kg            per unit n_l
    double dmu_lR_film_drho_lR = 0.0;     // (J/kg)/(kg/m^3)
    double dmu_lR_film_dp_conf = 0.0;     // (J/kg)/Pa
    double gate_value = 0.0;              // g in [0, 1]
    bool gate_open = false;              // g > 0
};

inline FilmPressureMicroPotentialData computeFilmPressureMicroPotential(
    double const Pi, double const p_conf, double const rho_lR,
    double const n_S, double const n_l, double const dmu_lR_vdw_dnl,
    double const gate_width_w, double const biot_b)
{
    (void)dmu_lR_vdw_dnl;  // reserved (the film delta's n_l-dependence is via the
                           // gate threshold only; vdW's own dmu_lR_dnl is folded
                           // separately at the call site).
    if (!(rho_lR > 0.0))
    {
        OGS_FATAL(
            "computeFilmPressureMicroPotential requires rho_lR > 0, got {:g}.",
            rho_lR);
    }
    if (!(n_S > 0.0))
    {
        OGS_FATAL(
            "computeFilmPressureMicroPotential requires n_S = 1 - phi_M > 0, "
            "got {:g}.",
            n_S);
    }
    FilmPressureMicroPotentialData out;

    // dPi/dn_l from the vdW chain: Pi = -rho_lR*mu_lR_vdw, density-agnostic form
    // (matches the assembly-site convention dPi_dnl = (Pi/mu_lR)*dmu_lR_dnl).
    double const mu_lR_vdw =
        (rho_lR > 0.0) ? -Pi / rho_lR : 0.0;  // mu_lR_vdw = -Pi/rho_lR
    double const dPi_dnl =
        (std::abs(mu_lR_vdw) > 1e-300) ? (Pi / mu_lR_vdw) * dmu_lR_vdw_dnl : 0.0;

    // REV-consistent gate threshold Pi_gate = phi_m*Pi = n_S*n_l*Pi.
    double const Pi_gate = n_S * n_l * Pi;
    double const dPi_gate_dnl = n_S * (Pi + n_l * dPi_dnl);
    double const x = p_conf - Pi_gate;  // gate argument [Pa]

    // C1 activation g(x) and dg/dx.
    double g = 0.0;
    double dg_dx = 0.0;
    if (!(gate_width_w > 0.0))
    {
        // Sharp fallback (w == 0): Heaviside step. dg/dx = 0 a.e. (the delta at
        // x = 0 is not represented — the constitutive choice of a sharp gate).
        g = (x >= 0.0) ? 1.0 : 0.0;
        dg_dx = 0.0;
    }
    else if (x <= 0.0)
    {
        g = 0.0;
        dg_dx = 0.0;
    }
    else if (x >= gate_width_w)
    {
        g = 1.0;
        dg_dx = 0.0;
    }
    else
    {
        double const t = x / gate_width_w;          // in (0, 1)
        g = t * t * (3.0 - 2.0 * t);                // smoothstep, C1
        dg_dx = 6.0 * t * (1.0 - t) / gate_width_w;  // g'(x)
    }

    out.gate_value = g;
    out.gate_open = (g > 0.0);

    // Film delta D = g * b * p_conf / rho_lR  (= +b*p_conf/rho_lR when gate open
    // and saturated; turns mu_lR -Pi/rho_lR into -(Pi - b*p_conf)/rho_lR).
    double const D = g * biot_b * p_conf / rho_lR;
    out.mu_lR_film_delta = D;
    // dx/dp_conf = +1, dx/dn_l = -dPi_gate/dn_l.
    out.dmu_lR_film_dp_conf =
        dg_dx * biot_b * p_conf / rho_lR + g * biot_b / rho_lR;
    out.dmu_lR_film_dnl =
        dg_dx * (-dPi_gate_dnl) * biot_b * p_conf / rho_lR;
    out.dmu_lR_film_drho_lR = -D / rho_lR;
    return out;
}

struct PotentialDrivenMassExchangeData
{
    double rho_l_hat = 0.0;
    double rho_L_hat = 0.0;

    double drho_l_hat_dmu_LR = 0.0;
    double drho_l_hat_dmu_lR = 0.0;
    double drho_l_hat_dalpha_M = 0.0;
};

// DSM dsm_micromacro sign convention:
// rho_l_hat = alpha_M * (mu_LR - mu_lR)
// rho_L_hat = -rho_l_hat
inline PotentialDrivenMassExchangeData computePotentialDrivenMassExchange(
    double const alpha_M, double const mu_LR, double const mu_lR)
{
    PotentialDrivenMassExchangeData out;
    out.rho_l_hat = alpha_M * (mu_LR - mu_lR);
    out.rho_L_hat = -out.rho_l_hat;
    out.drho_l_hat_dmu_LR = alpha_M;
    out.drho_l_hat_dmu_lR = -alpha_M;
    out.drho_l_hat_dalpha_M = mu_LR - mu_lR;
    return out;
}
}  // namespace ProcessLib::RichardsMechanics
