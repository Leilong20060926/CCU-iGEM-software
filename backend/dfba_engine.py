"""dFBA engine (formerly the standalone dfba_engine.py).

Ports the IPTG-induction + Etot-capacity logic from
dynamicFBA_T7simple004.m.
"""

import math

from cobra.flux_analysis import pfba

T7_ENZYME_RXNS = {
    "DXS": "DXS_T7_het",
    "IDI_fwd": "IDI_T7_IPP_to_DMAPP_het",
    "IDI_rev": "IDI_T7_DMAPP_to_IPP_het",
    "GPPS": "GPPS_S80F_het",
    "LS": "LIMS_MS_het",
}

CAPACITY_ENZYME_FOR_RXN = {
    "DXS_T7_het": "DXS",
    "IDI_T7_IPP_to_DMAPP_het": "IDI",
    "IDI_T7_DMAPP_to_IPP_het": "IDI",
    "GPPS_S80F_het": "GPPS",
    "LIMS_MS_het": "LS",
}

# The kcat x Etot capacity formula applies to all four T7 enzymes
# (DXS/IDI/GPPS/LS), exactly matching kcatEtotToFluxV9() in
# LimoneneCOBRA008_Read_V9_capacity_trace_scan.m:
#   capacityFlux = 3.6 * kcat[s^-1] * Etot[uM] * cellVolume[L/gDW].
#
# GPPS/LS capacities used to be mistakenly hardcoded to a fixed
# 20.0 mmol/gDW/h (not varying with the Etot input). That came from
# misreading a reference report generated at "one particular baseline Etot
# value" (GPPS/LIMS_MS_het Upper Bound = 20.00) as if it were a fixed
# design cap -- it was actually just "the result computed at one baseline
# Etot value," not a constant. That bug made the GPPS/LS Etot scans always
# flat (capacity stuck at 20 regardless of the scanned value), which
# disagreed with the reference chart (the V9 Etot one-at-a-time capacity
# scan), whose LS curve spans several orders of magnitude. This special
# case has been removed; all four enzymes now uniformly use kcat x Etot.

BURDEN_RXN = "T7_BURDEN"

# Mirrors the exclUptakeRxns default in dynamicFBA_T7simple004.m: these
# four reactions (gases/solvents) are excluded from dynamic concentration
# tracking -- their bounds stay fixed at the medium-setup value and are
# never touched by the per-step mass-balance update. Treating them as an
# ordinary "untracked substrate" virtual pool could exhaust that pool over
# a long simulation, wrongly restricting oxygen/CO2 uptake -- this is a
# real behavioral difference from a naive dFBA implementation.
DYNAMIC_TRACKING_EXCLUDED_RXNS = {"EX_co2_e", "EX_o2_e", "EX_h2o_e", "EX_h_e"}


def etot_to_capacity_flux(kcat_s, etot_uM, cell_volume_L_per_gDW):
    """capacityFlux [mmol/gDW/h] = 3.6 * kcat[s^-1] * Etot[uM] * cellVolume[L/gDW]"""
    vmax_mM_h = 3.6 * kcat_s * etot_uM
    return vmax_mM_h * cell_volume_L_per_gDW


def compute_enzyme_capacities(params):
    """Computes the upper-bound capacity (mmol/gDW/h, before the induction
    fraction is applied) for all four T7 enzymes, using the same
    kcat x Etot formula uniformly for DXS/IDI/GPPS/LS."""
    return {
        enzyme: etot_to_capacity_flux(
            params["kcat"][enzyme], params["Etot"][enzyme],
            params["cell_volume_L_per_gDW"])
        for enzyme in ("DXS", "IDI", "GPPS", "LS")
    }


def apply_t7_kinetic_rules(model, capacities, time_now, iptg_start_time,
                            pre_induction_expr, induced_expr,
                            induction_ramp_time, burden_base_lb, burden_max_lb):
    """Applies the IPTG-induction expression fraction and the T7 enzyme
    capacity upper bounds (mutates `model` in place)."""
    if time_now < iptg_start_time:
        phase = "preIPTG_locked"
        expr_factor = pre_induction_expr
    else:
        phase = "postIPTG_induced"
        if induction_ramp_time > 0:
            ramp = min(1, max(0, (time_now - iptg_start_time) / induction_ramp_time))
            expr_factor = pre_induction_expr + ramp * (induced_expr - pre_induction_expr)
        else:
            expr_factor = induced_expr
    expr_factor = max(0, expr_factor)

    for rxn_id, enzyme in CAPACITY_ENZYME_FOR_RXN.items():
        if rxn_id not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(rxn_id)
        capacity = max(0, capacities.get(enzyme, 0))
        rxn.lower_bound = 0
        rxn.upper_bound = capacity * expr_factor

    burden_lb = burden_base_lb + (burden_max_lb - burden_base_lb) * expr_factor
    if BURDEN_RXN in model.reactions:
        burden = model.reactions.get_by_id(BURDEN_RXN)
        burden.upper_bound = max(burden.upper_bound, burden_lb)
        burden.lower_bound = burden_lb

    return expr_factor, phase


def solve_growth_only(model, biomass_rxn):
    model.objective = biomass_rxn
    # Mirrors dynamicFBA_T7simple004.m's
    # optimizeCbModel(modelStep,'max','one'): every dFBA step should be
    # solved with pFBA (the 'one' flag), not plain static FBA.
    try:
        sol = pfba(model)
    except Exception:
        sol = model.optimize()
    return sol, "growthOnly", None


def solve_production_priority(model, biomass_rxn, prod_rxn, min_growth_frac,
                               fallback_growth_fracs, growth_floor_reference,
                               reference_mu_max, min_positive_mu,
                               stop_on_non_positive_growth):
    """Maximizes production flux subject to a biomass lower bound, retrying
    with a progressively relaxed lower bound if needed."""
    growth_fracs = list(dict.fromkeys([min_growth_frac, *fallback_growth_fracs]))

    biomass_bounds = model.reactions.get_by_id(biomass_rxn).bounds
    model.objective = biomass_rxn
    try:
        growth_sol = pfba(model)
    except Exception:
        growth_sol = model.optimize()
    if growth_sol.status != "optimal":
        model.reactions.get_by_id(biomass_rxn).bounds = biomass_bounds
        return growth_sol, "prodPriority_noCurrentGrowthReference", None, float("nan")

    available_mu_max = max(0.0, growth_sol.fluxes[biomass_rxn])
    reference_mu = available_mu_max if growth_floor_reference == "currentStepMax" \
        else reference_mu_max

    sol = growth_sol
    mode = "prodPriority_failed"
    growth_floor = None
    for frac in growth_fracs:
        floor_now = max(0.0, frac * reference_mu)
        model.reactions.get_by_id(biomass_rxn).lower_bound = floor_now
        model.objective = prod_rxn
        try:
            sol = pfba(model)
        except Exception:
            sol = model.optimize()
        if sol.status == "optimal":
            mu = sol.fluxes[biomass_rxn]
            acceptable = (mu > min_positive_mu) if stop_on_non_positive_growth \
                else (mu >= -abs(min_positive_mu))
            if acceptable:
                mode = f"prodPriority_frac_{frac:g}"
                growth_floor = floor_now
                break
    model.reactions.get_by_id(biomass_rxn).bounds = biomass_bounds
    return sol, mode, growth_floor, available_mu_max


def update_uptake_bounds(model, exchange_ids, concentrations, biomass,
                          time_step, original_uptake_capacity, uptake_allowed):
    for rxn_id in exchange_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        if not uptake_allowed.get(rxn_id, False):
            rxn.lower_bound = 0
            continue
        conc = concentrations.get(rxn_id, 0.0)
        bound = 0.0 if biomass <= 0 else conc / (biomass * time_step)
        bound = min(bound, 1000, original_uptake_capacity.get(rxn_id, 0.0))
        if abs(bound) < 1e-12:
            bound = 0.0
        rxn.lower_bound = -bound


def run_dfba(model, params):
    """Runs a single dFBA simulation and returns the time series."""
    capacities = compute_enzyme_capacities(params)

    biomass_rxn = params["biomass_rxn"]
    prod_rxn = params.get("prod_rxn", "EX_limonene_e")
    substrate_rxn = params.get("substrate_rxn", "EX_glyc_e")
    dt = params["time_step_h"]
    n_steps = params["n_steps"]

    exchange_ids = [r.id for r in model.exchanges
                    if r.id not in DYNAMIC_TRACKING_EXCLUDED_RXNS]
    original_uptake_capacity = {r_id: max(0.0, -model.reactions.get_by_id(r_id).lower_bound)
                                 for r_id in exchange_ids}
    uptake_allowed = {r_id: original_uptake_capacity[r_id] > 0 for r_id in exchange_ids}

    concentrations = {r_id: 0.0 for r_id in exchange_ids}
    concentrations[substrate_rxn] = params["init_glycerol_mM"]
    # Phosphate (EX_pi_e) is, like scan_eval.m treats it, also an explicitly
    # tracked substrate with its own initial concentration -- it must not
    # fall into the "untracked -> default 1000 mM pool" branch below.
    phosphate_rxn = "EX_pi_e"
    if phosphate_rxn in exchange_ids:
        concentrations[phosphate_rxn] = params["init_phosphate_mM"]
    explicit_substrates = {substrate_rxn, phosphate_rxn}
    for r_id in exchange_ids:
        if uptake_allowed[r_id] and r_id not in explicit_substrates and concentrations[r_id] == 0:
            concentrations[r_id] = params["non_tracked_conc_mM"]

    biomass = params["init_biomass_gDW_L"]

    time_vec = [0.0]
    biomass_vec = [biomass]
    substrate_vec = [concentrations[substrate_rxn]]
    limonene_vec = [concentrations.get(prod_rxn, 0.0)]
    expr_vec = [None]
    growth_rate_vec = [None]
    limonene_flux_vec = [None]
    t7_flux_log = {name: [None] for name in T7_ENZYME_RXNS}
    # phase (preIPTG_locked / postIPTG_induced) and mode (growthOnly /
    # prodPriority_frac_X / growthFallback etc.) were already computed each
    # step, but used to be discarded rather than actually stored in the
    # returned result -- the reference MATLAB log (cached_benchmark.log)
    # already has both columns, so this fills that gap in, to make it
    # easier to diagnose "did this step get fallen back to growth-only"
    # type questions.
    phase_vec = [None]
    mode_vec = [None]

    update_uptake_bounds(model, exchange_ids, concentrations, biomass, dt,
                          original_uptake_capacity, uptake_allowed)

    reference_mu_max = float("nan")
    stopped_reason = None

    # step_mode lets the user choose how the time step advances:
    #   "fixed" (default): a fixed step size dt = time_step_h, for exactly
    #     n_steps steps -- this fully mirrors dynamicFBA_T7simple004.m's
    #     behavior and matches the MATLAB reference numbers.
    #   "adaptive": the step size is scaled dynamically based on how fast
    #     the growth rate is changing (shrinks automatically during rapid
    #     change, e.g. right as induction ramps up; grows automatically
    #     when things are flat). The resulting numbers will *not* exactly
    #     equal the fixed-step version -- that's an intentional trade-off
    #     (covering the same total simulated time in fewer steps, for
    #     speed, at the cost of being a different numerical approximation).
    #     Choosing this mode means you don't need bit-for-bit alignment
    #     with MATLAB.
    step_mode = params.get("step_mode", "fixed")
    target_total_time_h = n_steps * dt
    dt_min = dt / 10.0
    dt_max = dt * 4.0
    adaptive_target_rel_change = 0.08
    max_adaptive_steps = n_steps * 8
    prev_mu = None
    step = 0

    while True:
        if step_mode == "adaptive":
            if time_vec[-1] >= target_total_time_h - 1e-9 or step >= max_adaptive_steps:
                break
            dt = min(dt, target_total_time_h - time_vec[-1])
        else:
            if step >= n_steps:
                break

        t_now = time_vec[-1]

        expr_factor, phase = apply_t7_kinetic_rules(
            model, capacities, t_now,
            params["iptg_start_time_h"], params["pre_induction_expr"],
            params["induced_expr"], params["induction_ramp_time_h"],
            params["burden_base_lb"], params["burden_max_lb"])

        is_induced = t_now >= params["iptg_start_time_h"]
        available_mu_max = float("nan")

        if is_induced and params.get("use_prod_priority", True) and expr_factor > 0:
            sol, mode, growth_floor, available_mu_max = solve_production_priority(
                model, biomass_rxn, prod_rxn, params["min_growth_frac"],
                params["fallback_growth_fracs"], params["growth_floor_reference"],
                reference_mu_max, params["min_positive_mu"],
                params["stop_on_non_positive_growth"])
        else:
            sol, mode, growth_floor = solve_growth_only(model, biomass_rxn)

        if sol.status != "optimal" and is_induced and params.get("allow_growth_fallback", True):
            sol, mode, growth_floor = solve_growth_only(model, biomass_rxn)
            mode = "growthFallback"

        if sol.status != "optimal":
            stopped_reason = f"no feasible solution at step {step}, t={t_now:.2f}h"
            break

        mu = sol.fluxes[biomass_rxn]
        if not math.isfinite(available_mu_max):
            available_mu_max = mu
        reference_mu_max = available_mu_max if math.isfinite(available_mu_max) else reference_mu_max

        if mu <= params["min_positive_mu"]:
            if params["stop_on_non_positive_growth"]:
                stopped_reason = f"no positive growth at step {step}, t={t_now:.2f}h"
                break
            mu = max(mu, 0.0)

        uptake_flux = {r_id: sol.fluxes.get(r_id, 0.0) for r_id in exchange_ids}

        # Analytical mass-balance update (matches the MATLAB version:
        # integrates over exponential growth).
        if abs(mu) < 1e-12:
            delta_conc = {r_id: uptake_flux[r_id] * biomass * dt for r_id in exchange_ids}
            new_biomass = biomass
        else:
            growth_integral = biomass * (math.exp(mu * dt) - 1) / mu
            delta_conc = {r_id: uptake_flux[r_id] * growth_integral for r_id in exchange_ids}
            new_biomass = biomass * math.exp(mu * dt)

        new_conc = {}
        for r_id in exchange_ids:
            c = concentrations[r_id] + delta_conc[r_id]
            new_conc[r_id] = 0.0 if abs(c) < 1e-12 or c < 0 else c

        biomass = new_biomass
        concentrations = new_conc
        t_next = t_now + dt

        time_vec.append(t_next)
        biomass_vec.append(biomass)
        substrate_vec.append(concentrations[substrate_rxn])
        limonene_vec.append(concentrations.get(prod_rxn, 0.0))
        expr_vec.append(expr_factor)
        growth_rate_vec.append(mu)
        limonene_flux_vec.append(sol.fluxes.get(prod_rxn, 0.0))
        phase_vec.append(phase)
        mode_vec.append(mode)
        for name, rxn_id in T7_ENZYME_RXNS.items():
            t7_flux_log[name].append(sol.fluxes.get(rxn_id))

        update_uptake_bounds(model, exchange_ids, concentrations, biomass, dt,
                              original_uptake_capacity, uptake_allowed)

        if step_mode == "adaptive":
            # Scale the *next* step's size based on this step's relative
            # change in growth rate: bigger change -> smaller next step
            # (e.g. right as induction ramps up fully); flatter change ->
            # larger next step.
            if prev_mu is not None:
                rel_change = abs(mu - prev_mu) / max(abs(prev_mu), 1e-9)
                factor = math.sqrt(adaptive_target_rel_change / max(rel_change, 1e-9))
                factor = min(2.0, max(0.5, factor))
                dt = min(dt_max, max(dt_min, dt * factor))
            prev_mu = mu

        step += 1

    return {
        "time_h": time_vec,
        "biomass_gDW_L": biomass_vec,
        "glycerol_mM": substrate_vec,
        "limonene_mM": limonene_vec,
        "expr_factor": expr_vec,
        "growth_rate_h": growth_rate_vec,
        "limonene_flux": limonene_flux_vec,
        "phase": phase_vec,
        "mode": mode_vec,
        "t7_fluxes": t7_flux_log,
        "capacities_mmol_gDW_h": capacities,
        "stopped_reason": stopped_reason,
    }