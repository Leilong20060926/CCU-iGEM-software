"""Etot auto-optimization (Nelder-Mead simplex, Nelder & Mead 1965).

Uses a handful of steady-state FBA solves instead of an exhaustive grid
scan, directly finding the Etot(DXS/IDI/GPPS/LS) combination that maximizes
limonene flux. Fully independent from the dFBA engine and doesn't affect
any previously-validated computation logic -- purely a new "find the
parameters for you" feature.
"""

from scipy.optimize import minimize

from .model_builder import apply_adjustable_medium
from .dfba_engine import compute_enzyme_capacities, apply_t7_kinetic_rules
from .model_cache import working_model


def _fba_limonene_flux_for_etot(etot_vec, base_params):
    """Given an Etot vector (DXS/IDI/GPPS/LS, uM), runs a single steady-state
    induced FBA and returns the limonene secretion flux. Any infeasible step
    returns 0 (so the optimizer treats that direction as a bad direction,
    rather than letting a single infeasible point abort the whole
    search)."""
    params = dict(base_params)
    params["Etot"] = {
        "DXS": max(0.0, etot_vec[0]), "IDI": max(0.0, etot_vec[1]),
        "GPPS": max(0.0, etot_vec[2]), "LS": max(0.0, etot_vec[3]),
    }
    model = working_model()
    apply_adjustable_medium(model, params)
    capacities = compute_enzyme_capacities(params)
    apply_t7_kinetic_rules(
        model, capacities, params["iptg_start_time_h"],
        params["iptg_start_time_h"], params["pre_induction_expr"],
        params["induced_expr"], params["induction_ramp_time_h"],
        params["burden_base_lb"], params["burden_max_lb"])

    biomass_rxn = params["biomass_rxn"]
    model.objective = biomass_rxn
    growth_sol = model.optimize()
    if growth_sol.status != "optimal":
        return 0.0
    mu_max = growth_sol.fluxes[biomass_rxn]
    model.reactions.get_by_id(biomass_rxn).lower_bound = params["min_growth_frac"] * mu_max
    model.objective = "EX_limonene_e"
    sol = model.optimize()
    if sol.status != "optimal":
        return 0.0
    return max(0.0, sol.fluxes.get("EX_limonene_e", 0.0))


def optimize_etot_nelder_mead(base_params, initial_etot, max_iter=200):
    """Uses Nelder-Mead to find the Etot combination that maximizes limonene
    flux. scipy's minimize only finds a minimum, so the internal objective
    function returns *negative* limonene flux, which is equivalent to
    maximizing it."""
    x0 = [initial_etot["DXS"], initial_etot["IDI"], initial_etot["GPPS"], initial_etot["LS"]]
    history = []

    def objective(x):
        flux = _fba_limonene_flux_for_etot(x, base_params)
        history.append({
            "etot": {
                "DXS": max(0.0, x[0]), "IDI": max(0.0, x[1]),
                "GPPS": max(0.0, x[2]), "LS": max(0.0, x[3]),
            },
            "limonene_flux": flux,
        })
        return -flux

    result = minimize(
        objective, x0, method="Nelder-Mead",
        options={"maxiter": max_iter, "xatol": 1e-3, "fatol": 1e-6, "adaptive": True})

    best_x = [max(0.0, v) for v in result.x]
    return {
        "best_etot": {"DXS": best_x[0], "IDI": best_x[1], "GPPS": best_x[2], "LS": best_x[3]},
        "best_limonene_flux": -result.fun,
        "n_evaluations": len(history),
        "converged": bool(result.success),
        "history": history,
    }