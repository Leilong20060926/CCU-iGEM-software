"""FSEOF -- Flux Scanning based on Enforced Objective Flux
(Choi et al. 2010, Biotechnol Bioeng; Park et al. 2012).

A published algorithm for finding overexpression/downregulation gene
targets. Approach: progressively raise the enforced lower bound on a
target reaction's flux (e.g. limonene secretion); at each step,
re-maximize growth rate and record the whole-model flux distribution. As
the enforced flux climbs, reactions whose flux keeps moving in the same
direction are candidate targets: consistently rising -> overexpression
candidate, consistently falling -> downregulation candidate.
"""


def fseof_scan(model, biomass_rxn, target_rxn, n_steps=10, enforce_frac=0.9):
    """Runs an FSEOF scan and returns the whole-model flux distribution at
    each enforced-flux step, for later classification."""
    with model:
        model.objective = biomass_rxn
        growth_sol = model.optimize()
        if growth_sol.status != "optimal":
            raise ValueError("Growth optimization is infeasible; cannot run FSEOF")
        mu_max = growth_sol.fluxes[biomass_rxn]
        v_start = max(0.0, growth_sol.fluxes.get(target_rxn, 0.0))

    with model:
        # Keep a small growth lower bound so the scan doesn't wander into a
        # degenerate "zero growth, all product" solution.
        model.reactions.get_by_id(biomass_rxn).lower_bound = 0.1 * mu_max
        model.objective = target_rxn
        target_sol = model.optimize()
        if target_sol.status != "optimal":
            raise ValueError("Target-reaction optimization is infeasible; cannot run FSEOF")
        v_theoretical_max = target_sol.fluxes[target_rxn]

    v_end = enforce_frac * v_theoretical_max
    if v_end <= v_start or n_steps < 2:
        raise ValueError("Too little headroom to enforce on the target reaction; "
                          "cannot run an FSEOF scan (check that this reaction actually "
                          "has room to increase further under the current parameters)")

    enforced_values = [
        v_start + (v_end - v_start) * i / (n_steps - 1) for i in range(n_steps)
    ]

    flux_matrix = {}
    for level in enforced_values:
        with model:
            model.reactions.get_by_id(target_rxn).lower_bound = level
            model.objective = biomass_rxn
            sol = model.optimize()
            if sol.status != "optimal":
                for rxn in model.reactions:
                    flux_matrix.setdefault(rxn.id, []).append(None)
                continue
            for rxn in model.reactions:
                flux_matrix.setdefault(rxn.id, []).append(float(sol.fluxes.get(rxn.id, 0.0)))

    return {
        "enforced_values": enforced_values,
        "v_start": v_start,
        "v_theoretical_max": v_theoretical_max,
        "flux_matrix": flux_matrix,
    }


def classify_fseof_targets(fseof_result, biomass_rxn, target_rxn, min_flux_magnitude=1e-6):
    """Classifies each reaction from an FSEOF scan result as 'up'
    (overexpression candidate) or 'down' (downregulation candidate).
    Method: check whether the step-to-step flux changes are almost all in
    the same direction (by default requiring at least 80% of steps to
    agree), then rank by "how consistent the trend is" and "the total
    first-to-last flux change." Excludes the target reaction itself, the
    biomass reaction, and any reaction whose flux is essentially zero
    throughout (meaningless noise)."""
    flux_matrix = fseof_result["flux_matrix"]
    results = []
    for rxn_id, fluxes in flux_matrix.items():
        if rxn_id in (biomass_rxn, target_rxn):
            continue
        clean = [f for f in fluxes if f is not None]
        if len(clean) < 2:
            continue
        if max(abs(f) for f in clean) < min_flux_magnitude:
            continue

        diffs = [clean[i + 1] - clean[i] for i in range(len(clean) - 1)]
        n_up = sum(1 for d in diffs if d > 1e-9)
        n_down = sum(1 for d in diffs if d < -1e-9)
        n_total = len(diffs)
        if n_total == 0:
            continue

        consistency = max(n_up, n_down) / n_total
        if consistency < 0.8:
            continue

        direction = "up" if n_up >= n_down else "down"
        results.append({
            "reaction": rxn_id,
            "direction": direction,
            "consistency": consistency,
            "delta_flux": clean[-1] - clean[0],
            "start_flux": clean[0],
            "end_flux": clean[-1],
        })

    results.sort(key=lambda r: (-r["consistency"], -abs(r["delta_flux"])))
    return results