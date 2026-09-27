"""Etot scans: one-at-a-time, 2D two-factor grid, N-factor grid, and
sampled (LHS / Sobol) multi-enzyme scanning.

The sampled approach differs from the exhaustive grid (etot_scan_multi):
instead of enumerating every combination of each enzyme's candidate
values, it draws a fixed-size, evenly-distributed set of sample points from
the whole bounds "volume" -- and the advantage grows with dimensionality.
Exhaustively scanning 4 enzymes at a meaningful resolution routinely needs
hundreds to thousands of grid points; sampling can capture the overall
trend with far fewer points (e.g. 50-100). Both methods are standard,
published techniques from the statistics literature, not something
invented for this project:
    LHS  : McKay, Conover & Beckman (1979), Technometrics 21(2), 239-245
    Sobol: Sobol' (1967), USSR Comput. Math. Math. Phys. 7(4), 86-112
Both are implemented directly via scipy.stats.qmc -- no custom algorithm
code was written for either.
"""

import itertools

from scipy.stats import qmc

from .worker_pool import run_params_grid


def etot_scan(model, base_params, enzyme, scan_values_uM):
    """Runs a one-at-a-time Etot scan on a single enzyme, holding every
    other enzyme at its baseline value. Each scanned value is independent
    of the others, so this is handed to run_params_grid for parallel
    computation.

    In addition to the usual final_limonene_mM/max_post_iptg_flux, this
    also derives three metrics from each result's full dFBA time series
    (trace) that mirror the V9 OAT capacity-scan diagnostic charts
    (mean post-induction flux, peak volumetric production rate, and the
    scanned enzyme's own peak capacity utilization). "Post-induction" is
    consistently defined as time_h >= iptg_start_time_h (matching how
    is_induced is determined elsewhere)."""
    params_list = []
    for value in scan_values_uM:
        p = dict(base_params)
        p["Etot"] = dict(base_params["Etot"])
        p["Etot"][enzyme] = value
        params_list.append(p)

    raw_results = run_params_grid(model, base_params, params_list)

    baseline_etot = base_params["Etot"].get(enzyme, 0.0)
    results = []
    for value, r, p in zip(scan_values_uM, raw_results, params_list):
        trace = r["trace"]
        iptg_t = p["iptg_start_time_h"]
        post_idx = [i for i, tv in enumerate(trace["time_h"]) if tv >= iptg_t]

        flux_series = trace["limonene_flux"]
        biomass_series = trace["biomass_gDW_L"]
        post_fluxes = [flux_series[i] for i in post_idx if flux_series[i] is not None]
        mean_post_flux = (sum(post_fluxes) / len(post_fluxes)) if post_fluxes else 0.0

        # Volumetric production rate (mM/h) = specific productivity
        # (mmol/gDW/h) x biomass (gDW/L). After unit conversion this comes
        # out as mmol/L/h = mM/h, reflecting how fast the whole culture is
        # actually accumulating product -- not the per-biomass-unit rate.
        vol_rates = [
            flux_series[i] * biomass_series[i]
            for i in post_idx if flux_series[i] is not None
        ]
        max_vol_rate = max(vol_rates) if vol_rates else 0.0

        capacity = trace.get("capacities_mmol_gDW_h", {}).get(enzyme, 0.0)
        t7 = trace.get("t7_fluxes", {})
        if enzyme == "IDI":
            # IDI is a reversible reaction split into two independent
            # directions, each with its own upper bound; either direction
            # may carry flux, so use whichever is larger as this enzyme's
            # current actual usage.
            fwd = t7.get("IDI_fwd", [])
            rev = t7.get("IDI_rev", [])
            enzyme_flux_series = [
                max(abs(fwd[i] or 0), abs(rev[i] or 0))
                for i in post_idx if i < len(fwd) and i < len(rev)
            ]
        else:
            series = t7.get(enzyme, [])
            enzyme_flux_series = [abs(series[i] or 0) for i in post_idx if i < len(series)]
        peak_utilization = (max(enzyme_flux_series) / capacity) \
            if enzyme_flux_series and capacity > 0 else 0.0

        results.append({
            "etot_uM": value,
            "etot_factor": (value / baseline_etot) if baseline_etot else None,
            "trace": trace,
            "final_limonene_mM": r["final_limonene_mM"],
            "max_post_iptg_flux": r["max_post_iptg_flux"],
            "mean_post_induction_flux_mmol_gDW_h": mean_post_flux,
            "max_volumetric_rate_mM_h": max_vol_rate,
            "peak_capacity_utilization": peak_utilization,
        })
    return results


def etot_scan_2d(model, base_params, enzyme_x, enzyme_y, values_x, values_y):
    """Two-factor grid scan: cross-combines the Etot grids for enzyme_x and
    enzyme_y, holding every other enzyme at its baseline value. Each grid
    point is independent, so it's handed to run_params_grid for parallel
    computation and then reshaped back into a 2D matrix in (y, x) order for
    the frontend's heatmap."""
    params_list = []
    for y_value in values_y:
        for x_value in values_x:
            p = dict(base_params)
            p["Etot"] = dict(base_params["Etot"])
            p["Etot"][enzyme_x] = x_value
            p["Etot"][enzyme_y] = y_value
            params_list.append(p)

    raw_results = run_params_grid(model, base_params, params_list)

    n_x = len(values_x)
    final_grid, flux_grid = [], []
    for row in range(len(values_y)):
        row_results = raw_results[row * n_x:(row + 1) * n_x]
        final_grid.append([r["final_limonene_mM"] for r in row_results])
        flux_grid.append([r["max_post_iptg_flux"] for r in row_results])

    return {
        "enzyme_x": enzyme_x,
        "enzyme_y": enzyme_y,
        "values_x": values_x,
        "values_y": values_y,
        "final_limonene_mM": final_grid,
        "max_post_iptg_flux": flux_grid,
    }


def etot_scan_multi(model, base_params, enzymes, values_list):
    """Cross-combined scan over 2-4 enzymes at once (a generalization of the
    old fixed-two-enzyme 2D scan). Returns a nested array whose dimension
    order matches `enzymes`: final_limonene_mM[i0][i1]...[iN], where i0
    indexes values_list[0], and so on. Each grid point is independent, so
    it's handed to run_params_grid for parallel computation -- 4 enzymes at
    5 values each is already 5^4=625 grid points, and the point count grows
    fast, so the frontend/endpoint rejects combinations that get too
    large."""
    grid_points = list(itertools.product(*values_list))
    params_list = []
    for point in grid_points:
        p = dict(base_params)
        p["Etot"] = dict(base_params["Etot"])
        for enzyme, value in zip(enzymes, point):
            p["Etot"][enzyme] = value
        params_list.append(p)

    raw_results = run_params_grid(model, base_params, params_list)
    final_flat = [r["final_limonene_mM"] for r in raw_results]
    flux_flat = [r["max_post_iptg_flux"] for r in raw_results]

    shape = [len(v) for v in values_list]

    def build_nested(flat, shape_remaining):
        if len(shape_remaining) == 1:
            return list(flat)
        size = shape_remaining[0]
        step = len(flat) // size
        return [build_nested(flat[i * step:(i + 1) * step], shape_remaining[1:])
                for i in range(size)]

    return {
        "enzymes": enzymes,
        "values_list": values_list,
        "final_limonene_mM": build_nested(final_flat, shape),
        "max_post_iptg_flux": build_nested(flux_flat, shape),
    }


def generate_qmc_samples(method, bounds, n_samples, seed=42):
    """method: "lhs" or "sobol". bounds: [(low, high), ...], one pair per
    enzyme. Returns n_samples sample points, each a float list of length
    len(bounds). The seed is fixed so re-running the same settings
    reproduces the same sample points (for reproducibility/debugging) --
    not for any security purpose."""
    d = len(bounds)
    if method == "lhs":
        sampler = qmc.LatinHypercube(d=d, seed=seed)
    elif method == "sobol":
        sampler = qmc.Sobol(d=d, seed=seed)
    else:
        raise ValueError(f"Unknown sampling method: {method} (must be lhs or sobol)")

    unit_samples = sampler.random(n=n_samples)
    lows = [b[0] for b in bounds]
    highs = [b[1] for b in bounds]
    scaled = qmc.scale(unit_samples, lows, highs)
    return scaled.tolist()


def etot_scan_sampled(model, base_params, enzymes, points):
    """points: each element is a set of Etot values of length len(enzymes)
    (drawn from LHS/Sobol sampling, not a regular grid, so it can't be
    reshaped into a matrix) -- returns a flat list where each record carries
    its own Etot combination plus the two result values."""
    params_list = []
    for point in points:
        p = dict(base_params)
        p["Etot"] = dict(base_params["Etot"])
        for enzyme, value in zip(enzymes, point):
            p["Etot"][enzyme] = value
        params_list.append(p)

    raw_results = run_params_grid(model, base_params, params_list)
    return [
        {
            "etot": dict(zip(enzymes, point)),
            "final_limonene_mM": r["final_limonene_mM"],
            "max_post_iptg_flux": r["max_post_iptg_flux"],
        }
        for point, r in zip(points, raw_results)
    ]