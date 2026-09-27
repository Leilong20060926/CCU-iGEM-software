"""Flask application instance and all HTTP routes (formerly the "4. Flask
API" section of the single-file main.py).
"""

import math
import os

from cobra.flux_analysis import pfba
from flask import Flask, jsonify, request, send_from_directory, render_template

from .config import BASE_DIR, ASSETS_DIR, FRONTEND_TEMPLATES_DIR, FRONTEND_STATIC_DIR
from .dfba_engine import T7_ENZYME_RXNS, run_dfba, compute_enzyme_capacities, apply_t7_kinetic_rules
from .model_builder import apply_adjustable_medium, metabolite_gross_production_rate, build_medium_pathway_report
from .model_cache import get_registry, get_base_model, working_model, merge_params
from .optimize import optimize_etot_nelder_mead
from .scans import etot_scan, etot_scan_2d, etot_scan_multi, etot_scan_sampled, generate_qmc_samples
from .surrogate import get_etot_surrogate, SURROGATE_PROVENANCE_OK, SURROGATE_PROVENANCE_LOAD_ERROR

# static_folder only points at the assets/ subfolder (logo, day/night-toggle
# icons, and other static image files) -- not the whole BASE_DIR -- to
# avoid accidentally exposing files like main.py / LIM009_params.csv /
# models/ that shouldn't be public. When adding image files, remember to
# put them in the assets/ folder next to main.py.
# template_folder points at frontend/templates/, so index.html can
# {% include %} the page-section partials under templates/partials/ --
# Flask/Jinja2 stitches them together into one page at request time, so the
# browser still just sees one HTML document at "/".
app = Flask(__name__, static_folder=ASSETS_DIR, static_url_path="/assets",
            template_folder=FRONTEND_TEMPLATES_DIR)


@app.route("/api/etot-surrogate", methods=["POST"])
def etot_surrogate_endpoint():
    """Quick estimate from the external MATLAB scan surrogate model (see the
    module docstring in surrogate.py -- already cross-validated against
    this tool's own dFBA, matching to 0.000000% relative error at multiple
    test points). Still deliberately kept separate from, and not a
    replacement for, live-COBRA paths like /api/fba; the response always
    carries a `provenance` field identifying its source, and the frontend
    must not omit displaying it."""
    predictor, load_error = get_etot_surrogate()
    if predictor is None:
        return jsonify({
            "error": f"Failed to load the surrogate model: {load_error}",
            "provenance": SURROGATE_PROVENANCE_LOAD_ERROR,
        }), 503

    body = request.get_json(force=True, silent=True) or {}
    etot = body.get("etot") or {}
    try:
        values = {k: float(etot.get(k)) for k in ("DXS", "IDI", "GPPS", "LS")}
    except (TypeError, ValueError):
        return jsonify({"error": "Please supply all four DXS/IDI/GPPS/LS values (uM)"}), 400

    try:
        result = predictor.predict(values)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400

    result["provenance"] = SURROGATE_PROVENANCE_OK
    result["provenance_note"] = "Uses a pre-fitted equation for an instant answer, without waiting for a live COBRA computation."
    return jsonify(result)


@app.route("/api/etot-surrogate-curve", methods=["POST"])
def etot_surrogate_curve_endpoint():
    """The surrogate model is fundamentally a piecewise-cubic equation that
    takes LS Etot as input and outputs a limonene prediction (DXS/IDI/GPPS
    only affect the domain check, not the actual output value -- see
    equations.md); this computes the entire LS response curve in one go so
    the frontend can render it as a log-log chart. Since it's pure math and
    doesn't need to run COBRA, the sample-point count can be set quite high
    while staying fast (millisecond-scale); default 150 points."""
    predictor, load_error = get_etot_surrogate()
    if predictor is None:
        return jsonify({
            "error": f"Failed to load the surrogate model: {load_error}",
            "provenance": SURROGATE_PROVENANCE_LOAD_ERROR,
        }), 503

    body = request.get_json(force=True, silent=True) or {}
    etot = body.get("etot") or {}
    try:
        fixed = {k: float(etot.get(k)) for k in ("DXS", "IDI", "GPPS")}
    except (TypeError, ValueError):
        return jsonify({"error": "Please supply all three DXS/IDI/GPPS values (uM)"}), 400

    n_points = int(body.get("n_points", 150))
    n_points = max(20, min(n_points, 500))

    ls_lo, ls_hi = predictor.domain["LS"]
    log_lo, log_hi = math.log10(ls_lo), math.log10(ls_hi)
    ls_values = [
        10 ** (log_lo + (log_hi - log_lo) * i / (n_points - 1))
        for i in range(n_points)
    ]
    # Use the exact boundary values at the ends, so floating-point noise
    # can't push the last point just outside the domain and get rejected
    # by the predictor.
    ls_values[0], ls_values[-1] = ls_lo, ls_hi

    points = []
    for ls in ls_values:
        try:
            result = predictor.predict({**fixed, "LS": ls})
            points.append({
                "LS_uM": ls,
                "limonene_mM": result["limonene_mM"],
                "termination": result["termination"],
            })
        except ValueError:
            continue

    return jsonify({
        "points": points,
        "termination_brackets_uM": predictor.termination_brackets,
        "domain_uM": predictor.domain,
        "fixed_etot": fixed,
        "provenance": SURROGATE_PROVENANCE_OK,
    })


@app.route("/")
def index():
    # index.html is a Jinja2 template that {% include %}s the per-section
    # partials in templates/partials/ (header, params panel, results panel,
    # info page, about page) -- render_template stitches them into one page.
    return render_template("index.html")


@app.route("/styles.css")
def frontend_styles():
    return send_from_directory(FRONTEND_STATIC_DIR, "styles.css")


@app.route("/js/<path:filename>")
def frontend_scripts(filename):
    return send_from_directory(os.path.join(FRONTEND_STATIC_DIR, "js"), filename)


@app.errorhandler(FileNotFoundError)
def handle_missing_model(err):
    return jsonify({"error": str(err)}), 503


@app.route("/api/health")
def health():
    return jsonify({"status": "ok"})


@app.route("/api/defaults")
def defaults():
    r = get_registry()
    d = r.to_dfba_defaults()
    d["scan_values"] = {
        enzyme: r.scan_values(enzyme) for enzyme in ("DXS", "IDI", "GPPS", "LS")
    }
    return jsonify(d)


@app.route("/api/fba", methods=["POST"])
def fba():
    """Steady-state FBA: first solve for max growth rate, then maximize
    limonene secretion with min_growth_frac as the growth lower bound."""
    body = request.get_json(force=True, silent=True) or {}
    params = merge_params(body.get("params"))
    induced = bool(body.get("induced", True))

    model = working_model()
    apply_adjustable_medium(model, params)
    capacities = compute_enzyme_capacities(params)
    time_now = params["iptg_start_time_h"] if induced else 0
    apply_t7_kinetic_rules(
        model, capacities, time_now,
        params["iptg_start_time_h"], params["pre_induction_expr"],
        params["induced_expr"], params["induction_ramp_time_h"],
        params["burden_base_lb"], params["burden_max_lb"])

    biomass_rxn = params["biomass_rxn"]
    model.objective = biomass_rxn
    growth_sol = model.optimize()
    if growth_sol.status != "optimal":
        return jsonify({"status": growth_sol.status,
                         "message": "Growth optimization is infeasible; check the medium or Etot settings"})

    mu_max = growth_sol.fluxes[biomass_rxn]
    model.reactions.get_by_id(biomass_rxn).lower_bound = params["min_growth_frac"] * mu_max
    model.objective = "EX_limonene_e"
    # Mirrors LimoneneCOBRA001.m's optimizeCbModel(model,'max','one'): 'one'
    # means "among all solutions tied for optimal (a common FBA
    # degeneracy), pick the one with the smallest sum of flux magnitudes
    # (L1 norm)" -- which is exactly what pFBA does. Plain model.optimize()
    # has no such tie-breaking rule -- the solver just returns whichever
    # one it happens to find, which can make individual reaction fluxes
    # (and therefore ATP/NADPH cofactor totals) disagree with the reference
    # results even though the objective value (limonene yield) matches.
    # Using pFBA instead keeps this aligned.
    try:
        prod_sol = pfba(model)
    except Exception:
        # pFBA can be infeasible in rare edge cases (e.g. the target flux
        # is itself 0); fall back to plain FBA instead of failing the whole
        # request.
        prod_sol = model.optimize()

    atp_rate = None
    nadph_rate = None
    if prod_sol.status == "optimal":
        atp_rate = metabolite_gross_production_rate(model, prod_sol, "atp_c")
        nadph_rate = metabolite_gross_production_rate(model, prod_sol, "nadph_c")

    return jsonify({
        "status": prod_sol.status,
        "induced": induced,
        "mu_max": mu_max,
        "growth_floor": params["min_growth_frac"] * mu_max,
        "growth_rate": prod_sol.fluxes.get(biomass_rxn) if prod_sol.status == "optimal" else None,
        "limonene_flux": prod_sol.fluxes.get("EX_limonene_e") if prod_sol.status == "optimal" else None,
        "enzyme_fluxes": {
            name: prod_sol.fluxes.get(rxn_id)
            for name, rxn_id in T7_ENZYME_RXNS.items()
        } if prod_sol.status == "optimal" else {},
        "capacities_mmol_gDW_h": capacities,
        "cofactors": {
            "atp_production_mmol_gDW_h": atp_rate,
            "nadph_production_mmol_gDW_h": nadph_rate,
        },
        "medium_pathway_report": build_medium_pathway_report(model, capacities),
    })


@app.route("/api/dfba", methods=["POST"])
def dfba():
    """A single dynamic FBA simulation; returns the full time series."""
    body = request.get_json(force=True, silent=True) or {}
    params = merge_params(body.get("params"))
    model = working_model()
    apply_adjustable_medium(model, params)
    result = run_dfba(model, params)
    return jsonify(result)


@app.route("/api/etot-scan", methods=["POST"])
def etot_scan_endpoint():
    """One-at-a-time Etot scan on a single enzyme (mirrors LimoneneCOBRA009's
    OAT scan)."""
    body = request.get_json(force=True, silent=True) or {}
    enzyme = body.get("enzyme", "LS")
    if enzyme not in ("DXS", "IDI", "GPPS", "LS"):
        return jsonify({"error": "enzyme must be one of DXS, IDI, GPPS, or LS"}), 400

    params = merge_params(body.get("params"))
    scan_values = body.get("scan_values") or get_registry().scan_values(enzyme)

    model = get_base_model()
    results = etot_scan(model, params, enzyme, scan_values)
    return jsonify({"enzyme": enzyme, "results": results})


@app.route("/api/etot-scan-2d", methods=["POST"])
def etot_scan_2d_endpoint():
    """Cross (two-factor) Etot scan over two enzymes; returns a grid result.
    Kept around so existing callers aren't broken; the newer multi-enzyme
    scan frontend always calls /api/etot-scan-multi instead."""
    body = request.get_json(force=True, silent=True) or {}
    enzyme_x = body.get("enzyme_x", "LS")
    enzyme_y = body.get("enzyme_y", "DXS")
    valid_enzymes = ("DXS", "IDI", "GPPS", "LS")
    if enzyme_x not in valid_enzymes or enzyme_y not in valid_enzymes:
        return jsonify({"error": "enzyme_x / enzyme_y must be one of DXS, IDI, GPPS, or LS"}), 400
    if enzyme_x == enzyme_y:
        return jsonify({"error": "enzyme_x and enzyme_y must be different"}), 400

    params = merge_params(body.get("params"))
    values_x = body.get("values_x") or get_registry().scan_values(enzyme_x)
    values_y = body.get("values_y") or get_registry().scan_values(enzyme_y)

    model = get_base_model()
    result = etot_scan_2d(model, params, enzyme_x, enzyme_y, values_x, values_y)
    return jsonify(result)


@app.route("/api/etot-scan-multi", methods=["POST"])
def etot_scan_multi_endpoint():
    """Cross-combined scan over 2-4 enzymes at once (generalizes the
    fixed-two-enzyme 2D scan).
    body: { enzymes: ["DXS","GPPS",...] (2-4 of them), values_list: [[...],[...],...],
            params: {...} }"""
    body = request.get_json(force=True, silent=True) or {}
    enzymes = body.get("enzymes") or []
    values_list = body.get("values_list") or []
    valid_enzymes = ("DXS", "IDI", "GPPS", "LS")

    if not (2 <= len(enzymes) <= 4):
        return jsonify({"error": "enzymes must contain 2 to 4 entries"}), 400
    if len(enzymes) != len(set(enzymes)):
        return jsonify({"error": "enzymes must not repeat the same enzyme"}), 400
    if any(e not in valid_enzymes for e in enzymes):
        return jsonify({"error": "enzymes must be one of DXS, IDI, GPPS, or LS"}), 400
    if len(values_list) != len(enzymes):
        return jsonify({"error": "values_list must have the same length as enzymes"}), 400

    params = merge_params(body.get("params"))
    resolved_values = []
    for enzyme, values in zip(enzymes, values_list):
        resolved_values.append(values or get_registry().scan_values(enzyme))

    total_points = 1
    for v in resolved_values:
        total_points *= max(1, len(v))
    if total_points > 2000:
        return jsonify({
            "error": f"Total grid points {total_points} is too many (limit 2000); "
                     "reduce the number of enzymes or the number of values per enzyme"}), 400

    model = get_base_model()
    result = etot_scan_multi(model, params, enzymes, resolved_values)
    return jsonify(result)


@app.route("/api/etot-scan-sampled", methods=["POST"])
def etot_scan_sampled_endpoint():
    """Samples 2-4 enzymes via LHS or Sobol instead of an exhaustive grid --
    when a meaningful resolution would need a huge number of grid points
    (especially with 3-4 enzymes), a fixed number of sample points captures
    the overall trend much faster.
    body: { enzymes: [...] (2-4 of them), bounds: [[low,high],...], method: "lhs"|"sobol",
            n_samples: 50, params: {...} }"""
    body = request.get_json(force=True, silent=True) or {}
    enzymes = body.get("enzymes") or []
    bounds = body.get("bounds") or []
    method = body.get("method", "lhs")
    n_samples = int(body.get("n_samples", 50))
    valid_enzymes = ("DXS", "IDI", "GPPS", "LS")

    if not (2 <= len(enzymes) <= 4):
        return jsonify({"error": "enzymes must contain 2 to 4 entries"}), 400
    if len(enzymes) != len(set(enzymes)):
        return jsonify({"error": "enzymes must not repeat the same enzyme"}), 400
    if any(e not in valid_enzymes for e in enzymes):
        return jsonify({"error": "enzymes must be one of DXS, IDI, GPPS, or LS"}), 400
    if len(bounds) != len(enzymes):
        return jsonify({"error": "bounds must have the same length as enzymes"}), 400
    if any(len(b) != 2 or b[1] <= b[0] for b in bounds):
        return jsonify({"error": "each enzyme's bounds must be [low, high] with high > low"}), 400
    if method not in ("lhs", "sobol"):
        return jsonify({"error": "method must be lhs or sobol"}), 400
    n_samples = max(5, min(n_samples, 500))

    params = merge_params(body.get("params"))
    points = generate_qmc_samples(method, bounds, n_samples)

    model = get_base_model()
    results = etot_scan_sampled(model, params, enzymes, points)
    return jsonify({"enzymes": enzymes, "bounds": bounds, "method": method, "results": results})


@app.route("/api/optimize-etot", methods=["POST"])
def optimize_etot_endpoint():
    """Uses Nelder-Mead to automatically find the Etot(DXS/IDI/GPPS/LS)
    combination that maximizes limonene flux, in place of manually
    scanning an exhaustive grid. Starts the search from the current Etot
    field values."""
    body = request.get_json(force=True, silent=True) or {}
    params = merge_params(body.get("params"))
    # A simultaneous 4D search (DXS/IDI/GPPS/LS) needs 300-500 evaluations
    # in practice to converge reliably; around 60 is often cut off before
    # convergence -- the default was raised, and the allowed upper bound
    # widened to match.
    max_iter = int(body.get("max_iter", 200))
    max_iter = max(50, min(max_iter, 600))

    result = optimize_etot_nelder_mead(params, params["Etot"], max_iter=max_iter)
    return jsonify(result)