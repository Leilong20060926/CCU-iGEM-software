"""Parallel-computation infrastructure.

Every grid cell of an Etot scan / 2D scan is an independent dFBA
simulation, so we run several of them at once across multiple processes --
roughly as many times faster as the machine has cores available.
"""

import concurrent.futures
import os

from .config import MODEL_PATH
from .model_builder import build_model, apply_medium009, detect_biomass_rxn, apply_adjustable_medium
from .dfba_engine import run_dfba
from .model_cache import get_registry

_worker_model = None  # Each worker process's own model object (not shared across processes).

# os.cpu_count() in some containerized environments (e.g. Render's
# resource-limited free tier) reports the *host's* core count, not what
# this specific container was actually allotted. Sizing the worker-process
# pool off that raw core count means every process loads its own full copy
# of the genome-scale model, multiplying memory usage by the worker count --
# that's an easy way to get OOM-killed on a small 512MB dyno, especially
# for the heaviest feature, multi-enzyme scanning. So this can be overridden
# via the LIMONENE_MAX_POOL_WORKERS environment variable; when unset, it
# falls back to the original os.cpu_count()-1 behavior (unaffected for
# local development, or on machines with genuinely enough memory).
# On memory-constrained deployments, set this explicitly to 1-2.
_env_max_workers = os.environ.get("LIMONENE_MAX_POOL_WORKERS")
if _env_max_workers:
    try:
        _MAX_POOL_WORKERS = max(1, int(_env_max_workers))
    except ValueError:
        _MAX_POOL_WORKERS = max(1, (os.cpu_count() or 2) - 1)
else:
    _MAX_POOL_WORKERS = max(1, (os.cpu_count() or 2) - 1)


def _init_pool_worker(model_path, medium_defaults):
    """ProcessPoolExecutor initializer: runs exactly once when each worker
    process starts. Loads the model and applies the medium setup here, so
    every task that worker later receives reuses the same model object
    (wrapped in `with model:` for temporary per-task changes). That way the
    parallel overhead is "read the model file once at startup," not "read
    the model file again for every grid cell" (which would be slower, not
    faster)."""
    global _worker_model
    model = build_model(model_path)
    model = apply_medium009(model, medium_defaults)
    model._biomass_rxn = detect_biomass_rxn(model)
    _worker_model = model


def _pool_dfba_task(params):
    """Runs one complete dFBA parameter set inside a worker process, and
    returns that set's summary result (including the full time series, for
    the frontend's trace charts)."""
    with _worker_model:
        apply_adjustable_medium(_worker_model, params)
        trace = run_dfba(_worker_model, params)
    fluxes = [f for f in trace["limonene_flux"] if f is not None]
    return {
        "trace": trace,
        "final_limonene_mM": trace["limonene_mM"][-1],
        "max_post_iptg_flux": max(fluxes) if fluxes else 0.0,
    }


def run_params_grid(model, base_params, params_list):
    """Runs several independent dFBA parameter sets in parallel.

    With only one parameter set, or when the machine doesn't appear to have
    multiple cores, runs sequentially in the current process instead --
    avoids the extra overhead of spinning up a process pool (loading the
    model file, starting processes) for cases where that overhead would
    make it slower than just running sequentially."""
    if len(params_list) <= 1 or _MAX_POOL_WORKERS <= 1:
        results = []
        with model:
            apply_adjustable_medium(model, base_params)
            for params in params_list:
                with model:
                    trace = run_dfba(model, params)
                fluxes = [f for f in trace["limonene_flux"] if f is not None]
                results.append({
                    "trace": trace,
                    "final_limonene_mM": trace["limonene_mM"][-1],
                    "max_post_iptg_flux": max(fluxes) if fluxes else 0.0,
                })
        return results

    n_workers = min(_MAX_POOL_WORKERS, len(params_list))
    medium_defaults = get_registry().to_dfba_defaults()
    with concurrent.futures.ProcessPoolExecutor(
            max_workers=n_workers,
            initializer=_init_pool_worker,
            initargs=(MODEL_PATH, medium_defaults)) as executor:
        return list(executor.map(_pool_dfba_task, params_list))