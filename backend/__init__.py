"""LimoSim backend package.

This package is the feature-organized successor to the original single-file
main.py. Each module below corresponds to one section of that file:

    config.py          Paths / environment configuration
    param_registry.py  LIM009_params.csv loader (formerly param_registry.py)
    model_builder.py    Model construction (formerly model_builder.py)
    dfba_engine.py       dFBA engine (formerly dfba_engine.py)
    optimize.py          Etot auto-optimization (Nelder-Mead simplex)
    fseof.py              FSEOF flux-scanning target search
    scans.py               Etot scan / 2D scan / multi-enzyme scan / QMC sampling
    worker_pool.py          Parallel computation infrastructure (ProcessPoolExecutor)
    surrogate.py             External MATLAB Etot surrogate model wrapper
    model_cache.py            Shared registry/model singletons + param merging
    app.py                     Flask app instance and all HTTP routes

`main.py`, at the project root, simply re-exports `app` from `backend.app`
so that both `python main.py` and `gunicorn main:app` keep working exactly
as before.
"""