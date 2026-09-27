"""
main.py
--------
LimoSim FBA / dFBA web backend -- thin entry point.

The implementation lives in the backend/ package, split one file per
feature (see backend/__init__.py for the map). This file just re-exports
the Flask `app` object so both of these keep working exactly as before:

    python main.py
    gunicorn main:app

Startup:
    pip install -r requirements.txt
    python main.py
    -> http://localhost:5001

Required files (same directory as main.py):
    LIM009_params.csv                  V9 parameter registry table
    index.html                         Frontend page
    models/iEC1356_Bl21DE3.mat         COBRA Toolbox strain model (add this yourself)

Optional (missing files don't affect anything else; only the "Etot quick
estimate" tab will return 503):
    limonene_predictor.py              Query script for the external surrogate model
                                        (standard library only, no extra dependencies)
    limonene_model.json                Coefficients and metadata for the external surrogate model
    These two files come from an external MATLAB scan pipeline that is
    completely independent of this tool; this tool has not independently
    verified whether the numbers in them really came from an actual MATLAB
    run -- spot-checking them yourself before deployment is recommended.

Speed:
    pip install highspy
    -> Once installed, automatically switches to the HiGHS LP solver
       (many times faster than cobrapy's default GLPK), helping FBA / dFBA /
       scans / 2D scans alike, and is enabled by default. If it can't be
       installed, this silently falls back to GLPK without breaking
       anything.
       If you need bit-for-bit alignment with MATLAB/GLPK results, set the
       environment variable LIMONENE_FAST_SOLVER=0 to disable it and force
       GLPK.
    Separately, each grid cell of an Etot scan or two-enzyme 2D scan is
    independent, so they're automatically parallelized across multiple CPU
    cores in the background (via ProcessPoolExecutor) -- no extra setup
    needed, and this doesn't affect the solved results.
"""

from backend.app import app

if __name__ == "__main__":
    app.run(port=5001, debug=True)