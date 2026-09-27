"""Paths and environment configuration shared by the rest of the backend.

BASE_DIR is the project root (one level up from this backend/ package),
which is where LIM009_params.csv, index.html, models/, and the optional
surrogate-model files (limonene_predictor.py, limonene_model.json) live.
"""

import os


def _default_model_path(base_dir):
    """Look in models/ for .json > .xml > .mat, in that priority order, and
    return the first one that exists."""
    models_dir = os.path.join(base_dir, "models")
    for ext in (".json", ".xml", ".mat"):
        candidate = os.path.join(models_dir, f"iEC1356_Bl21DE3{ext}")
        if os.path.exists(candidate):
            return candidate
    # If none exist, still return the .json path so the resulting error
    # message tells the user exactly which file to add.
    return os.path.join(models_dir, "iEC1356_Bl21DE3.json")


# This file lives at <project_root>/backend/config.py, so the project root
# is two directories up from here.
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_PATH = os.environ.get("LIMONENE_MODEL_PATH", _default_model_path(BASE_DIR))
PARAM_CSV = os.environ.get(
    "LIMONENE_PARAM_CSV", os.path.join(BASE_DIR, "LIM009_params.csv"))
ASSETS_DIR = os.path.join(BASE_DIR, "assets")

# The frontend lives in its own folder, split into a Jinja2 template
# (index.html + templates/partials/*.html, one file per page section) and
# static assets (styles.css + js/*.js), instead of being a single HTML file
# at the project root.
FRONTEND_DIR = os.path.join(BASE_DIR, "frontend")
FRONTEND_TEMPLATES_DIR = os.path.join(FRONTEND_DIR, "templates")
FRONTEND_STATIC_DIR = os.path.join(FRONTEND_DIR, "static")