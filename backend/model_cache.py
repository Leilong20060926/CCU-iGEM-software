"""Shared registry/base-model singletons, plus request-parameter merging.

Kept as its own module (rather than folded into app.py) because both
optimize.py and worker_pool.py need get_registry()/get_base_model() without
creating a circular import with app.py.
"""

import os

from .config import MODEL_PATH, PARAM_CSV
from .model_builder import build_model, apply_medium009, detect_biomass_rxn
from .param_registry import ParamRegistry

_registry = None
_base_model = None


def get_registry():
    global _registry
    if _registry is None:
        _registry = ParamRegistry(PARAM_CSV)
    return _registry


def get_base_model():
    global _base_model
    if _base_model is None:
        if not os.path.exists(MODEL_PATH):
            raise FileNotFoundError(
                f"Strain model file not found: {MODEL_PATH}\n"
                "Please place iEC1356_Bl21DE3.json (recommended) or .xml / "
                ".mat into the models/ folder, or set the "
                "LIMONENE_MODEL_PATH environment variable to point at the "
                "right path.")
        model = build_model(MODEL_PATH)
        model = apply_medium009(model, get_registry().to_dfba_defaults())
        model._biomass_rxn = detect_biomass_rxn(model)
        _base_model = model
    return _base_model


def working_model():
    return get_base_model().copy()


def merge_params(overrides):
    params = get_registry().to_dfba_defaults()
    if overrides:
        for key, value in overrides.items():
            if key == "Etot" and isinstance(value, dict):
                params["Etot"] = {**params["Etot"], **value}
            else:
                params[key] = value
    params["biomass_rxn"] = get_base_model()._biomass_rxn
    return params