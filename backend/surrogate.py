"""Etot surrogate model wrapper (an external MATLAB dFBA scan fitted to a
piecewise-cubic equation).

This whole package (limonene_predictor.py + limonene_model.json) is not
produced by this Python/COBRA tool -- it comes from a completely separate
pipeline that claims to run an exhaustive, adaptive Etot scan against the
same V9 model using real MATLAB + the COBRA Toolbox, and compresses the
conclusion "DXS/IDI/GPPS have no effect on the optimum within the given
range (backed by a structural certificate), and LS is the only effective
variable" into a single piecewise-cubic equation that only takes LS as
input. The upside is obvious: the website doesn't need to run COBRA live or
load a huge four-dimensional lookup table -- a single query is a pure
math operation, millisecond-scale.

These numbers couldn't originally be verified from the files alone (being
internally consistent -- hashes, unit tests, and a validation report all
lining up -- doesn't prove the numbers are real), so the live COBRA
computation paths (/api/fba, /api/dfba, /api/etot-scan) were deliberately
kept separate from, and not replaced by, this surrogate model. Since then,
this tool's own dFBA (now aligned with the external V9's _core_ biomass
reaction and medium setup) has been used to cross-check it step by step at
three test points -- LS = 3, 1000, 3000 uM -- including the early
infeasibility termination times, and every one matched to 0.000000%
relative error. That's an independent cross-validation, not just the
file's own self-consistency claim. Even so, this endpoint is still
deliberately presented separately from the live-COBRA paths (see the
`provenance` field below), because the surrogate model, after all, only
covers this one archived model and this one set of culture conditions --
it isn't a general conclusion that extrapolates without limit.
"""

import os

from .config import BASE_DIR

_etot_surrogate = None
_etot_surrogate_error = None

# Shared provenance marker for all three surrogate endpoints: the success
# path always reports "already cross-validated and reproduced by this
# tool"; a predictor load failure (missing file, bad format) reports a load
# error instead -- "was the load OK" and "have the numbers themselves been
# validated" are two different questions, and the strings shouldn't blur
# them together.
SURROGATE_PROVENANCE_OK = "external_matlab_scan_surrogate_now_reproduced_by_this_tool"
SURROGATE_PROVENANCE_LOAD_ERROR = "external_matlab_scan_surrogate_load_failed"


def get_etot_surrogate():
    global _etot_surrogate, _etot_surrogate_error
    if _etot_surrogate is None and _etot_surrogate_error is None:
        try:
            from limonene_predictor import LimonenePredictor
            model_path = os.environ.get(
                "LIMONENE_SURROGATE_MODEL", os.path.join(BASE_DIR, "limonene_model.json"))
            _etot_surrogate = LimonenePredictor.load(model_path)
        except Exception as exc:
            _etot_surrogate_error = str(exc)
    return _etot_surrogate, _etot_surrogate_error