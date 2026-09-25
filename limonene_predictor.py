"""Small, dependency-free Etot predictor for the archived limonene dFBA model.

Layer 1 checks all four concentrations and classifies horizon/early infeasibility.
Layer 2 evaluates a local cubic in log10(LS), then reverses the output transform.
This is a surrogate of one particular model and fermentation setup, not a
biological claim that DXS, IDI, and GPPS can never affect limonene production.
"""

from __future__ import annotations

import argparse
from bisect import bisect_right
import json
import math
from pathlib import Path
from typing import Mapping


DEFAULT_MODEL = Path(__file__).with_name("limonene_model.json")


class LimonenePredictor:
    """Load once and reuse for fast individual or batch predictions."""

    def __init__(self, model: Mapping, *, allow_unvalidated: bool = False):
        if model.get("schema_version") != 3:
            raise ValueError("Unsupported model schema_version")
        if model.get("input_contract") != "primary_LP_verified_no_growth_fallback_v3":
            raise ValueError("Model lacks the independently verified primary-LP input contract")
        conditions = model.get("conditions", {})
        if conditions.get("IPTG_h") not in (2.75, 4.0) or conditions.get("horizon_h") != 25.0:
            raise ValueError("Invalid or missing model simulation conditions")
        if not allow_unvalidated and model.get("acceptance", {}).get("status") != "validated":
            raise ValueError("Model has not passed independent validation; use only as a diagnostic candidate")
        self.model = dict(model)
        self.domain = self.model["domain_uM"]
        self.termination_rules = self.model["termination_rules"]
        self.termination_brackets = self.model["termination_brackets_uM"]
        self.knots = self.model["log10_LS_knots"]
        if len(self.knots) < 2 or any(
            not math.isfinite(x) for x in self.knots
        ) or any(a >= b for a, b in zip(self.knots, self.knots[1:])):
            raise ValueError("Model knots must be finite and strictly increasing")
        for name, spec in self.model["outputs"].items():
            if spec["transform"] not in ("identity", "log10"):
                raise ValueError(f"Unsupported output transform for {name}")
            if len(spec["coefficients"]) != len(self.knots) - 1 or any(
                len(coef) != 4 or any(not math.isfinite(v) for v in coef)
                for coef in spec["coefficients"]
            ):
                raise ValueError(f"Invalid cubic coefficients for {name}")

    @classmethod
    def load(cls, path: str | Path = DEFAULT_MODEL, *, allow_unvalidated: bool = False) -> "LimonenePredictor":
        with Path(path).open(encoding="utf-8") as stream:
            return cls(json.load(stream), allow_unvalidated=allow_unvalidated)

    def predict(self, etot_uM: Mapping[str, float]) -> dict:
        """Return last-feasible cumulative product and mandatory termination metadata.

        Required keys are DXS, IDI, GPPS and LS. Inputs outside the certified
        model domain, missing keys, nonnumeric values, NaN and infinity raise
        ValueError. No extrapolation or silent clipping is performed.
        """
        if set(etot_uM) != set(self.domain):
            raise ValueError("Supply exactly these Etot keys: DXS, IDI, GPPS, LS")
        values = {}
        for enzyme, (low, high) in self.domain.items():
            value = etot_uM[enzyme]
            if isinstance(value, (bool, str, bytes)):
                raise ValueError(f"{enzyme} Etot must be a finite number in uM")
            try:
                value = float(value)
            except (TypeError, ValueError, OverflowError) as exc:
                raise ValueError(f"{enzyme} Etot must be a finite number in uM") from exc
            if not math.isfinite(value) or not low <= value <= high:
                raise ValueError(f"{enzyme} Etot must lie in [{low}, {high}] uM")
            values[enzyme] = value
        x = math.log10(values["LS"])
        index = min(max(bisect_right(self.knots, x) - 1, 0), len(self.knots) - 2)
        delta = x - self.knots[index]
        predictions = self._termination_metadata(values["LS"])
        for name, spec in self.model["outputs"].items():
            a, b, c, d = spec["coefficients"][index]
            value = ((a * delta + b) * delta + c) * delta + d
            if spec["transform"] == "log10":
                value = 10.0 ** value
            predictions[name] = value
        return predictions

    def _termination_metadata(self, ls_uM: float) -> dict:
        result = {
            "induction_h": self.model["conditions"]["IPTG_h"],
            "metric": "cumulative_limonene_until_min_25h_or_first_model_infeasibility",
            "output_validity": "surrogate_estimate_of_last_feasible_cumulative_product",
            "termination": "boundary_uncertain",
            "completed_25h": None,
            "endpoint_time_h": None,
            "observed_terminal_time_h_range": None,
            "boundary_bracket_uM": None,
            "classification_basis": "training_LS_regimes_with_explicit_unresolved_brackets",
        }
        for rule in self.termination_rules:
            low, high = rule["LS_range_uM"]
            if low <= ls_uM <= high:
                horizon = rule["termination"] == "horizon"
                result.update(termination="horizon" if horizon else "model_infeasible_before_25h",
                              completed_25h=horizon, endpoint_time_h=25.0 if horizon else None,
                              observed_terminal_time_h_range=rule["observed_terminal_time_h_range"])
                return result
        for low, high in self.termination_brackets:
            if low < ls_uM < high:
                result["boundary_bracket_uM"] = [low, high]
                return result
        raise ValueError("Model termination rules do not cover the requested LS")

    def interval(self, ls_uM: float) -> dict:
        """Explain the selected decision interval without computing a prediction."""
        if isinstance(ls_uM, bool):
            raise ValueError("LS must be a finite number in uM")
        low, high = self.domain["LS"]
        if not math.isfinite(ls_uM) or not low <= ls_uM <= high:
            raise ValueError(f"LS Etot must lie in [{low}, {high}] uM")
        index = min(max(bisect_right(self.knots, math.log10(ls_uM)) - 1, 0), len(self.knots) - 2)
        return {
            "index": index,
            "LS_lower_uM": 10.0 ** self.knots[index],
            "LS_upper_uM": 10.0 ** self.knots[index + 1],
            "upper_inclusive": index == len(self.knots) - 2,
        }


def predict(etot_uM: Mapping[str, float], model_path: str | Path = DEFAULT_MODEL) -> dict:
    """Convenience API; reuse LimonenePredictor for high-volume web requests."""
    return LimonenePredictor.load(model_path).predict(etot_uM)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", type=Path, default=DEFAULT_MODEL)
    for enzyme in ("dxs", "idi", "gpps", "ls"):
        parser.add_argument(f"--{enzyme}", type=float, required=True, help="Etot in uM")
    args = parser.parse_args()
    try:
        predictor = LimonenePredictor.load(args.model)
        values = {e.upper(): getattr(args, e) for e in ("dxs", "idi", "gpps", "ls")}
        result = {"predictions": predictor.predict(values), "interval": predictor.interval(args.ls)}
    except (ValueError, KeyError, OSError) as exc:
        parser.error(str(exc))
    print(json.dumps(result, ensure_ascii=False, allow_nan=False))


if __name__ == "__main__":
    main()