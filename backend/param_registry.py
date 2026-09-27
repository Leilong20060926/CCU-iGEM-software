"""LIM009_params.csv loader (formerly the standalone param_registry.py)."""

import csv


def _to_bool(s):
    return str(s).strip().lower() in ("true", "1", "yes")


def _to_float_list(s):
    return [float(x) for x in str(s).split(";") if x.strip() != ""]


class ParamRegistry:
    """Loads LIM009_params.csv (the V9 parameter registry table)."""

    def __init__(self, csv_path):
        self._rows = {}
        with open(csv_path, newline="", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                self._rows[row["parameter_code"]] = row

    def scalar(self, code):
        return float(self._rows[code]["v9_value"])

    def scalar_or_default(self, code, default):
        """Like scalar(), but returns `default` instead of raising when the
        CSV doesn't have this parameter_code. Used for tunable parameters
        that aren't (yet) registered in LIM009_params.csv (e.g. the ATPM /
        DXPS / METAT / FPPS bounds)."""
        row = self._rows.get(code)
        if row is None:
            return default
        try:
            return float(row["v9_value"])
        except (KeyError, ValueError, TypeError):
            return default

    def string(self, code):
        return self._rows[code]["v9_value"]

    def boolean(self, code):
        return _to_bool(self._rows[code]["v9_value"])

    def vector(self, code):
        return _to_float_list(self._rows[code]["v9_value"])

    def scan_values(self, enzyme):
        """The scan_values field (uM) on the Etot.<enzyme> row."""
        row = self._rows[f"Etot.{enzyme}"]
        return _to_float_list(row["scan_values"])

    def to_dfba_defaults(self):
        r = self
        return {
            "kcat": {
                "DXS": r.scalar("kcat.DXS"),
                "IDI": r.scalar("kcat.IDI"),
                "GPPS": r.scalar("kcat.GPPS"),
                "LS": r.scalar("kcat.LS"),
            },
            "Etot": {
                "DXS": r.scalar("Etot.DXS"),
                "IDI": r.scalar("Etot.IDI"),
                "GPPS": r.scalar("Etot.GPPS"),
                "LS": r.scalar("Etot.LS"),
            },
            "cell_volume_L_per_gDW": r.scalar("conv.cellVolume_L_per_gDW"),
            "iptg_start_time_h": r.scalar("sim.iptgStartTime_h"),
            "pre_induction_expr": r.scalar("sim.preInductionExpr"),
            "induced_expr": r.scalar("sim.inducedExpr"),
            "induction_ramp_time_h": r.scalar("sim.inductionRampTime_h"),
            "init_glycerol_mM": r.scalar("sim.initGlycerol_mM"),
            # Mirrors scan_eval.m's dynamicFBA_cached(...,{'EX_glyc_e','EX_pi_e'},
            # [54.8;89],...): phosphate (EX_pi_e) is, just like glycerol, an
            # explicitly tracked substrate with its own initial concentration,
            # not an infinite untracked pool. sim.initPhosphate_mM has always
            # been present in the CSV (its own comment says "explicit initial
            # concentration tracked by dFBA"), but it was previously never
            # actually read/applied -- so phosphate was being treated as
            # "untracked" by our code and defaulted to a 1000 mM virtual pool,
            # 11x larger than the real 89 mM. That made phosphate effectively
            # never run out, which didn't match the real MATLAB trajectory
            # (phosphate is depleted at t~3.75h, and growth rate drops to 0
            # right after). Reading it here fixes that.
            "init_phosphate_mM": r.scalar_or_default("sim.initPhosphate_mM", 89.0),
            "init_biomass_gDW_L": r.scalar("sim.initBiomass_gDW_L"),
            "time_step_h": r.scalar("sim.timeStep_h"),
            "n_steps": int(r.scalar("sim.nSteps")),
            "min_growth_frac": r.scalar("sim.minGrowthFrac"),
            "fallback_growth_fracs": r.vector("sim.fallbackGrowthFracs"),
            "min_positive_mu": r.scalar("sim.minPositiveMu"),
            "growth_floor_reference": r.string("sim.growthFloorReference"),
            "burden_base_lb": r.scalar("sim.burdenBaseLB"),
            "burden_max_lb": r.scalar("sim.burdenMaxLB"),
            "non_tracked_conc_mM": r.scalar("sim.nonTrackedConc_mM"),
            "use_prod_priority": r.boolean("sim.useProdPriority"),
            "allow_growth_fallback": r.boolean(
                "sim.allowGrowthFallbackAfterProductionFailure"),
            "stop_on_non_positive_growth": r.boolean(
                "sim.stopOnNonPositiveGrowth"),
            "maxUptake.EX_glyc_e": r.scalar("maxUptake.EX_glyc_e"),
            "maxUptake.EX_o2_e": r.scalar("maxUptake.EX_o2_e"),
            "sim.richAminoAcidMaxUptake": r.scalar("sim.richAminoAcidMaxUptake"),
            # The following are native-pathway bounds the frontend can adjust.
            # If LIM009_params.csv doesn't have the matching parameter_code
            # yet, fall back to a default instead of raising.
            # Mirrors LimoneneCOBRA008_Read_V9.m: that script never touches
            # ATPM at all, meaning V9 just uses the model file's own native
            # lower bound (3.15 mmol/gDW/h in this iEC1356_Bl21DE3 model
            # file). The 4.0 that used to be here was copied from
            # LimoneneCOBRA001.m's assumption; V9 has no such override, and
            # removing it dropped our error vs. the external V9 MATLAB
            # results from 0.87% to 0.0000% (confirmed at multiple test
            # points).
            "atpm_lb": r.scalar_or_default("sim.atpmLB", 3.15),
            "dxps_ub": r.scalar_or_default("bounds.DXPS_ub", 30.0),
            "metat_ub": r.scalar_or_default("bounds.METAT_ub", 0.10),
            # Note: the original LimoneneCOBRA001.m script sets this to 0.05
            # (a knock-down, not a full block), which disagrees with the
            # 0.00 in a reference report you gave earlier. This follows the
            # script's 0.05 for now -- if 0.00 is actually correct, let us
            # know and it can be changed back.
            "fpps_ub": r.scalar_or_default("bounds.FPPS_ub", 0.05),
            # dFBA time-step mode: "fixed" (default, matches the MATLAB
            # reference results) or "adaptive" (adaptive step size --
            # trades exact bit-for-bit alignment with MATLAB for speed).
            "step_mode": "fixed",
        }