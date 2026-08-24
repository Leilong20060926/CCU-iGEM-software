"""
main.py
--------
LimoneneCOBRA FBA / dFBA 網頁後端（單檔版）。
整合了原本的 param_registry.py + model_builder.py + dfba_engine.py + app.py。

啟動方式：
    pip install -r requirements.txt
    python main.py
    → http://localhost:5001

必要檔案（放在 main.py 同一層目錄）：
    LIM009_params.csv                  V9 參數登錄表
    index.html                         前端網頁
    models/iEC1356_Bl21DE3.mat         COBRA Toolbox 菌株模型（需自行放入）
"""

import csv
import math
import os

import cobra
from cobra import Reaction, Metabolite
from flask import Flask, jsonify, request, send_from_directory


# ============================================================
# 1. 參數表載入器（原 param_registry.py）
# ============================================================

def _to_bool(s):
    return str(s).strip().lower() in ("true", "1", "yes")


def _to_float_list(s):
    return [float(x) for x in str(s).split(";") if x.strip() != ""]


class ParamRegistry:
    """載入 LIM009_params.csv (V9 參數登錄表)。"""

    def __init__(self, csv_path):
        self._rows = {}
        with open(csv_path, newline="", encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                self._rows[row["parameter_code"]] = row

    def scalar(self, code):
        return float(self._rows[code]["v9_value"])

    def string(self, code):
        return self._rows[code]["v9_value"]

    def boolean(self, code):
        return _to_bool(self._rows[code]["v9_value"])

    def vector(self, code):
        return _to_float_list(self._rows[code]["v9_value"])

    def scan_values(self, enzyme):
        """對應 Etot.<enzyme> 那一列的 scan_values 欄位（uM）。"""
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
        }


# ============================================================
# 2. 模型建構（原 model_builder.py）
#    移植 LimoneneCOBRA009.m 的 prepareModel009()
# ============================================================

HET_REACTIONS = {
    # rxn_id: (reactants, products, gene_rule)
    "DXS_T7_het": ("g3p_c + h_c + pyr_c", "co2_c + dxyl5p_c", "pET_T7_dxs"),
    "IDI_T7_IPP_to_DMAPP_het": ("ipdp_c", "dmpp_c", "pET_T7_idi"),
    "IDI_T7_DMAPP_to_IPP_het": ("dmpp_c", "ipdp_c", "pET_T7_idi"),
    "GPPS_S80F_het": ("dmpp_c + ipdp_c", "grdp_c + ppi_c", "pET_T7_ispA_S80F"),
    "LIMS_MS_het": ("grdp_c", "limonene_c + ppi_c", "pET_T7_msLS"),
    "PAIDS1_GPP_het": ("dmpp_c + ipdp_c", "grdp_c + ppi_c", "pET_T7_PaIDS1"),
    "PAIDS1_FPP_het": ("grdp_c + ipdp_c", "frdp_c + ppi_c", "pET_T7_PaIDS1"),
    "PAIDS1_GGPP_het": ("frdp_c + ipdp_c", "ggdp_c + ppi_c", "pET_T7_PaIDS1"),
    "T7_BURDEN": ("atp_c + h2o_c", "adp_c + pi_c + h_c", None),
}

ALL_T7_RXNS = list(HET_REACTIONS.keys())

TRACE_ION_RXNS = [
    "EX_fe2_e", "EX_fe3_e", "EX_mn2_e", "EX_zn2_e", "EX_cu2_e",
    "EX_cobalt2_e", "EX_mobd_e", "EX_ni2_e", "EX_sel_e", "EX_tungs_e",
]

RICH_MEDIUM_OPEN = [
    "EX_pi_e", "EX_so4_e", "EX_mg2_e", "EX_k_e", "EX_na1_e",
    "EX_ca2_e", "EX_cl_e", "EX_h2o_e", "EX_h_e", "EX_co2_e",
]

AA_EXCHANGES = [
    "EX_ala__L_e", "EX_arg__L_e", "EX_asn__L_e", "EX_asp__L_e",
    "EX_cys__L_e", "EX_gln__L_e", "EX_glu__L_e", "EX_gly_e",
    "EX_his__L_e", "EX_ile__L_e", "EX_leu__L_e", "EX_lys__L_e",
    "EX_met__L_e", "EX_phe__L_e", "EX_pro__L_e", "EX_ser__L_e",
    "EX_thr__L_e", "EX_trp__L_e", "EX_tyr__L_e", "EX_val__L_e",
]


def _ensure_metabolite(model, met_id, name, formula, charge, compartment):
    if met_id in model.metabolites:
        return
    met = Metabolite(met_id, formula=formula, name=name,
                      compartment=compartment, charge=charge)
    model.add_metabolites([met])


def _add_reaction_if_missing(model, rxn_id, reactants_str, products_str,
                              gene_rule, reversible=False):
    if rxn_id in model.reactions:
        return
    rxn = Reaction(rxn_id, name=rxn_id)
    model.add_reactions([rxn])
    arrow = "<=>" if reversible else "-->"
    rxn.build_reaction_from_string(f"{reactants_str} {arrow} {products_str}")
    if gene_rule:
        rxn.gene_reaction_rule = gene_rule


def _load_organism_model(model_path):
    """依副檔名選擇讀取方式。優先建議用 .json（cobrapy 原生格式，最穩），
    也支援 .xml/.sbml（SBML）與 .mat（COBRA Toolbox）。"""
    ext = os.path.splitext(model_path)[1].lower()
    if ext == ".json":
        return cobra.io.load_json_model(model_path)
    if ext in (".xml", ".sbml"):
        return cobra.io.read_sbml_model(model_path)
    if ext == ".mat":
        return cobra.io.load_matlab_model(model_path)
    raise ValueError(f"不支援的模型檔格式：{ext}（請用 .json / .xml / .mat）")


def build_model(model_path):
    """讀取原始模型並加上異源限烯烴路徑（對應 prepareModel009 前半段）。"""
    model = _load_organism_model(model_path)

    _ensure_metabolite(model, "limonene_c", "Limonene", "C10H16", 0, "c")
    _ensure_metabolite(model, "limonene_e", "Limonene", "C10H16", 0, "e")
    _ensure_metabolite(model, "ggdp_c", "Geranylgeranyl diphosphate",
                        "C20H33O7P2", -3, "c")

    for rxn_id, (reactants, products, gene_rule) in HET_REACTIONS.items():
        _add_reaction_if_missing(model, rxn_id, reactants, products, gene_rule)

    _add_reaction_if_missing(model, "LIMtex", "limonene_c", "limonene_e",
                              None, reversible=True)
    model.reactions.LIMtex.bounds = (0, 1000)

    _add_reaction_if_missing(model, "EX_limonene_e", "limonene_e", "",
                              None, reversible=True)
    model.reactions.EX_limonene_e.bounds = (0, 1000)

    # 所有 T7/異源反應先關閉，由每次模擬依 Etot 容量重新打開。
    for rxn_id in ALL_T7_RXNS:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).bounds = (0, 0)

    return model


def apply_medium009(model, defaults):
    """比照 prepareModel009() 的 TB 類培養基設定。"""
    for rxn in model.exchanges:
        if rxn.lower_bound < 0:
            rxn.lower_bound = 0

    model.reactions.EX_glyc_e.lower_bound = -abs(defaults["maxUptake.EX_glyc_e"])
    for rxn_id in RICH_MEDIUM_OPEN:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = -1000
    for rxn_id in TRACE_ION_RXNS:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = -1000
    model.reactions.EX_o2_e.lower_bound = -abs(defaults["maxUptake.EX_o2_e"])

    rich_aa_uptake = defaults["sim.richAminoAcidMaxUptake"]
    for rxn_id in AA_EXCHANGES:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = -abs(rich_aa_uptake)

    return model


# ============================================================
# 3. dFBA 引擎（原 dfba_engine.py）
#    移植 dynamicFBA_T7simple004.m 的 IPTG 誘導 + Etot 容量邏輯
# ============================================================

T7_ENZYME_RXNS = {
    "DXS": "DXS_T7_het",
    "IDI_fwd": "IDI_T7_IPP_to_DMAPP_het",
    "IDI_rev": "IDI_T7_DMAPP_to_IPP_het",
    "GPPS": "GPPS_S80F_het",
    "LS": "LIMS_MS_het",
}

CAPACITY_ENZYME_FOR_RXN = {
    "DXS_T7_het": "DXS",
    "IDI_T7_IPP_to_DMAPP_het": "IDI",
    "IDI_T7_DMAPP_to_IPP_het": "IDI",
    "GPPS_S80F_het": "GPPS",
    "LIMS_MS_het": "LS",
}

BURDEN_RXN = "T7_BURDEN"


def etot_to_capacity_flux(kcat_s, etot_uM, cell_volume_L_per_gDW):
    """capacityFlux [mmol/gDW/h] = 3.6 * kcat[s^-1] * Etot[uM] * cellVolume[L/gDW]"""
    vmax_mM_h = 3.6 * kcat_s * etot_uM
    return vmax_mM_h * cell_volume_L_per_gDW


def apply_t7_kinetic_rules(model, capacities, time_now, iptg_start_time,
                            pre_induction_expr, induced_expr,
                            induction_ramp_time, burden_base_lb, burden_max_lb):
    """套用 IPTG 誘導表現比例與 T7 酵素容量上限（就地修改 model）。"""
    if time_now < iptg_start_time:
        phase = "preIPTG_locked"
        expr_factor = pre_induction_expr
    else:
        phase = "postIPTG_induced"
        if induction_ramp_time > 0:
            ramp = min(1, max(0, (time_now - iptg_start_time) / induction_ramp_time))
            expr_factor = pre_induction_expr + ramp * (induced_expr - pre_induction_expr)
        else:
            expr_factor = induced_expr
    expr_factor = max(0, expr_factor)

    for rxn_id, enzyme in CAPACITY_ENZYME_FOR_RXN.items():
        if rxn_id not in model.reactions:
            continue
        rxn = model.reactions.get_by_id(rxn_id)
        capacity = max(0, capacities.get(enzyme, 0))
        rxn.lower_bound = 0
        rxn.upper_bound = capacity * expr_factor

    burden_lb = burden_base_lb + (burden_max_lb - burden_base_lb) * expr_factor
    if BURDEN_RXN in model.reactions:
        burden = model.reactions.get_by_id(BURDEN_RXN)
        burden.upper_bound = max(burden.upper_bound, burden_lb)
        burden.lower_bound = burden_lb

    return expr_factor, phase


def solve_growth_only(model, biomass_rxn):
    model.objective = biomass_rxn
    sol = model.optimize()
    return sol, "growthOnly", None


def solve_production_priority(model, biomass_rxn, prod_rxn, min_growth_frac,
                               fallback_growth_fracs, growth_floor_reference,
                               reference_mu_max, min_positive_mu,
                               stop_on_non_positive_growth):
    """在生物量下限限制下最大化生產通量，必要時逐步降低下限重試。"""
    growth_fracs = list(dict.fromkeys([min_growth_frac, *fallback_growth_fracs]))

    biomass_bounds = model.reactions.get_by_id(biomass_rxn).bounds
    model.objective = biomass_rxn
    growth_sol = model.optimize()
    if growth_sol.status != "optimal":
        model.reactions.get_by_id(biomass_rxn).bounds = biomass_bounds
        return growth_sol, "prodPriority_noCurrentGrowthReference", None, float("nan")

    available_mu_max = max(0.0, growth_sol.fluxes[biomass_rxn])
    reference_mu = available_mu_max if growth_floor_reference == "currentStepMax" \
        else reference_mu_max

    sol = growth_sol
    mode = "prodPriority_failed"
    growth_floor = None
    for frac in growth_fracs:
        floor_now = max(0.0, frac * reference_mu)
        model.reactions.get_by_id(biomass_rxn).lower_bound = floor_now
        model.objective = prod_rxn
        sol = model.optimize()
        if sol.status == "optimal":
            mu = sol.fluxes[biomass_rxn]
            acceptable = (mu > min_positive_mu) if stop_on_non_positive_growth \
                else (mu >= -abs(min_positive_mu))
            if acceptable:
                mode = f"prodPriority_frac_{frac:g}"
                growth_floor = floor_now
                break
    model.reactions.get_by_id(biomass_rxn).bounds = biomass_bounds
    return sol, mode, growth_floor, available_mu_max


def update_uptake_bounds(model, exchange_ids, concentrations, biomass,
                          time_step, original_uptake_capacity, uptake_allowed):
    for rxn_id in exchange_ids:
        rxn = model.reactions.get_by_id(rxn_id)
        if not uptake_allowed.get(rxn_id, False):
            rxn.lower_bound = 0
            continue
        conc = concentrations.get(rxn_id, 0.0)
        bound = 0.0 if biomass <= 0 else conc / (biomass * time_step)
        bound = min(bound, 1000, original_uptake_capacity.get(rxn_id, 0.0))
        if abs(bound) < 1e-12:
            bound = 0.0
        rxn.lower_bound = -bound


def run_dfba(model, params):
    """執行單次 dFBA 模擬，回傳時間序列。"""
    capacities = {
        enzyme: etot_to_capacity_flux(params["kcat"][enzyme], params["Etot"][enzyme],
                                       params["cell_volume_L_per_gDW"])
        for enzyme in ("DXS", "IDI", "GPPS", "LS")
    }

    biomass_rxn = params["biomass_rxn"]
    prod_rxn = params.get("prod_rxn", "EX_limonene_e")
    substrate_rxn = params.get("substrate_rxn", "EX_glyc_e")
    dt = params["time_step_h"]
    n_steps = params["n_steps"]

    exchange_ids = [r.id for r in model.exchanges]
    original_uptake_capacity = {r_id: max(0.0, -model.reactions.get_by_id(r_id).lower_bound)
                                 for r_id in exchange_ids}
    uptake_allowed = {r_id: original_uptake_capacity[r_id] > 0 for r_id in exchange_ids}

    concentrations = {r_id: 0.0 for r_id in exchange_ids}
    concentrations[substrate_rxn] = params["init_glycerol_mM"]
    for r_id in exchange_ids:
        if uptake_allowed[r_id] and r_id != substrate_rxn and concentrations[r_id] == 0:
            concentrations[r_id] = params["non_tracked_conc_mM"]

    biomass = params["init_biomass_gDW_L"]

    time_vec = [0.0]
    biomass_vec = [biomass]
    substrate_vec = [concentrations[substrate_rxn]]
    limonene_vec = [concentrations.get(prod_rxn, 0.0)]
    expr_vec = [None]
    growth_rate_vec = [None]
    limonene_flux_vec = [None]
    t7_flux_log = {name: [None] for name in T7_ENZYME_RXNS}

    update_uptake_bounds(model, exchange_ids, concentrations, biomass, dt,
                          original_uptake_capacity, uptake_allowed)

    reference_mu_max = float("nan")
    stopped_reason = None

    for step in range(n_steps):
        t_now = time_vec[-1]

        expr_factor, phase = apply_t7_kinetic_rules(
            model, capacities, t_now,
            params["iptg_start_time_h"], params["pre_induction_expr"],
            params["induced_expr"], params["induction_ramp_time_h"],
            params["burden_base_lb"], params["burden_max_lb"])

        is_induced = t_now >= params["iptg_start_time_h"]
        available_mu_max = float("nan")

        if is_induced and params.get("use_prod_priority", True) and expr_factor > 0:
            sol, mode, growth_floor, available_mu_max = solve_production_priority(
                model, biomass_rxn, prod_rxn, params["min_growth_frac"],
                params["fallback_growth_fracs"], params["growth_floor_reference"],
                reference_mu_max, params["min_positive_mu"],
                params["stop_on_non_positive_growth"])
        else:
            sol, mode, growth_floor = solve_growth_only(model, biomass_rxn)

        if sol.status != "optimal" and is_induced and params.get("allow_growth_fallback", True):
            sol, mode, growth_floor = solve_growth_only(model, biomass_rxn)
            mode = "growthFallback"

        if sol.status != "optimal":
            stopped_reason = f"no feasible solution at step {step}, t={t_now:.2f}h"
            break

        mu = sol.fluxes[biomass_rxn]
        if not math.isfinite(available_mu_max):
            available_mu_max = mu
        reference_mu_max = available_mu_max if math.isfinite(available_mu_max) else reference_mu_max

        if mu <= params["min_positive_mu"]:
            if params["stop_on_non_positive_growth"]:
                stopped_reason = f"no positive growth at step {step}, t={t_now:.2f}h"
                break
            mu = max(mu, 0.0)

        uptake_flux = {r_id: sol.fluxes.get(r_id, 0.0) for r_id in exchange_ids}

        # 解析解質量平衡更新（與 MATLAB 版一致：對指數成長積分）。
        if abs(mu) < 1e-12:
            delta_conc = {r_id: uptake_flux[r_id] * biomass * dt for r_id in exchange_ids}
            new_biomass = biomass
        else:
            growth_integral = biomass * (math.exp(mu * dt) - 1) / mu
            delta_conc = {r_id: uptake_flux[r_id] * growth_integral for r_id in exchange_ids}
            new_biomass = biomass * math.exp(mu * dt)

        new_conc = {}
        for r_id in exchange_ids:
            c = concentrations[r_id] + delta_conc[r_id]
            new_conc[r_id] = 0.0 if abs(c) < 1e-12 or c < 0 else c

        biomass = new_biomass
        concentrations = new_conc
        t_next = t_now + dt

        time_vec.append(t_next)
        biomass_vec.append(biomass)
        substrate_vec.append(concentrations[substrate_rxn])
        limonene_vec.append(concentrations.get(prod_rxn, 0.0))
        expr_vec.append(expr_factor)
        growth_rate_vec.append(mu)
        limonene_flux_vec.append(sol.fluxes.get(prod_rxn, 0.0))
        for name, rxn_id in T7_ENZYME_RXNS.items():
            t7_flux_log[name].append(sol.fluxes.get(rxn_id))

        update_uptake_bounds(model, exchange_ids, concentrations, biomass, dt,
                              original_uptake_capacity, uptake_allowed)

    return {
        "time_h": time_vec,
        "biomass_gDW_L": biomass_vec,
        "glycerol_mM": substrate_vec,
        "limonene_mM": limonene_vec,
        "expr_factor": expr_vec,
        "growth_rate_h": growth_rate_vec,
        "limonene_flux": limonene_flux_vec,
        "t7_fluxes": t7_flux_log,
        "capacities_mmol_gDW_h": capacities,
        "stopped_reason": stopped_reason,
    }


def etot_scan(model, base_params, enzyme, scan_values_uM):
    """對單一酵素做 one-at-a-time Etot 掃描，其餘酵素固定在基準值。"""
    results = []
    for value in scan_values_uM:
        params = dict(base_params)
        params["Etot"] = dict(base_params["Etot"])
        params["Etot"][enzyme] = value
        with model:
            trace = run_dfba(model, params)
        final_limonene = trace["limonene_mM"][-1]
        fluxes = [f for f in trace["limonene_flux"] if f is not None]
        max_flux = max(fluxes) if fluxes else 0.0
        results.append({
            "etot_uM": value,
            "trace": trace,
            "final_limonene_mM": final_limonene,
            "max_post_iptg_flux": max_flux,
        })
    return results


def etot_scan_2d(model, base_params, enzyme_x, enzyme_y, values_x, values_y):
    """雙酵素交叉掃描（two-factor grid）：enzyme_x/enzyme_y 各自的 Etot 網格組合，
    其餘酵素固定在基準值。回傳依 (y, x) 排列的網格矩陣，方便前端畫 heatmap。"""
    final_grid = []
    flux_grid = []
    for y_value in values_y:
        final_row = []
        flux_row = []
        for x_value in values_x:
            params = dict(base_params)
            params["Etot"] = dict(base_params["Etot"])
            params["Etot"][enzyme_x] = x_value
            params["Etot"][enzyme_y] = y_value
            with model:
                trace = run_dfba(model, params)
            final_row.append(trace["limonene_mM"][-1])
            fluxes = [f for f in trace["limonene_flux"] if f is not None]
            flux_row.append(max(fluxes) if fluxes else 0.0)
        final_grid.append(final_row)
        flux_grid.append(flux_row)
    return {
        "enzyme_x": enzyme_x,
        "enzyme_y": enzyme_y,
        "values_x": values_x,
        "values_y": values_y,
        "final_limonene_mM": final_grid,
        "max_post_iptg_flux": flux_grid,
    }


# ============================================================
# 4. Flask API（原 app.py）
# ============================================================

def _default_model_path(base_dir):
    """在 models/ 資料夾裡依序找 .json > .xml > .mat，回傳第一個存在的檔案。"""
    models_dir = os.path.join(base_dir, "models")
    for ext in (".json", ".xml", ".mat"):
        candidate = os.path.join(models_dir, f"iEC1356_Bl21DE3{ext}")
        if os.path.exists(candidate):
            return candidate
    # 都找不到時，回傳 .json 路徑，讓錯誤訊息明確指出要放哪個檔案。
    return os.path.join(models_dir, "iEC1356_Bl21DE3.json")


BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.environ.get("LIMONENE_MODEL_PATH", _default_model_path(BASE_DIR))
PARAM_CSV = os.environ.get(
    "LIMONENE_PARAM_CSV", os.path.join(BASE_DIR, "LIM009_params.csv"))

app = Flask(__name__, static_folder=BASE_DIR, static_url_path="")

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
                f"找不到菌株模型檔：{MODEL_PATH}\n"
                "請把 iEC1356_Bl21DE3.json（建議）或 .xml / .mat 放到 models/ 資料夾，"
                "或設定 LIMONENE_MODEL_PATH 環境變數指到正確路徑。")
        model = build_model(MODEL_PATH)
        model = apply_medium009(model, get_registry().to_dfba_defaults())
        biomass_rxns = [r.id for r in model.reactions if r.objective_coefficient != 0]
        model._biomass_rxn = biomass_rxns[0] if biomass_rxns else None
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


@app.route("/")
def index():
    return send_from_directory(BASE_DIR, "index.html")


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
    """穩態 FBA：先求最大生長速率，再以 min_growth_frac 為下限最大化限烯烴分泌。"""
    body = request.get_json(force=True, silent=True) or {}
    params = merge_params(body.get("params"))
    induced = bool(body.get("induced", True))

    model = working_model()
    capacities = {
        enzyme: params["kcat"][enzyme] * params["Etot"][enzyme] * 3.6
                * params["cell_volume_L_per_gDW"]
        for enzyme in ("DXS", "IDI", "GPPS", "LS")
    }
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
                         "message": "生長最佳化不可行，請檢查培養基或 Etot 設定"})

    mu_max = growth_sol.fluxes[biomass_rxn]
    model.reactions.get_by_id(biomass_rxn).lower_bound = params["min_growth_frac"] * mu_max
    model.objective = "EX_limonene_e"
    prod_sol = model.optimize()

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
    })


@app.route("/api/dfba", methods=["POST"])
def dfba():
    """單次動態 FBA 模擬，回傳完整時間序列。"""
    body = request.get_json(force=True, silent=True) or {}
    params = merge_params(body.get("params"))
    model = working_model()
    result = run_dfba(model, params)
    return jsonify(result)


@app.route("/api/etot-scan", methods=["POST"])
def etot_scan_endpoint():
    """對單一酵素做 Etot one-at-a-time 掃描（對應 LimoneneCOBRA009 的 OAT scan）。"""
    body = request.get_json(force=True, silent=True) or {}
    enzyme = body.get("enzyme", "LS")
    if enzyme not in ("DXS", "IDI", "GPPS", "LS"):
        return jsonify({"error": "enzyme 必須是 DXS、IDI、GPPS 或 LS 其中之一"}), 400

    params = merge_params(body.get("params"))
    scan_values = body.get("scan_values") or get_registry().scan_values(enzyme)

    model = get_base_model()
    results = etot_scan(model, params, enzyme, scan_values)
    return jsonify({"enzyme": enzyme, "results": results})


@app.route("/api/etot-scan-2d", methods=["POST"])
def etot_scan_2d_endpoint():
    """對兩個酵素做交叉（雙因子）Etot 掃描，回傳網格結果。"""
    body = request.get_json(force=True, silent=True) or {}
    enzyme_x = body.get("enzyme_x", "LS")
    enzyme_y = body.get("enzyme_y", "DXS")
    valid_enzymes = ("DXS", "IDI", "GPPS", "LS")
    if enzyme_x not in valid_enzymes or enzyme_y not in valid_enzymes:
        return jsonify({"error": "enzyme_x / enzyme_y 必須是 DXS、IDI、GPPS 或 LS 其中之一"}), 400
    if enzyme_x == enzyme_y:
        return jsonify({"error": "enzyme_x 與 enzyme_y 不可相同"}), 400

    params = merge_params(body.get("params"))
    values_x = body.get("values_x") or get_registry().scan_values(enzyme_x)
    values_y = body.get("values_y") or get_registry().scan_values(enzyme_y)

    model = get_base_model()
    result = etot_scan_2d(model, params, enzyme_x, enzyme_y, values_x, values_y)
    return jsonify(result)


if __name__ == "__main__":
    app.run(port=5001, debug=True)