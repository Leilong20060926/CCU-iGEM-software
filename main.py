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

加速：
    pip install highspy
    → 裝了就會自動改用 HiGHS 這個 LP solver（比 cobrapy 預設的 GLPK 快很多倍），
      對 FBA / dFBA / 掃描 / 2D 掃描全部都有幫助，預設就是啟用的。裝不到時會
      安靜地退回 GLPK，不影響任何功能。
      如果你需要跟 MATLAB／GLPK 的結果逐位對齊，設環境變數
      LIMONENE_FAST_SOLVER=0 可以停用、強制使用 GLPK。
    另外，Etot 掃描與雙酵素 2D 掃描的每一格彼此獨立，會自動用多個 CPU 核心平行
    計算（背景用 ProcessPoolExecutor），不需要額外設定，這個不影響求解結果。
"""

import concurrent.futures
import csv
import math
import os

import cobra
from cobra import Reaction, Metabolite
from cobra.flux_analysis import pfba
from flask import Flask, jsonify, request, send_from_directory
from scipy.optimize import minimize


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

    def scalar_or_default(self, code, default):
        """同 scalar()，但若 CSV 沒有這個 parameter_code 就回傳預設值，不拋例外。
        用於尚未登錄進 LIM009_params.csv 的可調參數（如 ATPM / DXPS / METAT / FPPS 邊界）。"""
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
            # 以下為可由前端調整的原生代謝路徑邊界，若 LIM009_params.csv 尚未登錄
            # 對應 parameter_code，就退回預設值（不會拋例外）。
            "atpm_lb": r.scalar_or_default("sim.atpmLB", 4.0),
            "dxps_ub": r.scalar_or_default("bounds.DXPS_ub", 30.0),
            "metat_ub": r.scalar_or_default("bounds.METAT_ub", 0.10),
            # 注意：LimoneneCOBRA001.m 原始腳本設的是 0.05（knock-down 但不完全
            # 阻斷），跟你之前給的參考報表數字 0.00 不一致，這裡先照腳本改成
            # 0.05——如果你確認 0.00 才是對的，請告訴我，我再改回去。
            "fpps_ub": r.scalar_or_default("bounds.FPPS_ub", 0.05),
            # dFBA 時間步長模式：使用者可選 "fixed"（預設，對照 MATLAB 參考結果）
            # 或 "adaptive"（自適應步長，換速度、數值不再逐位對齊 MATLAB）。
            "step_mode": "fixed",
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

# 對照 LimoneneCOBRA001.m 的 inorganic_rxns：多了 EX_nh4_e（銨根，主要氮源之一）。
# 原本這裡漏掉這一項，銨根會被前面「全部先歸零」那段鎖死在 0，等於完全沒有
# 這個氮源可用——這是造成生長速率跟參考版本兜不起來的原因之一。
RICH_MEDIUM_OPEN = [
    "EX_pi_e", "EX_nh4_e", "EX_so4_e", "EX_mg2_e", "EX_k_e", "EX_na1_e",
    "EX_ca2_e", "EX_cl_e", "EX_h2o_e", "EX_h_e", "EX_co2_e",
]

AA_EXCHANGES = [
    "EX_ala__L_e", "EX_arg__L_e", "EX_asn__L_e", "EX_asp__L_e",
    "EX_cys__L_e", "EX_gln__L_e", "EX_glu__L_e", "EX_gly_e",
    "EX_his__L_e", "EX_ile__L_e", "EX_leu__L_e", "EX_lys__L_e",
    "EX_met__L_e", "EX_phe__L_e", "EX_pro__L_e", "EX_ser__L_e",
    "EX_thr__L_e", "EX_trp__L_e", "EX_tyr__L_e", "EX_val__L_e",
]

# 對照 LimoneneCOBRA001.m 的 aa_ratios：BioShop TB（Tryptone 12g/L + Yeast Extract
# 24g/L）的胺基酸不是每種都給一樣的攝取上限，而是按這個比例分配一個「總量」。
# 順序必須跟 AA_EXCHANGES 一一對應（Ala, Arg, Asn, Asp, Cys, Gln, Glu, Gly, His,
# Ile, Leu, Lys, Met, Phe, Pro, Ser, Thr, Trp, Tyr, Val）。
_AA_RATIOS_RAW = [
    0.03, 0.03, 0.04, 0.07,
    0.005, 0.05, 0.15, 0.03,
    0.02, 0.05, 0.10, 0.08,
    0.025, 0.05, 0.10, 0.06,
    0.04, 0.01, 0.03, 0.07,
]
_AA_RATIOS_SUM = sum(_AA_RATIOS_RAW)
AA_RATIOS = dict(zip(AA_EXCHANGES, (r / _AA_RATIOS_SUM for r in _AA_RATIOS_RAW)))

MEDIUM_NAME = "BioShop Terrific Broth (Tryptone 12g/L, Yeast Extract 24g/L)"

# 對照 LimoneneCOBRA001.m：腳本明確寫死要用哪一個 biomass 反應，不是自動偵測。
# iEC1356 這類模型常常同時有好幾個 BIOMASS_ 開頭的變體（例如 _WT_ / _core_），
# 如果自動偵測（挑第一個 objective_coefficient != 0 的反應）跟腳本指定的不是
# 同一個，算出來的最大生長速率可以差到快 2 倍。優先用這個環境變數／固定 ID，
# 模型裡真的找不到才 fallback 回自動偵測。
PREFERRED_BIOMASS_RXN = os.environ.get(
    "LIMONENE_BIOMASS_RXN", "BIOMASS_Ec_iJO1366_WT_53p95M")


def detect_biomass_rxn(model):
    """優先使用 PREFERRED_BIOMASS_RXN（對齊 LimoneneCOBRA001.m 寫死的目標反應），
    模型裡沒有這個 ID 才 fallback 回『目標函數係數不為 0 的第一個反應』。"""
    if PREFERRED_BIOMASS_RXN in model.reactions:
        return PREFERRED_BIOMASS_RXN
    biomass_rxns = [r.id for r in model.reactions if r.objective_coefficient != 0]
    return biomass_rxns[0] if biomass_rxns else None


# 原生代謝路徑邊界可調的反應 ID。MEP_PATHWAY_RXNS 對照 LimoneneCOBRA001.m 的
# mep_rxns：原本我們只有 DXPS 一個反應吃到「MEP 路徑上限」這個參數，實際上
# MEPCT / MECDPS / IPDPI 這三個反應也要套用同一個上限值，只加 DXPS 會讓 MEP
# 路徑後段沒有跟著放寬，通量被卡住。
MEP_PATHWAY_RXNS = ("DXPS", "MEPCT", "MECDPS", "IPDPI")
NATIVE_PATHWAY_BOUNDS = ("DXPS", "METAT", "FPPS")


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

    _use_fast_solver(model)
    return model


def _use_fast_solver(model):
    """嘗試切換到比預設 GLPK 快很多的 LP solver（HiGHS）。

    預設會啟用（不需要再另外設環境變數）。之前刻意設成「要手動開啟」是因為
    當時還在跟 MATLAB／GLPK 逐位對照驗證數字，換 solver 在 LP 解不唯一
    （degenerate solution，pFBA 場景下很常見）時可能挑到不同的那一組解，
    個別反應通量會跟 GLPK 跑出來的參考結果有落差。現在對照驗證已經完成，
    改成預設啟用。如果你之後又需要跟 GLPK 逐位對齊，設環境變數
    LIMONENE_FAST_SOLVER=0 停用即可退回 GLPK。"""
    if os.environ.get("LIMONENE_FAST_SOLVER", "1").strip() in ("0", "false", "False"):
        return None
    for solver_name in ("highs", "cplex", "gurobi"):
        try:
            model.solver = solver_name
            return solver_name
        except Exception:
            continue
    return None


def apply_medium009(model, defaults):
    """比照 prepareModel009() 的 TB 類培養基設定（BioShop Terrific Broth Formulation）。"""
    for rxn in model.exchanges:
        if rxn.lower_bound < 0:
            rxn.lower_bound = 0

    # 明確封鎖葡萄糖攝取（主要碳源改用甘油）。
    if "EX_glc__D_e" in model.reactions:
        model.reactions.EX_glc__D_e.lower_bound = 0

    # 限制甘油攝取（主要碳源），數值來自 LIM009_params.csv 的 maxUptake.EX_glyc_e。
    model.reactions.EX_glyc_e.lower_bound = -abs(defaults["maxUptake.EX_glyc_e"])
    model.reactions.EX_glyc_e.upper_bound = 1000

    for rxn_id in RICH_MEDIUM_OPEN:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = -1000
    for rxn_id in TRACE_ION_RXNS:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = -1000

    # 好氧呼吸上限，數值來自 LIM009_params.csv 的 maxUptake.EX_o2_e。
    model.reactions.EX_o2_e.lower_bound = -abs(defaults["maxUptake.EX_o2_e"])
    model.reactions.EX_o2_e.upper_bound = 1000

    apply_amino_acid_ratios(model, defaults["sim.richAminoAcidMaxUptake"])

    # 對照 LimoneneCOBRA001.m 的「micro-allowance」：所有其餘還鎖在下限 0 的
    # 交換反應（沒被上面任何一組清單明確開放），除了葡萄糖跟限烯烴分泌以外，
    # 統一給一個極小的攝取下限 -0.01，避免某些沒被明確列出的微量代謝物
    # 造成 biomass 反應的前驅物缺口，使模型不可行或生長速率被異常壓低。
    for rxn in model.exchanges:
        if (rxn.lower_bound == 0
                and rxn.id not in ("EX_glc__D_e", "EX_limonene_e")):
            rxn.lower_bound = -0.01

    return model


def apply_amino_acid_ratios(model, total_aa_uptake):
    """對照 LimoneneCOBRA001.m 的 BioShop TB 胺基酸比例限制：不是每種胺基酸都給
    同一個攝取上限，而是把 total_aa_uptake（mmol/gDW/h，代表 Tryptone/Yeast
    Extract 提供的有機氮源總量）按 AA_RATIOS 的比例分配到 20 種胺基酸交換反應。"""
    total = abs(total_aa_uptake)
    for rxn_id, ratio in AA_RATIOS.items():
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = -total * ratio


def apply_mep_pathway_bound(model, ub):
    """對照 LimoneneCOBRA001.m 的 mep_rxns：DXPS/MEPCT/MECDPS/IPDPI 這四個反應
    共用同一個 MEP 路徑上限，只放寬 DXPS 一個會讓路徑後段卡住。"""
    for rxn_id in MEP_PATHWAY_RXNS:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).upper_bound = ub


def apply_adjustable_medium(model, params):
    """套用可由前端調整的培養基攝取限制與原生路徑邊界（每次請求時套用在傳入的
    working model / cobra context 上，不會影響共用的 base model 快取）。
    對應可調參數：甘油/氧氣攝取上限、豐富培養基胺基酸攝取總量、ATPM 維持能量
    下限，以及 MEP 路徑（DXPS/MEPCT/MECDPS/IPDPI）／METAT／FPPS 的上限
    （模型中不存在的反應會安全略過）。"""
    if "EX_glyc_e" in model.reactions:
        model.reactions.EX_glyc_e.lower_bound = -abs(params["maxUptake.EX_glyc_e"])
        model.reactions.EX_glyc_e.upper_bound = 1000
    if "EX_o2_e" in model.reactions:
        model.reactions.EX_o2_e.lower_bound = -abs(params["maxUptake.EX_o2_e"])
        model.reactions.EX_o2_e.upper_bound = 1000

    apply_amino_acid_ratios(model, params["sim.richAminoAcidMaxUptake"])

    if "ATPM" in model.reactions:
        model.reactions.ATPM.lower_bound = params["atpm_lb"]
    apply_mep_pathway_bound(model, params["dxps_ub"])
    if "METAT" in model.reactions:
        model.reactions.METAT.upper_bound = params["metat_ub"]
    if "FPPS" in model.reactions:
        model.reactions.FPPS.upper_bound = params["fpps_ub"]

    return model


def metabolite_gross_production_rate(model, solution, met_id):
    """計算某代謝物在目前解下的『總生產速率』：加總所有對該代謝物淨貢獻為正
    （即生產方向）的反應通量 × 化學計量係數。用於 ATP / NADPH 的 cofactor balance 分析。
    若模型中沒有該代謝物 ID，回傳 None。"""
    if met_id not in model.metabolites:
        return None
    met = model.metabolites.get_by_id(met_id)
    total = 0.0
    for rxn in met.reactions:
        coeff = rxn.metabolites[met]
        flux = solution.fluxes.get(rxn.id, 0.0)
        contribution = coeff * flux
        if contribution > 0:
            total += contribution
    return total


def build_medium_pathway_report(model, capacities):
    """組出培養基與路徑邊界報告用的資料（供前端渲染成文字報表）。
    找不到的反應回傳 None，前端會顯示為 N/A，不會因模型缺少該反應而出錯。
    GPPS/LIMS_MS_het 這兩個直接讀模型上「套用完誘導規則之後」的實際上限
    （而不是套用前的 kcat×Etot 理論容量），這樣報表數字才會跟真正拿去解 LP
    的邊界一致——如果只讀理論容量，遇到 GPPS_FIXED_CAPACITY 這種固定值覆蓋
    的情況，報表會跟實際生效的上限對不起來。"""
    def bounds_of(rxn_id):
        if rxn_id in model.reactions:
            r = model.reactions.get_by_id(rxn_id)
            return [r.lower_bound, r.upper_bound]
        return None

    def ub_of(rxn_id):
        b = bounds_of(rxn_id)
        return b[1] if b else None

    def lb_of(rxn_id):
        b = bounds_of(rxn_id)
        return b[0] if b else None

    return {
        "medium_name": MEDIUM_NAME,
        "glycerol_bounds": bounds_of("EX_glyc_e"),
        "oxygen_bounds": bounds_of("EX_o2_e"),
        "atpm_lb": lb_of("ATPM"),
        "dxps_ub": ub_of("DXPS"),
        "metat_ub": ub_of("METAT"),
        "fpps_ub": ub_of("FPPS"),
        "gpps_ub": ub_of("GPPS_S80F_het") if "GPPS_S80F_het" in model.reactions else capacities.get("GPPS"),
        "lims_ub": ub_of("LIMS_MS_het") if "LIMS_MS_het" in model.reactions else capacities.get("LS"),
    }


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

# GPPS/LIMS_MS_het 改用固定的酵素容量上限（mmol/gDW/h），不再用 kcat×Etot
# 動態計算——對齊使用者提供的參考模型（GPPS/LIMS_MS_het Upper Bound = 20.00）。
# 誘導比例（expr_factor）依然套用，只是「誘導完全開啟時的上限」固定在這裡，
# 不受 Etot(GPPS) / Etot(LS) 這兩個輸入欄位影響。DXS/IDI 維持原本 kcat×Etot 邏輯。
FIXED_ENZYME_CAPACITY_MMOL_GDW_H = {
    "GPPS": 20.0,
    "LS": 20.0,
}

BURDEN_RXN = "T7_BURDEN"

# 對照 dynamicFBA_T7simple004.m 的 exclUptakeRxns 預設值：這四個反應（氣體／溶劑）
# 不納入動態濃度追蹤，邊界維持培養基設定時的固定值，不會被每步的質量平衡更新去動。
# 如果誤把它們當一般「未追蹤受質」給虛擬池，長時間模擬下這個池可能被耗盡，
# 造成氧氣／二氧化碳等攝取被不該有的限制——這是跟原始 dFBA 引擎行為的真實差異。
DYNAMIC_TRACKING_EXCLUDED_RXNS = {"EX_co2_e", "EX_o2_e", "EX_h2o_e", "EX_h_e"}


def etot_to_capacity_flux(kcat_s, etot_uM, cell_volume_L_per_gDW):
    """capacityFlux [mmol/gDW/h] = 3.6 * kcat[s^-1] * Etot[uM] * cellVolume[L/gDW]"""
    vmax_mM_h = 3.6 * kcat_s * etot_uM
    return vmax_mM_h * cell_volume_L_per_gDW


def compute_enzyme_capacities(params):
    """算出四個 T7 酵素的容量上限（mmol/gDW/h，尚未乘上誘導比例）。
    GPPS/LS 直接套用 FIXED_ENZYME_CAPACITY_MMOL_GDW_H 的固定值，不用 kcat×Etot；
    這裡統一處理，確保 API 回傳的 capacities_mmol_gDW_h 跟報表、跟實際套用到
    模型上的上限三處數字一致，不會各說各話。"""
    capacities = {}
    for enzyme in ("DXS", "IDI", "GPPS", "LS"):
        if enzyme in FIXED_ENZYME_CAPACITY_MMOL_GDW_H:
            capacities[enzyme] = FIXED_ENZYME_CAPACITY_MMOL_GDW_H[enzyme]
        else:
            capacities[enzyme] = etot_to_capacity_flux(
                params["kcat"][enzyme], params["Etot"][enzyme],
                params["cell_volume_L_per_gDW"])
    return capacities


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
        if enzyme in FIXED_ENZYME_CAPACITY_MMOL_GDW_H:
            capacity = FIXED_ENZYME_CAPACITY_MMOL_GDW_H[enzyme]
        else:
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
    # 對照 dynamicFBA_T7simple004.m 的 optimizeCbModel(modelStep,'max','one')：
    # dFBA 每一步的求解本來就該用 pFBA（'one' 旗標），不是只有靜態 FBA。
    try:
        sol = pfba(model)
    except Exception:
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
    try:
        growth_sol = pfba(model)
    except Exception:
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
        try:
            sol = pfba(model)
        except Exception:
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


# ============================================================
# 3a. Etot 自動最佳化（Nelder-Mead simplex，Nelder & Mead 1965）
#     用少量幾次穩態 FBA 求解取代窮舉網格掃描，直接找出讓限烯烴通量最大化的
#     Etot(DXS/IDI/GPPS/LS) 組合。跟 dFBA 引擎完全獨立、不影響任何已驗證過
#     的計算邏輯，純粹是新增的「幫你找參數」功能。
# ============================================================

def _fba_limonene_flux_for_etot(etot_vec, base_params):
    """給定一組 Etot（DXS/IDI/GPPS/LS，µM），跑一次穩態誘導 FBA，回傳限烯烴
    分泌通量。任何一步不可行都回傳 0（讓最佳化演算法把這個方向當作壞方向，
    不會讓整個搜尋因為單一不可行點而中斷）。"""
    params = dict(base_params)
    params["Etot"] = {
        "DXS": max(0.0, etot_vec[0]), "IDI": max(0.0, etot_vec[1]),
        "GPPS": max(0.0, etot_vec[2]), "LS": max(0.0, etot_vec[3]),
    }
    model = working_model()
    apply_adjustable_medium(model, params)
    capacities = compute_enzyme_capacities(params)
    apply_t7_kinetic_rules(
        model, capacities, params["iptg_start_time_h"],
        params["iptg_start_time_h"], params["pre_induction_expr"],
        params["induced_expr"], params["induction_ramp_time_h"],
        params["burden_base_lb"], params["burden_max_lb"])

    biomass_rxn = params["biomass_rxn"]
    model.objective = biomass_rxn
    growth_sol = model.optimize()
    if growth_sol.status != "optimal":
        return 0.0
    mu_max = growth_sol.fluxes[biomass_rxn]
    model.reactions.get_by_id(biomass_rxn).lower_bound = params["min_growth_frac"] * mu_max
    model.objective = "EX_limonene_e"
    sol = model.optimize()
    if sol.status != "optimal":
        return 0.0
    return max(0.0, sol.fluxes.get("EX_limonene_e", 0.0))


def optimize_etot_nelder_mead(base_params, initial_etot, max_iter=200):
    """用 Nelder-Mead 找出讓限烯烴通量最大化的 Etot 組合。
    scipy 的 minimize 只會找最小值，所以內部目標函式回傳「負的限烯烴通量」，
    等於在找最大值。"""
    x0 = [initial_etot["DXS"], initial_etot["IDI"], initial_etot["GPPS"], initial_etot["LS"]]
    history = []

    def objective(x):
        flux = _fba_limonene_flux_for_etot(x, base_params)
        history.append({
            "etot": {
                "DXS": max(0.0, x[0]), "IDI": max(0.0, x[1]),
                "GPPS": max(0.0, x[2]), "LS": max(0.0, x[3]),
            },
            "limonene_flux": flux,
        })
        return -flux

    result = minimize(
        objective, x0, method="Nelder-Mead",
        options={"maxiter": max_iter, "xatol": 1e-3, "fatol": 1e-6, "adaptive": True})

    best_x = [max(0.0, v) for v in result.x]
    return {
        "best_etot": {"DXS": best_x[0], "IDI": best_x[1], "GPPS": best_x[2], "LS": best_x[3]},
        "best_limonene_flux": -result.fun,
        "n_evaluations": len(history),
        "converged": bool(result.success),
        "history": history,
    }


def run_dfba(model, params):
    """執行單次 dFBA 模擬，回傳時間序列。"""
    capacities = compute_enzyme_capacities(params)

    biomass_rxn = params["biomass_rxn"]
    prod_rxn = params.get("prod_rxn", "EX_limonene_e")
    substrate_rxn = params.get("substrate_rxn", "EX_glyc_e")
    dt = params["time_step_h"]
    n_steps = params["n_steps"]

    exchange_ids = [r.id for r in model.exchanges
                    if r.id not in DYNAMIC_TRACKING_EXCLUDED_RXNS]
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

    # step_mode 讓使用者選擇時間步長怎麼走：
    #   "fixed"（預設）：固定步長 dt = time_step_h，跑滿 n_steps 步——完全對照
    #     dynamicFBA_T7simple004.m 的行為，數字跟 MATLAB 參考結果一致。
    #   "adaptive"：步長會依生長速率的變化率動態縮放（變化大時自動縮小步長，
    #     像誘導剛拉滿那個瞬間；變化平緩時自動放大步長），數字「不會」完全
    #     等於固定步長版本——這是刻意的取捨（用更少步數跑到同樣的模擬總時長，
    #     換取速度，數值上是另一種近似），選這個模式代表你不需要跟 MATLAB
    #     逐位對齊。
    step_mode = params.get("step_mode", "fixed")
    target_total_time_h = n_steps * dt
    dt_min = dt / 10.0
    dt_max = dt * 4.0
    adaptive_target_rel_change = 0.08
    max_adaptive_steps = n_steps * 8
    prev_mu = None
    step = 0

    while True:
        if step_mode == "adaptive":
            if time_vec[-1] >= target_total_time_h - 1e-9 or step >= max_adaptive_steps:
                break
            dt = min(dt, target_total_time_h - time_vec[-1])
        else:
            if step >= n_steps:
                break

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

        if step_mode == "adaptive":
            # 依這一步生長速率的相對變化量，決定「下一步」的步長：變化越大，
            # 下一步縮得越小（例如誘導剛拉滿那瞬間）；變化越平緩，下一步放大。
            if prev_mu is not None:
                rel_change = abs(mu - prev_mu) / max(abs(prev_mu), 1e-9)
                factor = math.sqrt(adaptive_target_rel_change / max(rel_change, 1e-9))
                factor = min(2.0, max(0.5, factor))
                dt = min(dt_max, max(dt_min, dt * factor))
            prev_mu = mu

        step += 1

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


# ============================================================
# 3b. FSEOF — Flux Scanning based on Enforced Objective Flux
#     （Choi et al. 2010, Biotechnol Bioeng；Park et al. 2012）
#     公開發表的過量表現／下調基因標的搜尋演算法。做法：
#     逐步拉高目標反應（例如限烯烴分泌）的強制通量下限，每一步都重新最大化
#     生長速率、記錄全模型的通量分布；隨著強制通量一路拉高，通量「持續同向
#     變化」的反應就是候選標的：持續上升 → 過量表現候選，持續下降 → 下調候選。
# ============================================================

def fseof_scan(model, biomass_rxn, target_rxn, n_steps=10, enforce_frac=0.9):
    """執行 FSEOF 掃描，回傳每一步強制通量下的全模型通量分布，供後續分類。"""
    with model:
        model.objective = biomass_rxn
        growth_sol = model.optimize()
        if growth_sol.status != "optimal":
            raise ValueError("生長最佳化不可行，無法執行 FSEOF")
        mu_max = growth_sol.fluxes[biomass_rxn]
        v_start = max(0.0, growth_sol.fluxes.get(target_rxn, 0.0))

    with model:
        # 保留一點點生長下限，避免掃到「完全不長、只拚產物」的退化解。
        model.reactions.get_by_id(biomass_rxn).lower_bound = 0.1 * mu_max
        model.objective = target_rxn
        target_sol = model.optimize()
        if target_sol.status != "optimal":
            raise ValueError("目標反應最佳化不可行，無法執行 FSEOF")
        v_theoretical_max = target_sol.fluxes[target_rxn]

    v_end = enforce_frac * v_theoretical_max
    if v_end <= v_start or n_steps < 2:
        raise ValueError("目標反應的可強制拉高空間太小，無法進行 FSEOF 掃描"
                          "（請確認目前參數下這個反應本來就有機會再提高）")

    enforced_values = [
        v_start + (v_end - v_start) * i / (n_steps - 1) for i in range(n_steps)
    ]

    flux_matrix = {}
    for level in enforced_values:
        with model:
            model.reactions.get_by_id(target_rxn).lower_bound = level
            model.objective = biomass_rxn
            sol = model.optimize()
            if sol.status != "optimal":
                for rxn in model.reactions:
                    flux_matrix.setdefault(rxn.id, []).append(None)
                continue
            for rxn in model.reactions:
                flux_matrix.setdefault(rxn.id, []).append(float(sol.fluxes.get(rxn.id, 0.0)))

    return {
        "enforced_values": enforced_values,
        "v_start": v_start,
        "v_theoretical_max": v_theoretical_max,
        "flux_matrix": flux_matrix,
    }


def classify_fseof_targets(fseof_result, biomass_rxn, target_rxn, min_flux_magnitude=1e-6):
    """依 FSEOF 掃描結果分類每個反應：'up'（過量表現候選）或 'down'（下調候選）。
    判斷方式：看逐步之間的通量變化是不是幾乎都同一個方向（預設要求至少 80% 的
    步驟同向），並用「趨勢一致程度」與「首尾通量變化量」排序。排除目標反應本身、
    biomass 反應，以及全程通量幾乎都是 0 的反應（沒有實質意義的雜訊）。"""
    flux_matrix = fseof_result["flux_matrix"]
    results = []
    for rxn_id, fluxes in flux_matrix.items():
        if rxn_id in (biomass_rxn, target_rxn):
            continue
        clean = [f for f in fluxes if f is not None]
        if len(clean) < 2:
            continue
        if max(abs(f) for f in clean) < min_flux_magnitude:
            continue

        diffs = [clean[i + 1] - clean[i] for i in range(len(clean) - 1)]
        n_up = sum(1 for d in diffs if d > 1e-9)
        n_down = sum(1 for d in diffs if d < -1e-9)
        n_total = len(diffs)
        if n_total == 0:
            continue

        consistency = max(n_up, n_down) / n_total
        if consistency < 0.8:
            continue

        direction = "up" if n_up >= n_down else "down"
        results.append({
            "reaction": rxn_id,
            "direction": direction,
            "consistency": consistency,
            "delta_flux": clean[-1] - clean[0],
            "start_flux": clean[0],
            "end_flux": clean[-1],
        })

    results.sort(key=lambda r: (-r["consistency"], -abs(r["delta_flux"])))
    return results


def etot_scan(model, base_params, enzyme, scan_values_uM):
    """對單一酵素做 one-at-a-time Etot 掃描，其餘酵素固定在基準值。
    每個掃描值彼此獨立（互不影響），交給 run_params_grid 平行計算。"""
    params_list = []
    for value in scan_values_uM:
        p = dict(base_params)
        p["Etot"] = dict(base_params["Etot"])
        p["Etot"][enzyme] = value
        params_list.append(p)

    raw_results = run_params_grid(model, base_params, params_list)
    return [{"etot_uM": value, **r} for value, r in zip(scan_values_uM, raw_results)]


def etot_scan_2d(model, base_params, enzyme_x, enzyme_y, values_x, values_y):
    """雙酵素交叉掃描（two-factor grid）：enzyme_x/enzyme_y 各自的 Etot 網格組合，
    其餘酵素固定在基準值。每個網格點彼此獨立，交給 run_params_grid 平行計算，
    再依 (y, x) reshape 回二維矩陣，方便前端畫 heatmap。"""
    params_list = []
    for y_value in values_y:
        for x_value in values_x:
            p = dict(base_params)
            p["Etot"] = dict(base_params["Etot"])
            p["Etot"][enzyme_x] = x_value
            p["Etot"][enzyme_y] = y_value
            params_list.append(p)

    raw_results = run_params_grid(model, base_params, params_list)

    n_x = len(values_x)
    final_grid, flux_grid = [], []
    for row in range(len(values_y)):
        row_results = raw_results[row * n_x:(row + 1) * n_x]
        final_grid.append([r["final_limonene_mM"] for r in row_results])
        flux_grid.append([r["max_post_iptg_flux"] for r in row_results])

    return {
        "enzyme_x": enzyme_x,
        "enzyme_y": enzyme_y,
        "values_x": values_x,
        "values_y": values_y,
        "final_limonene_mM": final_grid,
        "max_post_iptg_flux": flux_grid,
    }


# ---- 平行運算基礎設施：Etot 掃描／2D 掃描的每一格都是獨立的 dFBA 模擬，---
# ---- 用多個 process 同時算，機器有幾顆核心大致就能快幾倍。 -------------

_worker_model = None  # 每個 worker process 自己的模型物件（process 之間不共用）
_MAX_POOL_WORKERS = max(1, (os.cpu_count() or 2) - 1)


def _init_pool_worker(model_path, medium_defaults):
    """ProcessPoolExecutor 的 initializer：每個 worker process 啟動時只執行「一次」，
    在這裡把模型讀進來、套用培養基設定，之後這個 worker 收到的每一個任務都重複使用
    同一個模型物件（用 with model: 包住做暫時性修改）。這樣平行的成本只有『開頭讀一次
    模型檔』，不會變成『每一格都重新讀一次模型檔』（那樣反而更慢，得不償失）。"""
    global _worker_model
    model = build_model(model_path)
    model = apply_medium009(model, medium_defaults)
    model._biomass_rxn = detect_biomass_rxn(model)
    _worker_model = model


def _pool_dfba_task(params):
    """在 worker process 裡跑「一組」完整 dFBA 參數，回傳這組的摘要結果（含完整
    時間序列，供前端畫軌跡圖）。"""
    global _worker_model
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
    """平行執行多組獨立的 dFBA 參數。
    只有 1 組、或機器偵測不到多核心時，直接在目前 process 裡照原本方式跑，
    避免產生行程池的額外開銷（讀模型檔、啟動 process）反而比循序執行更慢。"""
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
ASSETS_DIR = os.path.join(BASE_DIR, "assets")

# static_folder 只指到 assets/ 這個子資料夾（logo、日夜切換圖示這類靜態圖檔），
# 不是整個 BASE_DIR——避免把 main.py／LIM009_params.csv／models/ 這些不該公開
# 的檔案又意外對外開放。放圖檔時記得放進 main.py 同層的 assets/ 資料夾。
app = Flask(__name__, static_folder=ASSETS_DIR, static_url_path="/assets")

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
                         "message": "生長最佳化不可行，請檢查培養基或 Etot 設定"})

    mu_max = growth_sol.fluxes[biomass_rxn]
    model.reactions.get_by_id(biomass_rxn).lower_bound = params["min_growth_frac"] * mu_max
    model.objective = "EX_limonene_e"
    # 對照 LimoneneCOBRA001.m 的 optimizeCbModel(model,'max','one')：'one' 代表
    # 在所有跟最佳解一樣好的解裡（FBA 常見退化解），挑「通量總和（L1 範數）最小」
    # 的那一組——這正是 pFBA 在做的事。單純 model.optimize() 沒有這個篩選規則，
    # solver 撿到哪組算哪組，會導致目標值（限烯烴產量）一樣，但個別反應通量
    # （進而 ATP/NADPH 這類 cofactor 加總）跟參考結果對不起來。改用 pFBA 對齊。
    try:
        prod_sol = pfba(model)
    except Exception:
        # pFBA 在極少數邊界情況下可能不可行（例如目標通量本身就是 0），
        # 退回普通 FBA，不讓整個請求失敗。
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
    """單次動態 FBA 模擬，回傳完整時間序列。"""
    body = request.get_json(force=True, silent=True) or {}
    params = merge_params(body.get("params"))
    model = working_model()
    apply_adjustable_medium(model, params)
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


@app.route("/api/optimize-etot", methods=["POST"])
def optimize_etot_endpoint():
    """用 Nelder-Mead 自動找出讓限烯烴通量最大化的 Etot(DXS/IDI/GPPS/LS) 組合，
    取代手動窮舉網格掃描。以目前 Etot 欄位的值當作搜尋起點。"""
    body = request.get_json(force=True, silent=True) or {}
    params = merge_params(body.get("params"))
    # 4 維同時搜尋（DXS/IDI/GPPS/LS）實測需要 300~500 次評估才會穩定收斂，
    # 60 次左右常常還沒收斂就被截斷——預設值調高，範圍上限也跟著放寬。
    max_iter = int(body.get("max_iter", 200))
    max_iter = max(50, min(max_iter, 600))

    result = optimize_etot_nelder_mead(params, params["Etot"], max_iter=max_iter)
    return jsonify(result)


if __name__ == "__main__":
    app.run(port=5001, debug=True)