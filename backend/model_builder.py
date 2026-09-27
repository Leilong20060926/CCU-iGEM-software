"""Model construction (formerly the standalone model_builder.py).

Ports the first half of LimoneneCOBRA009.m's prepareModel009().
"""

import os

import cobra
from cobra import Reaction, Metabolite

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

# Mirrors the inorganic_rxns list in LimoneneCOBRA008_Read_V9.m. This used to
# also open EX_nh4_e (ammonium) by mistake -- V9's medium setup does not
# include that; it relies solely on the 20 amino acids for organic nitrogen.
# That stray ammonium opening, combined with the wrong biomass reaction
# below, together caused our tool's numbers to disagree with the external
# V9 MATLAB results.
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

# This ratio table (relative share of each amino acid) matches the earlier
# LimoneneCOBRA001.m, whose logic split a single total uptake across amino
# acids by ratio. V9 instead gives every amino acid the same flat uptake
# cap (see apply_amino_acid_ratios()), so this table is no longer used --
# it's kept only as a historical record, in case the 001.m-style logic is
# ever needed again.
_AA_RATIOS_RAW_LEGACY_001 = [
    0.03, 0.03, 0.04, 0.07,
    0.005, 0.05, 0.15, 0.03,
    0.02, 0.05, 0.10, 0.08,
    0.025, 0.05, 0.10, 0.06,
    0.04, 0.01, 0.03, 0.07,
]

MEDIUM_NAME = "BioShop Terrific Broth (Tryptone 12g/L, Yeast Extract 24g/L)"

# Mirrors LimoneneCOBRA008_Read_V9.m: this uses the auto-detected _core_
# reaction rather than a hardcoded _WT_. Reason: cross-checked by reloading
# the model directly -- under V9's simplified medium (see apply_medium009
# below), _WT_'s growth rate is 0 (it's missing two precursors,
# adenosylcobalamin and colipa_e, which this simplified medium doesn't
# supply), while _core_ is the version that actually grows under this
# medium and is what V9 actually uses.
PREFERRED_BIOMASS_RXN = os.environ.get(
    "LIMONENE_BIOMASS_RXN", "BIOMASS_Ec_iJO1366_core_53p95M")


def detect_biomass_rxn(model):
    """Prefer PREFERRED_BIOMASS_RXN (aligned with the reaction hardcoded in
    LimoneneCOBRA001.m); fall back to "the first reaction with a nonzero
    objective coefficient" if the model doesn't contain that ID."""
    if PREFERRED_BIOMASS_RXN in model.reactions:
        return PREFERRED_BIOMASS_RXN
    biomass_rxns = [r.id for r in model.reactions if r.objective_coefficient != 0]
    return biomass_rxns[0] if biomass_rxns else None


# Reaction IDs whose native-pathway bounds are user-adjustable.
# MEP_PATHWAY_RXNS mirrors LimoneneCOBRA001.m's mep_rxns: we used to only
# apply the "MEP pathway upper bound" parameter to DXPS, but MEPCT / MECDPS /
# IPDPI also need the same bound applied -- opening only DXPS leaves the
# rest of the MEP pathway bottlenecked, and flux downstream stays capped.
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
    """Pick the loader by file extension. .json (cobrapy's native format) is
    the most robust and recommended; .xml/.sbml (SBML) and .mat (COBRA
    Toolbox) are also supported."""
    ext = os.path.splitext(model_path)[1].lower()
    if ext == ".json":
        return cobra.io.load_json_model(model_path)
    if ext in (".xml", ".sbml"):
        return cobra.io.read_sbml_model(model_path)
    if ext == ".mat":
        return cobra.io.load_matlab_model(model_path)
    raise ValueError(f"Unsupported model file format: {ext} (use .json / .xml / .mat)")


def build_model(model_path):
    """Load the raw strain model and add the heterologous limonene pathway
    (mirrors the first half of prepareModel009())."""
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

    # All T7/heterologous reactions start closed; each simulation reopens
    # them according to that step's Etot capacity.
    for rxn_id in ALL_T7_RXNS:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).bounds = (0, 0)

    _use_fast_solver(model)
    return model


def _use_fast_solver(model):
    """Try switching to a much faster LP solver (HiGHS) than the GLPK
    default.

    Enabled by default -- no environment variable needed. It used to
    require opting in manually because, at the time, we were still doing
    bit-for-bit validation against MATLAB/GLPK: a different solver can pick
    a different vertex when the LP solution isn't unique (a degenerate
    solution, which is common under pFBA), causing individual reaction
    fluxes to diverge from the GLPK reference results. Now that this
    cross-validation is complete, it defaults to enabled. If you need
    bit-for-bit alignment with GLPK again later, set the environment
    variable LIMONENE_FAST_SOLVER=0 to disable it."""
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
    """Applies the TB-style medium settings from prepareModel009() (BioShop
    Terrific Broth Formulation)."""
    for rxn in model.exchanges:
        if rxn.lower_bound < 0:
            rxn.lower_bound = 0

    # Explicitly block glucose uptake (the main carbon source is glycerol).
    if "EX_glc__D_e" in model.reactions:
        model.reactions.EX_glc__D_e.lower_bound = 0

    # Cap glycerol uptake (the main carbon source); value comes from
    # LIM009_params.csv's maxUptake.EX_glyc_e.
    model.reactions.EX_glyc_e.lower_bound = -abs(defaults["maxUptake.EX_glyc_e"])
    model.reactions.EX_glyc_e.upper_bound = 1000

    for rxn_id in RICH_MEDIUM_OPEN:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = -1000
    for rxn_id in TRACE_ION_RXNS:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = -1000

    # Aerobic respiration cap; value comes from LIM009_params.csv's
    # maxUptake.EX_o2_e.
    model.reactions.EX_o2_e.lower_bound = -abs(defaults["maxUptake.EX_o2_e"])
    model.reactions.EX_o2_e.upper_bound = 1000

    apply_amino_acid_ratios(model, defaults["sim.richAminoAcidMaxUptake"])

    # There used to be a "micro-allowance" rule here mirroring
    # LimoneneCOBRA001.m (any exchange reaction not explicitly opened gets a
    # tiny -0.01 uptake lower bound). After re-checking
    # LimoneneCOBRA008_Read_V9.m (the medium-setup script V9 actually uses),
    # that rule simply isn't there -- and removing it made our computed
    # mu_max (2.3744851841360686) match the value computed from the real
    # modelDyn.mat (2.3744851841360806, differing only in the 13th decimal
    # place, i.e. floating-point noise). Keeping that rule made the numbers
    # disagree, so it was removed.

    return model


def apply_amino_acid_ratios(model, per_aa_uptake):
    """Mirrors LimoneneCOBRA008_Read_V9.m's openRichAminoAcidsV9(): all 20
    amino acids each get the same flat uptake upper bound (not a single
    total split by ratio), matching the official note on the
    sim.richAminoAcidMaxUptake row in LIM009_params.csv -- "approximates the
    amino acids supplied by tryptone/yeast extract with a single uptake
    cap." The function name stays apply_amino_acid_ratios to avoid touching
    its other call sites; the AA_RATIOS table above is no longer used
    here."""
    value = abs(per_aa_uptake)
    for rxn_id in AA_EXCHANGES:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).lower_bound = -value


def apply_mep_pathway_bound(model, ub):
    """Mirrors LimoneneCOBRA001.m's mep_rxns: DXPS/MEPCT/MECDPS/IPDPI share
    the same MEP-pathway upper bound; loosening only DXPS leaves the rest of
    the pathway bottlenecked."""
    for rxn_id in MEP_PATHWAY_RXNS:
        if rxn_id in model.reactions:
            model.reactions.get_by_id(rxn_id).upper_bound = ub


def apply_adjustable_medium(model, params):
    """Applies the medium-uptake limits and native-pathway bounds that the
    frontend lets users adjust (applied per-request on the working model /
    cobra context passed in -- never touches the shared base-model cache).
    Covered parameters: glycerol/oxygen uptake caps, total rich-medium amino
    acid uptake, the ATPM maintenance-energy lower bound, and the MEP
    pathway (DXPS/MEPCT/MECDPS/IPDPI) / METAT / FPPS upper bounds (any
    reaction missing from the model is safely skipped)."""
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
    """Computes a metabolite's total 'gross production rate' under the
    current solution: sum, over every reaction, of (stoichiometric
    coefficient x flux) restricted to the reactions where that product is
    net-positive (i.e. actually producing it). Used for ATP / NADPH
    cofactor-balance analysis. Returns None if the model doesn't contain
    this metabolite ID."""
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
    """Assembles the data behind the medium & pathway-bounds report (the
    frontend renders it as text). Reactions not found in the model return
    None, which the frontend shows as N/A rather than erroring out. GPPS /
    LIMS_MS_het read their *actual* upper bound directly off the model
    after the induction rules have been applied (not the theoretical
    kcat x Etot capacity before induction), so the report matches the bound
    that's actually handed to the LP solver -- reading the theoretical
    capacity alone would disagree with reality whenever something like a
    GPPS_FIXED_CAPACITY override is in effect."""
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