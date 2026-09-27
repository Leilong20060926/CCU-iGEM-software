function bindLangToggle() {
  document.querySelectorAll(".lang-toggle__btn").forEach((btn) => {
    if (!btn.dataset.lang) return; // Skip the simple/advanced toggle button -- that isn't a language switch
    btn.addEventListener("click", () => {
      if (btn.dataset.lang === state.lang) return;
      state.lang = btn.dataset.lang;
      applyStaticI18n();
      updateModeUI();
      // If the status bar currently shows ready/checking, re-render it in the
      // new language; keep error messages in their original language to avoid confusion.
      if (els.modelStatus.querySelector(".dot--ok")) setStatus("ok", t("ready"));
      else if (els.modelStatus.querySelector(".dot--pending")) setStatus("pending", t("checking"));
      // The preset dropdown's "-- Select a saved preset --" option is built
      // dynamically via innerHTML, not a static data-i18n element, so
      // applyStaticI18n() can't reach it -- it needs a separate redraw on
      // language switch. Remember which preset is currently selected and
      // restore that selection after the redraw.
      const selectedPreset = els.presetSelect.value;
      refreshPresetSelect();
      els.presetSelect.value = selectedPreset;
      rerenderAllResults();
    });
  });
}

/* ================= Day/night theme toggle ================= */
const THEME_STORAGE_KEY = "limoneneCobra.theme.v1";

function initTheme() {
  let saved = null;
  try { saved = localStorage.getItem(THEME_STORAGE_KEY); } catch (err) { /* ignore read failure */ }
  const theme = saved === "light" ? "light" : "dark";
  document.documentElement.dataset.theme = theme;
  updateThemeIcon(theme);
}

function bindThemeToggle() {
  document.querySelectorAll(".theme-toggle").forEach((btn) => {
    btn.addEventListener("click", () => {
      const next = document.documentElement.dataset.theme === "light" ? "dark" : "light";
      document.documentElement.dataset.theme = next;
      updateThemeIcon(next);
      try { localStorage.setItem(THEME_STORAGE_KEY, next); } catch (err) { /* ignore write failure */ }
      rerenderAllResults();
    });
  });
}

/* ================= Simple/advanced mode toggle ================= */
function bindUiModeToggle() {
  document.querySelectorAll("#uiModeToggle .lang-toggle__btn").forEach((btn) => {
    btn.addEventListener("click", () => {
      if (btn.dataset.uimode === state.uiMode) return;
      document.querySelectorAll("#uiModeToggle .lang-toggle__btn").forEach((b) => {
        b.classList.toggle("is-active", b === btn);
      });
      state.uiMode = btn.dataset.uimode;
      updateUiMode();
    });
  });
}

function updateUiMode() {
  const isSimple = state.uiMode === "simple";
  // Simple mode: keeps only the quick-estimate sub-mode under the "Etot
  // Analysis" tab, hiding the mode-tab row, advanced medium settings, and
  // the analysis-method dropdown -- COBRA jargon and advanced options never
  // appear on screen; the user just fills in Etot and clicks to get a number.
  els.modeTabs.hidden = isSimple;
  if (els.mediumFieldGroup) els.mediumFieldGroup.hidden = isSimple;
  if (els.scanSubModeField) els.scanSubModeField.hidden = isSimple;
  // Plain-language explanations only show in simple mode; advanced users
  // already know COBRA terminology and don't need these hints.
  if (els.simpleIntroPanel) els.simpleIntroPanel.hidden = !isSimple;
  if (els.etotSimpleHint) els.etotSimpleHint.hidden = !isSimple;
  if (els.surrogateSimpleHint) els.surrogateSimpleHint.hidden = !isSimple;
  if (isSimple) {
    state.mode = "scan";
    state.scanSubMode = "surrogate";
    const scanSubModeSelect = document.getElementById("scanSubMode");
    if (scanSubModeSelect) scanSubModeSelect.value = "surrogate";
  }
  // When switching back to advanced mode, sync the tab row's "current tab"
  // highlight back to state.mode, so it doesn't keep showing whichever tab
  // was stale from before switching to simple mode.
  els.modeTabs.querySelectorAll(".mode-tab").forEach((btn) => {
    btn.classList.toggle("is-active", btn.dataset.mode === state.mode);
  });
  updateModeUI();
}

function updateThemeIcon(theme) {
  // The dark theme shows a sun icon (click to switch to day), the light
  // theme shows a moon icon (click to switch to night). Icon files live in
  // assets/sun.svg, assets/moon.svg (the assets/ folder next to main.py).
  // The homepage and both sub-pages (Related Info / About Us) each have
  // their own toggle button; update them all together.
  const src = theme === "light" ? "/assets/moon.svg" : "/assets/sun.svg";
  document.getElementById("themeToggleIcon").src = src;
  document.querySelectorAll(".theme-toggle-icon-mirror").forEach((img) => { img.src = src; });
}

/* ================= Page navigation (console / related info / about us) ================= */
function bindInfoPages() {
  const mainView = document.getElementById("mainView");
  const infoPage = document.getElementById("infoPage");
  const aboutPage = document.getElementById("aboutPage");

  const showPage = (page) => {
    mainView.hidden = true;
    infoPage.hidden = true;
    aboutPage.hidden = true;
    page.hidden = false;
    window.scrollTo(0, 0);
  };
  const showMain = () => {
    infoPage.hidden = true;
    aboutPage.hidden = true;
    mainView.hidden = false;
    window.scrollTo(0, 0);
  };

  document.getElementById("navInfoBtn").addEventListener("click", () => showPage(infoPage));
  document.getElementById("navAboutBtn").addEventListener("click", () => showPage(aboutPage));
  document.getElementById("infoBackBtn").addEventListener("click", showMain);
  document.getElementById("aboutBackBtn").addEventListener("click", showMain);
}

/* Charts, plain-language summaries, and the technical report all generate
   their text and colors dynamically via t() / cssVar() rather than static
   data-i18n attributes, so applyStaticI18n() never touches them. After a
   language or theme switch, redraw them from the cached most-recent result
   instead of re-hitting the API. Any mode with no cached result is skipped. */
function rerenderAllResults() {
  if (state.lastResults.fba) renderFbaAll(state.lastResults.fba);
  if (state.dfbaHistory.length) renderDfbaCharts();
  if (state.scanDisplayMode === "all" && state.lastResults.scanAll) {
    renderScanAllCharts(state.lastResults.scanAll);
  } else if (state.lastResults.scan) {
    renderScanCharts(state.lastResults.scan);
  }
  if (state.lastResults.scan2d) renderScan2dCharts(state.lastResults.scan2d);
  if (state.lastResults.optimize) renderOptimizeResults(state.lastResults.optimize);
  if (state.lastResults.surrogate) {
    renderSurrogateResults(state.lastResults.surrogate, state.lastResults.surrogateCurve, state.lastResults.surrogateLS);
  }
}

function setStatus(kind, text) {
  const dotClass = kind === "ok" ? "dot--ok" : kind === "err" ? "dot--err" : "dot--pending";
  els.modelStatus.innerHTML = `<span class="dot ${dotClass}"></span> ${text}`;
}

function populateDefaults(d) {
  for (const enzyme of ["DXS", "IDI", "GPPS", "LS"]) {
    document.getElementById(`etot-${enzyme}`).value = d.Etot[enzyme];
  }
  document.getElementById("iptgTime").value = d.iptg_start_time_h;
  document.getElementById("minGrowthFrac").value = d.min_growth_frac;
  document.getElementById("totalTime").value = (d.time_step_h * d.n_steps).toFixed(2);
  document.getElementById("scanValues").value = d.scan_values["LS"].join(", ");

  document.getElementById("glycerolUptake").value = Math.abs(d["maxUptake.EX_glyc_e"]);
  document.getElementById("oxygenUptake").value = Math.abs(d["maxUptake.EX_o2_e"]);
  document.getElementById("aaUptake").value = Math.abs(d["sim.richAminoAcidMaxUptake"]);
  document.getElementById("atpmLb").value = d.atpm_lb;
  document.getElementById("dxpsUb").value = d.dxps_ub;
  document.getElementById("metatUb").value = d.metat_ub;
  document.getElementById("fppsUb").value = d.fpps_ub;

  document.getElementById("scanEnzyme").addEventListener("change", (e) => {
    document.getElementById("scanValues").value = d.scan_values[e.target.value].join(", ");
    updateScanEnzymeDisabledStates();
  });

  document.getElementById("scanValuesX").value = d.scan_values["LS"].join(", ");
  document.getElementById("scanValuesY").value = d.scan_values["DXS"].join(", ");
  document.getElementById("scanValuesZ").value = d.scan_values["IDI"].join(", ");
  document.getElementById("scanValuesW").value = d.scan_values["GPPS"].join(", ");

  for (const axis of ["X", "Y", "Z", "W"]) {
    document.getElementById(`scanEnzyme${axis}`).addEventListener("change", (e) => {
      document.getElementById(`scanValues${axis}`).value = d.scan_values[e.target.value].join(", ");
      updateScanEnzymeDisabledStates();
    });
  }

  document.getElementById("scanMultiCount").addEventListener("change", () => {
    updateScanMultiFieldsVisibility();
    updateScanEnzymeDisabledStates();
    updateModeUI();
  });
  document.getElementById("scanStrategy").addEventListener("change", updateScanStrategyUI);
  document.getElementById("scanSubMode").addEventListener("change", (e) => {
    state.scanSubMode = e.target.value;
    updateModeUI();
  });

  updateScanMultiFieldsVisibility();
  updateScanStrategyUI();
  updateScanEnzymeDisabledStates();
}

/* Whichever enzyme is currently selected for scanning (single- or
   multi-enzyme mode) has no effect if you also fill in its value in the
   "total enzyme concentration Etot" field -- the scan sequence overrides it
   -- which used to confuse users into thinking "nothing I change does
   anything." This grays out and disables the corresponding Etot input and
   adds a hint next to its label, so it's obvious at a glance that the
   scan setting has taken over that enzyme. */
function updateScanEnzymeDisabledStates() {
  const activeEnzymes = new Set();
  // Locking fields only makes sense under the "manual scan" sub-mode
  // (filling them in gets overwritten by the scan sequence anyway); the
  // "auto-optimize" and "quick estimate" sub-modes both consume all four
  // fields at once, with no notion of "which one is currently being
  // scanned." This used to not check scanSubMode, so switching to either of
  // those two sub-modes would misread the scanEnzyme dropdown's leftover
  // value (LS by default) as "currently being scanned" and lock it.
  if (state.mode === "scan" && state.scanSubMode === "manual") {
    activeEnzymes.add(document.getElementById("scanEnzyme").value);
  } else if (state.mode === "scan2d") {
    const count = Number(document.getElementById("scanMultiCount").value) || 2;
    const axes = ["X", "Y", "Z", "W"].slice(0, count);
    axes.forEach((axis) => activeEnzymes.add(document.getElementById(`scanEnzyme${axis}`).value));
  }
  for (const enzyme of ["DXS", "IDI", "GPPS", "LS"]) {
    const input = document.getElementById(`etot-${enzyme}`);
    const hintId = `etotDisabledHint-${enzyme}`;
    let hint = document.getElementById(hintId);
    if (activeEnzymes.has(enzyme)) {
      input.disabled = true;
      if (!hint) {
        hint = document.createElement("span");
        hint.id = hintId;
        hint.className = "hint";
        hint.style.display = "block";
        hint.textContent = t("scanEnzymeDisabledHint");
        input.insertAdjacentElement("afterend", hint);
      }
    } else {
      input.disabled = false;
      if (hint) hint.remove();
    }
  }
}

/* The number of scanned enzymes (2/3/4) decides whether the 3rd and 4th enzyme fields are shown. */
function updateScanMultiFieldsVisibility() {
  const count = Number(document.getElementById("scanMultiCount").value) || 2;
  els.scanEnzymeZField.hidden = count < 3;
  els.scanValuesZField.hidden = count < 3;
  els.scanEnzymeWField.hidden = count < 4;
  els.scanValuesWField.hidden = count < 4;
}

/* The scan strategy (exhaustive grid / LHS / Sobol) decides whether the
   "sample count" field is shown and swaps the hint text, and relabels the
   "value sequence" input fields to remind the user: in non-grid mode, this
   field is read as a [min, max] range, not an enumerated list of candidate
   values. */
function updateScanStrategyUI() {
  const strategy = document.getElementById("scanStrategy").value;
  const isSampled = strategy !== "grid";
  els.scanNSamplesField.hidden = !isSampled;
  document.getElementById("scanMultiGridHint").hidden = isSampled;

  const hintKey = strategy === "lhs" ? "scanStrategyLhsHint"
    : strategy === "sobol" ? "scanStrategySobolHint" : "scanStrategyGridHint";
  document.getElementById("scanStrategyHint").textContent = t(hintKey);

  for (const axis of ["X", "Y", "Z", "W"]) {
    const labelEl = document.getElementById(`scanValues${axis}LabelEl`);
    if (!labelEl) continue;
    const baseKey = { X: "scanValuesXLabel", Y: "scanValuesYLabel", Z: "scanValuesZLabel", W: "scanValuesWLabel" }[axis];
    labelEl.textContent = isSampled ? `${t(baseKey)}${t("scanValueRangeSuffix")}` : t(baseKey);
  }
}

function bindModeTabs() {
  els.modeTabs.querySelectorAll(".mode-tab").forEach((btn) => {
    btn.addEventListener("click", () => {
      els.modeTabs.querySelectorAll(".mode-tab").forEach((b) => b.classList.remove("is-active"));
      btn.classList.add("is-active");
      state.mode = btn.dataset.mode;
      updateModeUI();
    });
  });
}

function updateModeUI() {
  const cfg = modeConfig()[state.mode];
  els.runBtnLabel.textContent = cfg.label;

  const isScan = state.mode === "scan";
  const isOptimize = isScan && state.scanSubMode === "optimize";
  const isSurrogate = isScan && state.scanSubMode === "surrogate";
  els.scanFieldGroup.hidden = !isScan;
  els.scan2dFieldGroup.hidden = state.mode !== "scan2d";
  els.scanManualFields.hidden = isScan && (isOptimize || isSurrogate);
  els.scanOptimizeFields.hidden = !isOptimize;
  els.scanSurrogateFields.hidden = !isSurrogate;
  // The quick-estimate sub-mode is a pure math equation (only consuming
  // Etot), completely unrelated to parameters like induction time or the
  // growth-rate floor that only dFBA/FBA use -- the whole field group is
  // hidden, so users don't mistakenly think changing those numbers affects
  // the quick estimate's result.
  els.inductionFieldGroup.hidden = isSurrogate;
  els.inducedField.hidden = state.mode !== "fba";
  els.durationField.hidden = state.mode === "fba" || isOptimize || isSurrogate;
  // step_mode only makes sense for modes that actually run the dFBA engine
  // internally: dFBA itself, and Etot scan / 2D scan (each grid point
  // internally runs its own dFBA too); the auto-optimize and quick-estimate
  // sub-modes are not dFBA under the hood (they're steady-state FBA and the
  // external surrogate model, respectively), so both are excluded here.
  els.stepModeField.hidden = isOptimize || isSurrogate
    || !["dfba", "scan", "scan2d"].includes(state.mode);
  els.compareField.hidden = state.mode !== "dfba";
  els.compareClearRow.hidden = state.mode !== "dfba" || state.dfbaHistory.length === 0;

  for (const key of ["fbaResults", "dfbaResults", "scanResults", "scan2dResults", "optimizeResults", "surrogateResults"]) {
    document.getElementById(key).hidden = key !== cfg.resultBlock;
  }
  if (state.mode === "scan2d") updateScanMultiFieldsVisibility();
  updateScanEnzymeDisabledStates();
  els.statusHint.textContent = "";
  els.statusHint.className = "hint";
}

function readEtot() {
  return {
    DXS: Number(document.getElementById("etot-DXS").value),
    IDI: Number(document.getElementById("etot-IDI").value),
    GPPS: Number(document.getElementById("etot-GPPS").value),
    LS: Number(document.getElementById("etot-LS").value),
  };
}

function readCommonParams() {
  const totalTime = Number(document.getElementById("totalTime").value);
  const timeStep = state.defaults ? state.defaults.time_step_h : 0.25;
  return {
    Etot: readEtot(),
    iptg_start_time_h: Number(document.getElementById("iptgTime").value),
    min_growth_frac: Number(document.getElementById("minGrowthFrac").value),
    n_steps: Math.max(1, Math.round(totalTime / timeStep)),
    "maxUptake.EX_glyc_e": Number(document.getElementById("glycerolUptake").value),
    "maxUptake.EX_o2_e": Number(document.getElementById("oxygenUptake").value),
    "sim.richAminoAcidMaxUptake": Number(document.getElementById("aaUptake").value),
    atpm_lb: Number(document.getElementById("atpmLb").value),
    dxps_ub: Number(document.getElementById("dxpsUb").value),
    metat_ub: Number(document.getElementById("metatUb").value),
    fpps_ub: Number(document.getElementById("fppsUb").value),
    step_mode: document.getElementById("stepMode").value,
  };
}