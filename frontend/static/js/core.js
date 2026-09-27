/* ================= state / elements ================= */
const state = {
  mode: "fba",
  lang: "en",
  uiMode: "simple",   // "simple" (general users, only the quick estimate) / "advanced" (COBRA users, full feature set)
  defaults: null,
  compare: false,
  dfbaHistory: [],       // { runNo, time_h, limonene_mM, biomass_gDW_L, glycerol_mM, growth_rate_h, limonene_flux }
  lastResults: {},        // Most recent result per mode, used for CSV export: { fba, dfba, scan, scan2d }
  scan2dSliceIndices: [],  // For multi-enzyme scans (3-4 enzymes), the currently selected slice index for dimensions 3 and 4
  scanDisplayMode: "single",  // "single" (one enzyme) or "all" (all four enzymes overlaid), decides
                               // which data to use when redrawing on language/theme switch or CSV export
  scanSubMode: "surrogate",  // "manual" (manual scan) / "optimize" (auto-optimize) /
                             // "surrogate" (quick estimate, the default, since it's the
                             // fastest thing a user sees when switching to "Etot Analysis" --
                             // millisecond-scale instant response) --
                           // the two sub-modes under the "Etot Analysis" tab
};

const els = {
  modeTabs: document.getElementById("modeTabs"),
  mediumFieldGroup: document.getElementById("mediumFieldGroup"),
  scanSubModeField: document.getElementById("scanSubModeField"),
  simpleIntroPanel: document.getElementById("simpleIntroPanel"),
  etotSimpleHint: document.getElementById("etotSimpleHint"),
  surrogateSimpleHint: document.getElementById("surrogateSimpleHint"),
  modelStatus: document.getElementById("modelStatus"),
  runBtn: document.getElementById("runBtn"),
  runBtnLabel: document.getElementById("runBtnLabel"),
  statusHint: document.getElementById("statusHint"),
  scanFieldGroup: document.getElementById("scanFieldGroup"),
  scan2dFieldGroup: document.getElementById("scan2dFieldGroup"),
  scanSurrogateFields: document.getElementById("scanSurrogateFields"),
  inductionFieldGroup: document.getElementById("inductionFieldGroup"),
  scanNSamplesField: document.getElementById("scanNSamplesField"),
  scanEnzymeZField: document.getElementById("scanEnzymeZField"),
  scanValuesZField: document.getElementById("scanValuesZField"),
  scanEnzymeWField: document.getElementById("scanEnzymeWField"),
  scanValuesWField: document.getElementById("scanValuesWField"),
  scan2dSliceControls: document.getElementById("scan2dSliceControls"),
  scanManualFields: document.getElementById("scanManualFields"),
  scanOptimizeFields: document.getElementById("scanOptimizeFields"),
  stepModeField: document.getElementById("stepModeField"),
  inducedField: document.getElementById("inducedField"),
  durationField: document.getElementById("durationField"),
  fbaResults: document.getElementById("fbaResults"),
  dfbaResults: document.getElementById("dfbaResults"),
  scanResults: document.getElementById("scanResults"),
  scan2dResults: document.getElementById("scan2dResults"),
  compareField: document.getElementById("compareField"),
  compareToggle: document.getElementById("compareToggle"),
  compareClearRow: document.getElementById("compareClearRow"),
  presetSelect: document.getElementById("presetSelect"),
  presetName: document.getElementById("presetName"),
};

function modeConfig() {
  const scanResultBlocks = { manual: "scanResults", optimize: "optimizeResults", surrogate: "surrogateResults" };
  const scanLabels = { manual: t("runScan"), optimize: t("runOptimize"), surrogate: t("runSurrogate") };
  return {
    fba: { label: t("runFba"), resultBlock: "fbaResults" },
    dfba: { label: t("runDfba"), resultBlock: "dfbaResults" },
    // scan now covers three sub-modes -- "manual scan," "auto-optimize," and
    // "quick estimate" -- which result block and run-button label to use
    // depends on state.scanSubMode's current value.
    scan: {
      label: scanLabels[state.scanSubMode] || scanLabels.manual,
      resultBlock: scanResultBlocks[state.scanSubMode] || scanResultBlocks.manual,
    },
    scan2d: { label: t("runScan2dN", document.getElementById("scanMultiCount").value), resultBlock: "scan2dResults" },
  };
}

/* Reads the current theme's (day/night) CSS variable value, so Plotly chart
   colors follow the theme. Called fresh before every draw instead of being
   cached as a constant, so a redraw after a theme switch picks up the new
   colors. */
function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

function plotlyLayoutBase() {
  return {
    paper_bgcolor: "transparent",
    plot_bgcolor: "transparent",
    font: { family: "JetBrains Mono, monospace", color: cssVar("--text-muted"), size: 11 },
    margin: { l: 55, r: 20, t: 40, b: 45 },
    legend: { orientation: "h", y: -0.2 },
    xaxis: { gridcolor: cssVar("--border-soft"), zerolinecolor: cssVar("--border") },
    yaxis: { gridcolor: cssVar("--border-soft"), zerolinecolor: cssVar("--border") },
  };
}

function chartTitleColor() {
  return cssVar("--text");
}

function accentColors() {
  // All five colors come from the NoFold logo: gold / olive green / light
  // blue / grain orange-brown / bright gold. The gold tones (colors 1 and 5)
  // use --accent-bright / --accent2-bright, which don't darken under the
  // day theme, so charts keep the logo's original vivid gold in both
  // day and night mode.
  return [
    cssVar("--accent-bright"),
    cssVar("--leaf"),
    cssVar("--steel"),
    cssVar("--danger"),
    cssVar("--accent2-bright"),
  ];
}

/* Fixed color-cycling rule: gold / green / blue / red / purple, in that
   order (the order accentColors() returns). Any bar chart with "several
   side-by-side categories in one chart" (e.g. one bar per enzyme, one bar
   per reaction) should pull colors from this function in sequence -- never
   use a single color for the whole chart, since that reads as "all one
   color" and makes it impossible to tell which bar is which category. A
   single-value metric (a gauge, or one bar comparing two numbers) instead
   assigns a semantic color directly (see each chart's own comments). */
function categoricalColors(n) {
  const palette = accentColors();
  return Array.from({ length: n }, (_, i) => palette[i % palette.length]);
}

// NOTE: init() is defined here but deliberately NOT called here. It calls
// functions defined in later-loaded files (initTheme, bindLangToggle, etc.
// in ui-bindings.js), and separate <script src> files only make their
// functions available globally once that file has actually finished
// executing -- unlike a single inline script, where hoisting makes every
// function available from the start. So the call is placed at the very end
// of presets-export-utils.js (the last file loaded), after every other
// script has already run.
async function init() {
  initTheme();
  applyStaticI18n();
  bindLangToggle();
  bindThemeToggle();
  bindUiModeToggle();
  bindInfoPages();
  bindModeTabs();
  bindPresetControls();
  bindCompareControls();
  bindExportControls();
  els.runBtn.addEventListener("click", runCurrentMode);
  updateUiMode();
  refreshPresetSelect();
  try {
    const res = await fetch(`${API_BASE}/api/defaults`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    state.defaults = await res.json();
    populateDefaults(state.defaults);
    applyUrlParams();
    setStatus("ok", t("ready"));
  } catch (err) {
    setStatus("err", t("connectError", err.message));
  }
}