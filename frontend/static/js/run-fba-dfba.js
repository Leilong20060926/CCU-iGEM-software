async function runCurrentMode() {
  els.runBtn.disabled = true;
  els.statusHint.className = "hint";
  // "Accelerated by parallel computation" is only true for Etot scan (manual)
  // / 2D scan (each grid point is independent, farmed out to multiple CPU
  // cores via ProcessPoolExecutor); FBA/dFBA only run once, auto-optimize's
  // Nelder-Mead searches sequentially, and the quick estimate doesn't run
  // COBRA at all -- these all used to show the same "parallel-accelerated"
  // hint, which was inaccurate, so it now varies by mode.
  const computingKey = (state.mode === "scan2d" || (state.mode === "scan" && state.scanSubMode === "manual"))
    ? "computingParallel"
    : state.mode === "scan" && state.scanSubMode === "surrogate"
      ? "computingSurrogate"
      : "computingSimple";
  els.statusHint.textContent = t(computingKey);
  try {
    if (state.mode === "fba") await runFba();
    else if (state.mode === "dfba") await runDfba();
    else if (state.mode === "scan") {
      if (state.scanSubMode === "optimize") await runOptimize();
      else if (state.scanSubMode === "surrogate") await runSurrogate();
      else await runScan();
    } else await runScan2d();
    els.statusHint.textContent = t("done");
  } catch (err) {
    els.statusHint.className = "hint hint--error";
    els.statusHint.textContent = t("errorPrefix", err.message);
  } finally {
    els.runBtn.disabled = false;
  }
}

async function postJson(path, body) {
  const res = await fetch(`${API_BASE}${path}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  const data = await res.json();
  if (!res.ok || data.error) throw new Error(data.error || `HTTP ${res.status}`);
  return data;
}

/* ---------------- FBA ---------------- */
async function runFba() {
  const induced = document.getElementById("induced").value === "true";
  const params = readCommonParams();
  const data = await postJson("/api/fba", { params, induced });
  state.lastResults.fba = data;

  if (data.status !== "optimal") {
    els.statusHint.className = "hint hint--error";
    els.statusHint.textContent = data.message || t("solverStatus", data.status);
  }

  renderFbaAll(data);
}

/* Centralizes all FBA result rendering (metric cards, enzyme flux chart,
   plain-language summary, charted report, text report) in one function, so
   "right after running" and "redrawing from cache after a language switch"
   can share the same logic. */
function renderFbaAll(data) {
  document.getElementById("mMuMax").textContent = fmt(data.mu_max, 4) + " h⁻¹";
  document.getElementById("mMu").textContent = fmt(data.growth_rate, 4) + " h⁻¹";
  document.getElementById("mLimFlux").textContent = fmt(data.limonene_flux, 5) + " mmol/gDW/h";

  renderFbaEnzymeChart(data);
  renderFbaReport(data);
  renderFbaSummary(data);
  renderFbaCharts(data);
  renderDataTable("fba", data);
}

function renderFbaEnzymeChart(data) {
  const enzymes = Object.keys(data.enzyme_fluxes || {});
  const values = enzymes.map((k) => data.enzyme_fluxes[k] ?? 0);
  Plotly.newPlot("fbaEnzymeChart", [{
    type: "bar",
    x: enzymes,
    y: values,
    marker: { color: categoricalColors(enzymes.length) },
  }], {
    ...plotlyLayoutBase(),
    title: { text: t("enzymeFluxChart"), font: { size: 13, color: chartTitleColor() } },
  }, { displayModeBar: false, responsive: true });
}

/* ---------------- FBA plain-language summary ---------------- */
function renderFbaSummary(data) {
  const el = document.getElementById("fbaSummaryLine");
  if (data.status !== "optimal" || data.mu_max === null || data.mu_max === undefined) {
    el.textContent = t("summaryInfeasible");
    return;
  }
  const pct = data.mu_max > 0 ? fmt((data.growth_rate / data.mu_max) * 100, 1) : "0.0";
  el.textContent = `${t("summaryGrowthPart", pct)} ${t("summaryYieldPart", fmt(data.limonene_flux, 4))}`;
}

/* ---------------- FBA charted report ---------------- */
function renderFbaCharts(data) {
  const mp = data.medium_pathway_report || {};
  const cof = data.cofactors || {};

  // Growth-rate gauge: the needle = actual growth rate during production, the red line = theoretical max growth rate.
  const muMax = data.mu_max || 0;
  const gaugeMax = Math.max(muMax * 1.15, 0.01);
  Plotly.newPlot("fbaGrowthChart", [{
    type: "indicator",
    mode: "gauge+number",
    value: data.growth_rate ?? 0,
    number: { suffix: " /h", valueformat: ".3f", font: { color: chartTitleColor(), size: 24 } },
    gauge: {
      axis: {
        range: [0, gaugeMax],
        tickcolor: cssVar("--text-muted"),
        tickfont: { color: cssVar("--text-muted"), size: 10 },
        // Explicitly place a tick mark at the exact theoretical-max value,
        // so there's a number right next to the red line to check against,
        // instead of eyeballing alignment.
        tickmode: "array",
        tickvals: [0, muMax / 2, muMax],
        ticktext: ["0", fmt(muMax / 2, 2), fmt(muMax, 3)],
      },
      bar: { color: accentColors()[0] },
      bgcolor: "transparent",
      borderwidth: 1,
      bordercolor: cssVar("--border"),
      threshold: {
        line: { color: accentColors()[3], width: 3 },
        thickness: 0.85,
        value: muMax,
      },
    },
  }], {
    ...plotlyLayoutBase(),
    margin: { l: 30, r: 30, t: 50, b: 20 },
    title: { text: t("growthGaugeTitle"), font: { size: 12.5, color: chartTitleColor() } },
  }, { displayModeBar: false, responsive: true });

  // Medium uptake caps (glycerol/oxygen, absolute value for easy comparison).
  const oxygenAbs = Math.abs(mp.oxygen_bounds?.[0] ?? 0);
  const glycerolAbs = Math.abs(mp.glycerol_bounds?.[0] ?? 0);
  Plotly.newPlot("fbaMediumChart", [{
    type: "bar",
    orientation: "h",
    y: [t("oxygenUptakeLabel"), t("glycerolUptakeLabel")],
    x: [oxygenAbs, glycerolAbs],
    text: [fmt(oxygenAbs, 2), fmt(glycerolAbs, 2)],
    textposition: "outside",
    cliponaxis: false,
    marker: { color: [accentColors()[2], accentColors()[1]] },
  }], {
    ...plotlyLayoutBase(),
    title: { text: t("mediumChartTitle"), font: { size: 13, color: chartTitleColor() } },
    xaxis: {
      ...plotlyLayoutBase().xaxis,
      type: "linear",
      title: "mmol/gDW/h",
      range: [0, Math.max(oxygenAbs, glycerolAbs, 1) * 1.4],
    },
    yaxis: {
      ...plotlyLayoutBase().yaxis,
      type: "category",
      automargin: true,
    },
    margin: { l: 140, r: 60, t: 40, b: 40 },
  }, { displayModeBar: false, responsive: true });

  // Enzyme/pathway bound comparison: values span 3 orders of magnitude
  // (0.003 ~ 30), so a log Y axis is used; a missing reaction (no matching
  // ID in the model, value is null) skips its bar and is labeled N/A under
  // that category instead.
  const pathwayLabels = ["DXPS", "METAT", "FPPS", "GPPS", "LIMS_MS_het"];
  const pathwayKeys = ["dxps_ub", "metat_ub", "fpps_ub", "gpps_ub", "lims_ub"];
  const pathwayValues = pathwayKeys.map((k) => mp[k]);
  const pathwayLogVals = pathwayValues.map((v) => (v === null || v === undefined || v <= 0 ? null : v));
  const pathwayText = pathwayValues.map((v) => (v === null || v === undefined ? "" : fmt(v, 3)));
  const pathwayNaIndices = pathwayValues
    .map((v, i) => (v === null || v === undefined ? i : null))
    .filter((i) => i !== null);

  // The log axis only shows integer-power ticks (1, 10, 100 or 0.1, 0.01...), not minor ticks like 2 or 5, to avoid overcrowding.
  const posVals = pathwayLogVals.filter((v) => v !== null && v > 0);
  let tickvals = [1];
  if (posVals.length) {
    const minExp = Math.floor(Math.log10(Math.min(...posVals)));
    const maxExp = Math.ceil(Math.log10(Math.max(...posVals)));
    tickvals = [];
    for (let e = minExp; e <= maxExp; e++) tickvals.push(Math.pow(10, e));
  }
  const ticktext = tickvals.map((v) => (v >= 1 ? String(v) : v.toString()));

  Plotly.newPlot("fbaPathwayChart", [{
    type: "bar",
    x: pathwayLabels,
    y: pathwayLogVals,
    text: pathwayText,
    textposition: "outside",
    cliponaxis: false,
    marker: { color: categoricalColors(pathwayLabels.length) },
  }], {
    ...plotlyLayoutBase(),
    title: { text: t("pathwayChartTitle"), font: { size: 13, color: chartTitleColor() } },
    margin: { l: 60, r: 20, t: 60, b: 45 },
    xaxis: { ...plotlyLayoutBase().xaxis, type: "category" },
    yaxis: {
      ...plotlyLayoutBase().yaxis,
      type: "log",
      title: t("pathwayAxisTitle"),
      tickmode: "array",
      tickvals,
      ticktext,
    },
    annotations: pathwayNaIndices.map((i) => ({
      x: pathwayLabels[i], y: 0.04, yref: "paper", yanchor: "bottom",
      text: t("reportNA"), showarrow: false,
      font: { color: cssVar("--text-muted"), size: 11 },
    })),
  }, { displayModeBar: false, responsive: true });

  // Intracellular ATP/NADPH production-rate comparison (short labels avoid clipped X-axis text).
  const atpVal = cof.atp_production_mmol_gDW_h ?? 0;
  const nadphVal = cof.nadph_production_mmol_gDW_h ?? 0;
  Plotly.newPlot("fbaCofactorChart", [{
    type: "bar",
    x: [t("cofactorAtpShort"), t("cofactorNadphShort")],
    y: [atpVal, nadphVal],
    text: [fmt(atpVal, 3), fmt(nadphVal, 3)],
    textposition: "outside",
    cliponaxis: false,
    marker: { color: [accentColors()[0], accentColors()[2]] },
  }], {
    ...plotlyLayoutBase(),
    title: { text: t("cofactorChartTitle"), font: { size: 13, color: chartTitleColor() } },
    margin: { l: 55, r: 20, t: 50, b: 45 },
    xaxis: { ...plotlyLayoutBase().xaxis, type: "category", tickangle: 0 },
    yaxis: {
      ...plotlyLayoutBase().yaxis,
      type: "linear",
      title: "mmol/gDW/h",
      range: [0, Math.max(atpVal, nadphVal, 1) * 1.2],
    },
  }, { displayModeBar: false, responsive: true });

  // Safety net: right after a container is un-hidden or before layout has
  // settled, size measurements can be off -- force a resize on the next
  // frame to prevent the axes and bars from misrendering.
  requestAnimationFrame(() => {
    ["fbaGrowthChart", "fbaMediumChart", "fbaPathwayChart", "fbaCofactorChart"].forEach((id) => {
      const el = document.getElementById(id);
      if (el && el.data) Plotly.Plots.resize(el);
    });
  });
}

/* ---------------- FBA text report ---------------- */
function fmtBounds(b) {
  if (!b) return t("reportNA");
  return `[${fmt(b[0], 2)}, ${fmt(b[1], 2)}]`;
}

function fmtOrNA(v, digits) {
  if (v === null || v === undefined || Number.isNaN(v)) return t("reportNA");
  return fmt(v, digits);
}

function renderFbaReport(data) {
  const mp = data.medium_pathway_report || {};
  const cof = data.cofactors || {};
  const lines = [
    t("reportTitle", mp.medium_name || t("reportNA")),
    "",
    t("reportMediumHeader"),
    t("reportMediumLine", mp.medium_name || t("reportNA")),
    `${t("reportGlycerolLine")}: ${fmtBounds(mp.glycerol_bounds)} mmol/gDW/h`,
    `${t("reportOxygenLine")}: ${fmtBounds(mp.oxygen_bounds)} mmol/gDW/h`,
    `${t("reportAtpmLine")}: ${fmtOrNA(mp.atpm_lb, 2)} mmol/gDW/h (Lower Bound)`,
    "",
    t("reportPathwayHeader"),
    `${t("reportDxpsLine")} = ${fmtOrNA(mp.dxps_ub, 2)}`,
    `${t("reportMetatLine")} = ${fmtOrNA(mp.metat_ub, 2)}`,
    `${t("reportFppsLine")} = ${fmtOrNA(mp.fpps_ub, 2)}`,
    `${t("reportGppsLine")} = ${fmtOrNA(mp.gpps_ub, 2)}`,
    `${t("reportLimsLine")} = ${fmtOrNA(mp.lims_ub, 2)}`,
    "",
    t("reportGrowthHeader"),
    `${t("reportMuMaxLine")}: ${fmtOrNA(data.mu_max, 4)} h^-1`,
    `${t("reportMuLine")}: ${fmtOrNA(data.growth_rate, 4)} h^-1`,
    `${t("reportYieldLine")}: ${fmtOrNA(data.limonene_flux, 4)} mmol/gDW/h`,
    "",
    t("reportCofactorHeader"),
    `${t("reportAtpLine")}: ${fmtOrNA(cof.atp_production_mmol_gDW_h, 4)} mmol/gDW/h`,
    `${t("reportNadphLine")}: ${fmtOrNA(cof.nadph_production_mmol_gDW_h, 4)} mmol/gDW/h`,
  ];
  document.getElementById("fbaReportBlock").textContent = lines.join("\n");
}

/* ---------------- dFBA ---------------- */
async function runDfba() {
  const params = readCommonParams();
  const data = await postJson("/api/dfba", { params });
  state.lastResults.dfba = data;

  if (!els.compareToggle.checked) state.dfbaHistory = [];
  state.dfbaHistory.push({ runNo: state.dfbaHistory.length + 1, ...data });
  els.compareClearRow.hidden = state.dfbaHistory.length === 0;

  renderDfbaCharts(data);

  if (data.stopped_reason) {
    els.statusHint.className = "hint hint--error";
    els.statusHint.textContent = t("stoppedEarly", data.stopped_reason);
  }
}

/* Keeps dFBA chart rendering as its own function, so a language switch can
   redraw from state.dfbaHistory without re-hitting the API. latestData is
   only used in single-run (non-overlay) mode, to draw the full three-axis
   curve. */
function renderDfbaCharts(latestData) {
  const data = latestData || state.dfbaHistory[state.dfbaHistory.length - 1];
  if (!data) return;

  if (state.dfbaHistory.length <= 1) {
    // Single run: the original full three-axis curve.
    Plotly.newPlot("dfbaConcChart", [
      { x: data.time_h, y: data.limonene_mM, name: t("limoneneLegend"), line: { color: accentColors()[0], width: 2.5 } },
      { x: data.time_h, y: data.biomass_gDW_L, name: t("biomassLegend"), yaxis: "y2", line: { color: accentColors()[1], width: 2 } },
      { x: data.time_h, y: data.glycerol_mM, name: t("glycerolLegend"), yaxis: "y3", line: { color: accentColors()[2], width: 1.5, dash: "dot" } },
    ], {
      ...plotlyLayoutBase(),
      title: { text: t("concChartTitle"), font: { size: 13, color: chartTitleColor() } },
      // Right margin widened to 60px (default is only 20px): this chart has
      // a titled secondary axis on the right (biomass gDW/L), and 20px
      // isn't enough room for the title text with its unit -- Plotly would
      // be forced to stack the characters vertically into the narrow gap,
      // which reads as garbled text. The 3 legend entries also have fairly
      // long labels (with units), so the legend is pushed further down too,
      // to avoid overlapping the X-axis title.
      margin: { ...plotlyLayoutBase().margin, r: 60, b: 80 },
      legend: { ...plotlyLayoutBase().legend, y: -0.3 },
      xaxis: { ...plotlyLayoutBase().xaxis, title: t("timeAxis") },
      yaxis: { ...plotlyLayoutBase().yaxis, title: t("limoneneAxis") },
      yaxis2: { overlaying: "y", side: "right", showgrid: false, title: t("biomassAxis") },
      yaxis3: { overlaying: "y", side: "right", position: 1, showgrid: false, visible: false },
    }, { displayModeBar: false, responsive: true });

    Plotly.newPlot("dfbaFluxChart", [
      { x: data.time_h, y: data.limonene_flux, name: t("limFluxLegend"), line: { color: accentColors()[0], width: 2 } },
      { x: data.time_h, y: data.growth_rate_h, name: t("growthRateLegend"), yaxis: "y2", line: { color: accentColors()[1], width: 2 } },
    ], {
      ...plotlyLayoutBase(),
      title: { text: t("fluxChartTitle"), font: { size: 13, color: chartTitleColor() } },
      margin: { ...plotlyLayoutBase().margin, r: 50 },
      xaxis: { ...plotlyLayoutBase().xaxis, title: t("timeAxis") },
      yaxis: { ...plotlyLayoutBase().yaxis, title: "mmol/gDW/h" },
      yaxis2: { overlaying: "y", side: "right", showgrid: false, title: "h⁻¹" },
    }, { displayModeBar: false, responsive: true });
  } else {
    // Overlay-comparison mode: multiple runs overlay only the limonene concentration/flux curves, to avoid several dual-Y-axis sets interfering with each other.
    // The legend entry count grows with each run (one more line per run),
    // making it more prone to overlapping the X-axis title than a
    // fixed-entry-count chart -- the bottom margin is likewise widened and
    // the legend pushed further down.
    const legendMargin = { margin: { ...plotlyLayoutBase().margin, b: 90 }, legend: { ...plotlyLayoutBase().legend, y: -0.35 } };
    const concTraces = state.dfbaHistory.map((run, i) => ({
      x: run.time_h, y: run.limonene_mM,
      name: `${t("compareRunLabel", run.runNo)} · ${t("limoneneLegend")}`,
      line: { color: accentColors()[i % accentColors().length], width: 2 },
    }));
    Plotly.newPlot("dfbaConcChart", concTraces, {
      ...plotlyLayoutBase(),
      title: { text: t("concChartTitle"), font: { size: 13, color: chartTitleColor() } },
      ...legendMargin,
      xaxis: { ...plotlyLayoutBase().xaxis, title: t("timeAxis") },
      yaxis: { ...plotlyLayoutBase().yaxis, title: t("limoneneAxis") },
    }, { displayModeBar: false, responsive: true });

    const fluxTraces = state.dfbaHistory.map((run, i) => ({
      x: run.time_h, y: run.limonene_flux,
      name: `${t("compareRunLabel", run.runNo)} · ${t("limFluxLegend")}`,
      line: { color: accentColors()[i % accentColors().length], width: 2 },
    }));
    Plotly.newPlot("dfbaFluxChart", fluxTraces, {
      ...plotlyLayoutBase(),
      title: { text: t("fluxChartTitle"), font: { size: 13, color: chartTitleColor() } },
      ...legendMargin,
      xaxis: { ...plotlyLayoutBase().xaxis, title: t("timeAxis") },
      yaxis: { ...plotlyLayoutBase().yaxis, title: "mmol/gDW/h" },
    }, { displayModeBar: false, responsive: true });
  }

  renderDfbaT7FluxChart(data);
  renderDfbaYieldSummary();
  renderDataTable("dfba", data);
}

/* Each enzyme's (DXS/IDI forward-reverse/GPPS/LS) own flux over time --
   the backend's /api/dfba has always computed t7_fluxes, it just wasn't
   plotted before. In overlay-comparison mode this only draws the most
   recent run (overlaying many runs x 5 enzymes would be too cluttered to
   read); single-run mode behaves the same as the other charts. */
function renderDfbaT7FluxChart(data) {
  const t7 = data.t7_fluxes || {};
  const seriesDefs = [
    ["DXS", "t7FluxDXS"],
    ["IDI_fwd", "t7FluxIDIFwd"],
    ["IDI_rev", "t7FluxIDIRev"],
    ["GPPS", "t7FluxGPPS"],
    ["LS", "t7FluxLS"],
  ];
  const traces = seriesDefs
    .filter(([key]) => t7[key])
    .map(([key, labelKey], i) => ({
      x: data.time_h,
      y: t7[key],
      name: t(labelKey),
      line: { color: accentColors()[i % accentColors().length], width: 2 },
    }));
  Plotly.newPlot("dfbaT7FluxChart", traces, {
    ...plotlyLayoutBase(),
    title: { text: t("t7FluxChartTitle"), font: { size: 13, color: chartTitleColor() } },
    // The 5-line horizontal legend (DXS/IDIx2/GPPS/LS) has fairly long
    // labels, and competing with the X-axis title "Time (h)" for the same
    // vertical space would overlap them -- the bottom margin is widened and
    // the legend pushed further down to keep them apart.
    margin: { ...plotlyLayoutBase().margin, b: 90 },
    legend: { ...plotlyLayoutBase().legend, y: -0.35 },
    xaxis: { ...plotlyLayoutBase().xaxis, title: t("timeAxis") },
    yaxis: { ...plotlyLayoutBase().yaxis, title: "mmol/gDW/h" },
  }, { displayModeBar: false, responsive: true });
}

// Limonene's (C10H16) molecular weight, for converting mM -> g/L: g/L = mM * MW(g/mol) / 1000.
const LIMONENE_MW_G_PER_MOL = 136.24;

/* "Total yield" summary: since this is a batch culture (no continuous
   feed/harvest), total yield is simply the final limonene concentration at
   the end of the simulation. The metric card shows the endpoint number
   from the most recent run (mM + converted to g/L); the bar chart lists
   the endpoint concentration for every run currently in the overlay
   comparison, making it easy to compare at a glance which parameter set
   yields the most -- in single-run mode this is just one bar. */
function renderDfbaYieldSummary() {
  const history = state.dfbaHistory;
  if (!history.length) return;
  const latest = history[history.length - 1];
  const n = latest.time_h.length;
  const finalTiterMM = latest.limonene_mM[n - 1];
  const finalTiterGL = finalTiterMM * LIMONENE_MW_G_PER_MOL / 1000;
  const finalBiomass = latest.biomass_gDW_L[n - 1];
  const totalTime = latest.time_h[n - 1];

  document.getElementById("mDfbaFinalTiter").textContent =
    `${fmt(finalTiterMM, 3)} mM (${fmt(finalTiterGL, 3)} g/L)`;
  document.getElementById("mDfbaFinalBiomass").textContent = `${fmt(finalBiomass, 3)} gDW/L`;
  document.getElementById("mDfbaSimTime").textContent = `${fmt(totalTime, 2)} h`;

  const runLabels = history.map((run) => t("compareRunLabel", run.runNo));
  const finalTiters = history.map((run) => run.limonene_mM[run.limonene_mM.length - 1]);
  const finalTitersGL = finalTiters.map((v) => v * LIMONENE_MW_G_PER_MOL / 1000);

  Plotly.newPlot("dfbaYieldChart", [{
    type: "bar",
    x: runLabels,
    y: finalTiters,
    text: finalTitersGL.map((v, i) => `${fmt(finalTiters[i], 2)} mM<br>${fmt(v, 3)} g/L`),
    textposition: "outside",
    cliponaxis: false,
    marker: { color: categoricalColors(runLabels.length) },
  }], {
    ...plotlyLayoutBase(),
    title: { text: t("dfbaYieldChartTitle"), font: { size: 13, color: chartTitleColor() } },
    margin: { l: 55, r: 20, t: 60, b: 45 },
    xaxis: { ...plotlyLayoutBase().xaxis, type: "category" },
    yaxis: { ...plotlyLayoutBase().yaxis, title: t("limoneneAxis") },
  }, { displayModeBar: false, responsive: true });
}

/* ================= Raw data table ================= */
// The raw data table for each mode uses the same data as exportCsv(), just
// rendered as an HTML table viewable directly on the page, instead of
// having to download a file to see the numbers. The table itself uses
// Number.prototype.toFixed to control digit count, avoiding a flood of
// long floating-point tails.
function buildTableHtml(headers, rows) {
  const thead = `<thead><tr>${headers.map((h) => `<th>${h}</th>`).join("")}</tr></thead>`;
  const tbody = `<tbody>${rows.map((r) =>
    `<tr>${r.map((v) => `<td>${v}</td>`).join("")}</tr>`).join("")}</tbody>`;
  return `<table class="data-table">${thead}${tbody}</table>`;
}

function renderDataTable(mode, data) {
  const containerId = {
    fba: "fbaDataTable", dfba: "dfbaDataTable", scan: "scanDataTable",
    scan2d: "scan2dDataTable", optimize: "optimizeDataTable",
  }[mode];
  const container = document.getElementById(containerId);
  if (!container) return;

  let html = "";
  if (mode === "fba") {
    const rows = Object.entries(data.enzyme_fluxes || {}).map(([k, v]) => [k, fmt(v, 5)]);
    rows.push(["mu_max (h⁻¹)", fmt(data.mu_max, 4)]);
    rows.push(["growth_rate (h⁻¹)", fmt(data.growth_rate, 4)]);
    rows.push(["limonene_flux (mmol/gDW/h)", fmt(data.limonene_flux, 5)]);
    const cof = data.cofactors || {};
    rows.push(["ATP (mmol/gDW/h)", fmt(cof.atp_production_mmol_gDW_h, 4)]);
    rows.push(["NADPH (mmol/gDW/h)", fmt(cof.nadph_production_mmol_gDW_h, 4)]);
    html = buildTableHtml([t("rawDataMetricCol"), t("rawDataValueCol")], rows);
  } else if (mode === "dfba") {
    const t7 = data.t7_fluxes || {};
    const headers = ["t (h)", "phase", "solve mode", "biomass (gDW/L)", "glycerol (mM)", "limonene (mM)",
      "μ (h⁻¹)", "lim.flux", "DXS", "IDI_fwd", "IDI_rev", "GPPS", "LS"];
    const rows = data.time_h.map((tv, i) => [
      fmt(tv, 2), data.phase?.[i] ?? "—", data.mode?.[i] ?? "—",
      fmt(data.biomass_gDW_L[i], 4), fmt(data.glycerol_mM[i], 3),
      fmt(data.limonene_mM[i], 4), fmt(data.growth_rate_h[i], 4), fmt(data.limonene_flux[i], 5),
      fmt(t7.DXS && t7.DXS[i], 5), fmt(t7.IDI_fwd && t7.IDI_fwd[i], 5),
      fmt(t7.IDI_rev && t7.IDI_rev[i], 5), fmt(t7.GPPS && t7.GPPS[i], 5), fmt(t7.LS && t7.LS[i], 5),
    ]);
    html = buildTableHtml(headers, rows);
  } else if (mode === "scan") {
    const headers = ["Etot (µM)", "Etot factor", "final limonene (mM)", "max post-IPTG flux (mmol/gDW/h)",
      "mean post-induction flux (mmol/gDW/h)", "max volumetric rate (mM/h)", "peak utilization (v/UB)"];
    const rows = data.results.map((r) => [
      fmt(r.etot_uM, 3), fmt(r.etot_factor, 3), fmt(r.final_limonene_mM, 2), fmt(r.max_post_iptg_flux, 2),
      fmt(r.mean_post_induction_flux_mmol_gDW_h, 4), fmt(r.max_volumetric_rate_mM_h, 4),
      fmt(r.peak_capacity_utilization, 4),
    ]);
    html = buildTableHtml(headers, rows);
  } else if (mode === "scan2d") {
    if (data.method) {
      // Sampled results: each row is one sample point, listing every enzyme's Etot value plus the two result values directly.
      const headers = [...data.enzymes.map((e) => `${e}_Etot_uM`), "final_limonene_mM", "max_post_iptg_flux"];
      const rows = data.results.map((r) => [
        ...data.enzymes.map((e) => fmt(r.etot[e], 3)),
        fmt(r.final_limonene_mM, 3), fmt(r.max_post_iptg_flux, 4),
      ]);
      html = buildTableHtml(headers, rows);
    } else {
      const enzymeX = data.enzymes[0];
      const enzymeY = data.enzymes[1];
      const { xVals, yVals, zFinal } = extractScan2dSlice(data, state.scan2dSliceIndices || []);
      const headers = [`${enzymeY} \\ ${enzymeX}`, ...xVals.map((v) => fmt(v, 3))];
      const rows = yVals.map((y, i) => [
        fmt(y, 3), ...zFinal[i].map((v) => fmt(v, 3)),
      ]);
      html = `<p style="padding:8px 12px;margin:0;font-size:11px;color:var(--text-muted);">${t("rawDataScan2dNote")}</p>` +
        buildTableHtml(headers, rows);
    }
  } else if (mode === "optimize") {
    const headers = ["#", "DXS", "IDI", "GPPS", "LS", "limonene_flux (mmol/gDW/h)"];
    const rows = (data.history || []).map((h, i) => [
      i + 1, fmt(h.etot.DXS, 3), fmt(h.etot.IDI, 3), fmt(h.etot.GPPS, 3), fmt(h.etot.LS, 3),
      fmt(h.limonene_flux, 5),
    ]);
    html = buildTableHtml(headers, rows);
  }
  container.innerHTML = html;
}

/* ---------------- Etot scan ---------------- */