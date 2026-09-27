async function runScan() {
  const enzyme = document.getElementById("scanEnzyme").value;
  const scanValues = document.getElementById("scanValues").value
    .split(",").map((s) => Number(s.trim())).filter((n) => !Number.isNaN(n));
  const params = readCommonParams();
  const data = await postJson("/api/etot-scan", { enzyme, scan_values: scanValues, params });
  state.lastResults.scan = { enzyme, results: data.results, params };
  state.scanDisplayMode = "single";
  renderScanCharts(state.lastResults.scan);
}

/* Plain-language hint: which enzymes this scan is holding fixed (and at
   what value), and which one is being scanned -- to avoid the earlier
   confusion where, say, DXS/LS were already saturated and the user
   couldn't tell why the numbers weren't changing. */
function renderScanFixedSummary(elId, scanningEnzymes, params) {
  const fixedParts = ["DXS", "IDI", "GPPS", "LS"]
    .filter((e) => !scanningEnzymes.includes(e))
    .map((e) => `${e}=${fmt(params.Etot[e], 3)} µM`);
  const scanningLabel = scanningEnzymes.join(" × ");
  const text = `${t("scanScanningSummaryPrefix")}${scanningLabel}　${t("scanFixedSummaryPrefix")}${fixedParts.join("、")}`;
  document.getElementById(elId).textContent = text;
}

/* If a log-axis chart just lets Plotly auto-scale to "tightly fit the
   data's own min/max," then when the real spread is only a few decimal
   places (e.g. DXS/GPPS/IDI were never the bottleneck to begin with, so
   changing them barely matters), the chart gets stretched into a
   steep-looking line that misleadingly suggests a big change. This instead
   forces the log axis to span at least MIN_LOG_DECADES orders of
   magnitude: a line with genuinely little change now looks flat, and only
   a line that truly spans many orders of magnitude shows a visible slope
   -- so "how steep the line looks" faithfully reflects "how big the actual
   difference is," and different enzymes/metrics can be fairly compared
   against each other. */
const MIN_LOG_DECADES = 4;
function computeLogAxisRange(values) {
  const positives = values.filter((v) => Number.isFinite(v) && v > 0);
  if (!positives.length) return undefined;
  const logMin = Math.log10(Math.min(...positives));
  const logMax = Math.log10(Math.max(...positives));
  const span = logMax - logMin;
  if (span >= MIN_LOG_DECADES) {
    // The data itself already spans enough orders of magnitude; just leave a small margin instead of forcing an expansion.
    const pad = Math.max(0.1, span * 0.05);
    return [logMin - pad, logMax + pad];
  }
  const center = (logMin + logMax) / 2;
  const half = MIN_LOG_DECADES / 2;
  return [center - half, center + half];
}

function renderScanCharts(scanData) {
  const enzyme = scanData.enzyme;
  renderScanFixedSummary("scanSummaryLine", [enzyme], scanData.params || readCommonParams());

  const trajTraces = scanData.results.map((r, i) => ({
    x: r.trace.time_h,
    y: r.trace.limonene_mM,
    name: `Etot = ${r.etot_uM} µM`,
    line: { color: accentColors()[i % accentColors().length], width: 2 },
  }));
  Plotly.newPlot("scanTrajChart", trajTraces, {
    ...plotlyLayoutBase(),
    title: { text: t("scanTrajTitle", enzyme), font: { size: 13, color: chartTitleColor() } },
    // Usually a 5-line (5 Etot values) horizontal legend, the same
    // label-length problem as the T7 flux chart: competing with the X-axis
    // title for the same vertical space overlaps them, so the bottom
    // margin is widened and the legend pushed further down.
    margin: { ...plotlyLayoutBase().margin, b: 90 },
    legend: { ...plotlyLayoutBase().legend, y: -0.35 },
    xaxis: { ...plotlyLayoutBase().xaxis, title: t("timeAxis") },
    yaxis: { ...plotlyLayoutBase().yaxis, title: t("extracellularLimAxis") },
  }, { displayModeBar: false, responsive: true });

  // Mirrors the V9 OAT capacity-scan four-panel diagnostic layout: the X
  // axis is uniformly "Etot fold-change relative to baseline," all
  // log-log; the four metrics each get their own chart instead of being
  // crammed into one chart with dual Y axes -- so you can see, all at
  // once, four different angles: "cumulative yield," "mean flux,"
  // "volumetric production rate," and "how saturated is this enzyme itself
  // (v/UB)." The last one in particular is a very direct diagnostic:
  // utilization sitting near 1 long-term means this enzyme is near full
  // load and is the bottleneck; staying low long-term means this enzyme
  // still has headroom and something else (not this enzyme) is holding
  // back yield.
  const etotFactors = scanData.results.map((r) => r.etot_factor ?? r.etot_uM);
  const xAxisLog = { ...plotlyLayoutBase().xaxis, title: t("etotFactorAxis"), type: "log" };

  const logChart = (containerId, yKey, titleKey, yTitle, color) => {
    const yVals = scanData.results.map((r) => r[yKey]);
    Plotly.newPlot(containerId, [{
      x: etotFactors,
      y: yVals,
      type: "scatter",
      mode: "lines+markers",
      line: { color },
      marker: { size: 7 },
    }], {
      ...plotlyLayoutBase(),
      title: { text: t(titleKey), font: { size: 13, color: chartTitleColor() } },
      xaxis: xAxisLog,
      yaxis: { ...plotlyLayoutBase().yaxis, title: yTitle, type: "log", range: computeLogAxisRange(yVals) },
    }, { displayModeBar: false, responsive: true });
  };

  logChart("scanAccumChart", "final_limonene_mM", "scanAccumTitle", "mM", accentColors()[0]);
  logChart("scanMeanFluxChart", "mean_post_induction_flux_mmol_gDW_h", "scanMeanFluxTitle",
    "mmol/gDW/h", accentColors()[1]);
  logChart("scanVolRateChart", "max_volumetric_rate_mM_h", "scanVolRateTitle",
    "mM/h", accentColors()[2]);
  logChart("scanUtilizationChart", "peak_capacity_utilization", "scanUtilizationTitle",
    "v / UB", accentColors()[3]);

  renderDataTable("scan", scanData);
}

/* Mirrors the reference chart (the V9 Etot one-at-a-time capacity scan):
   that chart actually overlays the four results from "scanning DXS/IDI/
   GPPS/LS each independently, one at a time" onto the same set of 4
   panels, rather than looking at just one enzyme. This calls /api/etot-scan
   once per enzyme in sequence (each using its own default value sequence
   from the registry table LIM009_params.csv, not whatever is currently
   typed into the scan-settings field), then overlays the results once
   they're back. The trace chart (concentration vs. time) is not drawn in
   this mode -- adding up each enzyme's own value-sequence point count
   would produce 15-20 overlapping lines, which would be unreadable anyway
   -- and this also matches the reference chart, which has no trace chart
   either, only these 4 diagnostic panels. */
async function runScanAllEnzymes() {
  const enzymes = ["DXS", "IDI", "GPPS", "LS"];
  els.statusHint.textContent = t("scanAllEnzymesRunning");

  const params = readCommonParams();
  const allResults = {};
  for (const enzyme of enzymes) {
    const data = await postJson("/api/etot-scan", { enzyme, params });
    allResults[enzyme] = data.results;
  }

  state.lastResults.scanAll = { enzymes, allResults, params };
  state.scanDisplayMode = "all";
  renderScanAllCharts(state.lastResults.scanAll);
}

function renderScanAllCharts(scanAllData) {
  const { enzymes, allResults } = scanAllData;

  document.getElementById("scanSummaryLine").textContent =
    `${t("scanScanningSummaryPrefix")}${enzymes.join(", ")}`;

  // The trace chart is cleared in this mode, to avoid four enzymes x 5 lines each overlapping into an unreadable mess.
  Plotly.newPlot("scanTrajChart", [], {
    ...plotlyLayoutBase(),
    title: { text: t("scanAllEnzymesBtn"), font: { size: 13, color: chartTitleColor() } },
    xaxis: { ...plotlyLayoutBase().xaxis, title: t("timeAxis") },
    yaxis: { ...plotlyLayoutBase().yaxis, title: t("extracellularLimAxis") },
  }, { displayModeBar: false, responsive: true });

  const logChartMulti = (containerId, yKey, titleKey, yTitle) => {
    const traces = enzymes.map((enzyme, i) => {
      const results = allResults[enzyme];
      return {
        x: results.map((r) => r.etot_factor ?? r.etot_uM),
        y: results.map((r) => r[yKey]),
        name: enzyme,
        type: "scatter",
        mode: "lines+markers",
        line: { color: accentColors()[i % accentColors().length] },
        marker: { size: 7 },
      };
    });
    // The 4 lines (4 enzymes) share the same Y-axis range, so their slopes
    // can be fairly compared -- if each enzyme computed its own range
    // independently, an otherwise-flat line could get magnified into
    // looking like it changes.
    const allYVals = traces.flatMap((tr) => tr.y);
    Plotly.newPlot(containerId, traces, {
      ...plotlyLayoutBase(),
      title: { text: t(titleKey), font: { size: 13, color: chartTitleColor() } },
      xaxis: { ...plotlyLayoutBase().xaxis, title: t("etotFactorAxis"), type: "log" },
      yaxis: { ...plotlyLayoutBase().yaxis, title: yTitle, type: "log", range: computeLogAxisRange(allYVals) },
      showlegend: true,
    }, { displayModeBar: false, responsive: true });
  };

  logChartMulti("scanAccumChart", "final_limonene_mM", "scanAccumTitle", "mM");
  logChartMulti("scanMeanFluxChart", "mean_post_induction_flux_mmol_gDW_h", "scanMeanFluxTitle", "mmol/gDW/h");
  logChartMulti("scanVolRateChart", "max_volumetric_rate_mM_h", "scanVolRateTitle", "mM/h");
  logChartMulti("scanUtilizationChart", "peak_capacity_utilization", "scanUtilizationTitle", "v / UB");

  // The raw data table instead lists the flattened full results for all four enzymes, handled separately from the single-enzyme mode's table.
  const container = document.getElementById("scanDataTable");
  if (container) {
    const headers = ["enzyme", "etot_uM", "etot_factor", "final_limonene_mM",
      "mean_post_induction_flux_mmol_gDW_h", "max_volumetric_rate_mM_h", "peak_capacity_utilization"];
    const rows = [];
    for (const enzyme of enzymes) {
      for (const r of allResults[enzyme]) {
        rows.push([enzyme, fmt(r.etot_uM, 3), fmt(r.etot_factor, 3), fmt(r.final_limonene_mM, 4),
          fmt(r.mean_post_induction_flux_mmol_gDW_h, 6), fmt(r.max_volumetric_rate_mM_h, 6),
          fmt(r.peak_capacity_utilization, 4)]);
      }
    }
    container.innerHTML = buildTableHtml(headers, rows);
  }
}

document.getElementById("scanAllEnzymesBtn").addEventListener("click", async () => {
  const btn = document.getElementById("scanAllEnzymesBtn");
  btn.disabled = true;
  try {
    await runScanAllEnzymes();
    els.statusHint.textContent = t("done");
  } catch (err) {
    els.statusHint.className = "hint hint--error";
    els.statusHint.textContent = t("errorPrefix", err.message);
  } finally {
    btn.disabled = false;
  }
});

/* ---------------- Multi-enzyme Etot cross scan (2-4 enzymes) ---------------- */
async function runScan2d() {
  const count = Number(document.getElementById("scanMultiCount").value) || 2;
  const axes = ["X", "Y", "Z", "W"].slice(0, count);
  const enzymes = axes.map((axis) => document.getElementById(`scanEnzyme${axis}`).value);
  if (new Set(enzymes).size !== enzymes.length) throw new Error(t("scan2dSameEnzymeError"));

  const valuesList = axes.map((axis) => document.getElementById(`scanValues${axis}`).value
    .split(",").map((s) => Number(s.trim())).filter((n) => !Number.isNaN(n)));

  const params = readCommonParams();
  const strategy = document.getElementById("scanStrategy").value;

  if (strategy === "grid") {
    const data = await postJson("/api/etot-scan-multi", { enzymes, values_list: valuesList, params });
    data.params = params;
    state.lastResults.scan2d = data;
    state.scan2dSliceIndices = new Array(Math.max(0, enzymes.length - 2)).fill(0);
    renderScan2dCharts(data);
    return;
  }

  // LHS/Sobol sampling mode: the "value sequence" field instead only takes the min/max as a range, rather than an enumerated list.
  const bounds = valuesList.map((vals) => [Math.min(...vals), Math.max(...vals)]);
  const nSamples = Number(document.getElementById("scanNSamples").value) || 50;
  const data = await postJson("/api/etot-scan-sampled", {
    enzymes, bounds, method: strategy, n_samples: nSamples, params,
  });
  data.params = params;
  state.lastResults.scan2d = data;
  renderScan2dCharts(data);
}

/* Extracts the "enzymes[0] x enzymes[1]" 2D slice from the N-dimensional
   nested array, holding the remaining dimensions (enzymes[2], enzymes[3],
   if present) fixed at the index given by sliceIndices. */
function extractScan2dSlice(data, sliceIndices) {
  const xVals = data.values_list[0];
  const yVals = data.values_list[1];
  const zFinal = [];
  const zFlux = [];
  for (let iy = 0; iy < yVals.length; iy++) {
    const rowFinal = [];
    const rowFlux = [];
    for (let ix = 0; ix < xVals.length; ix++) {
      const idx = [ix, iy, ...sliceIndices];
      let nodeFinal = data.final_limonene_mM;
      let nodeFlux = data.max_post_iptg_flux;
      for (const i of idx) { nodeFinal = nodeFinal[i]; nodeFlux = nodeFlux[i]; }
      rowFinal.push(nodeFinal);
      rowFlux.push(nodeFlux);
    }
    zFinal.push(rowFinal);
    zFlux.push(rowFlux);
  }
  return { xVals, yVals, zFinal, zFlux };
}

/* With 3-4 enzymes, the extra dimensions (enzymes[2]/enzymes[3]) can't be
   drawn into a 2D heatmap, so a dropdown instead lets the user pick a fixed
   value as a "slice" -- picking one doesn't re-hit the API, it just
   re-slices the already-fetched full result and redraws the chart. */
function renderScan2dSliceControls(data) {
  const container = els.scan2dSliceControls;
  container.innerHTML = "";
  if (data.enzymes.length <= 2) return;

  data.enzymes.slice(2).forEach((enzyme, dimOffset) => {
    const dim = dimOffset + 2;
    const label = document.createElement("label");
    label.textContent = `${t("scan2dSliceLabel")}${enzyme}`;
    label.style.marginRight = "6px";
    label.style.color = "var(--text-muted)";
    label.style.fontSize = "12.5px";
    label.style.fontFamily = "var(--font-mono)";
    label.style.display = "flex";
    label.style.alignItems = "center";
    label.style.gap = "6px";

    const select = document.createElement("select");
    data.values_list[dim].forEach((v, i) => {
      const opt = document.createElement("option");
      opt.value = i;
      opt.textContent = `${v}`;
      if (i === state.scan2dSliceIndices[dim - 2]) opt.selected = true;
      select.appendChild(opt);
    });
    select.addEventListener("change", () => {
      state.scan2dSliceIndices[dim - 2] = Number(select.value);
      renderScan2dCharts(state.lastResults.scan2d);
    });

    label.appendChild(select);
    container.appendChild(label);
  });
}

/* Dispatcher: if data has data.method ("lhs"/"sobol") it's a sampled
   result (scattered points, not a regular grid); otherwise it's an
   exhaustive grid result -- the two have completely different data
   structures and are drawn differently. */
function renderScan2dCharts(data) {
  if (data.method) {
    renderScan2dSampledCharts(data);
  } else {
    renderScan2dGridCharts(data);
  }
}

function renderScan2dGridCharts(data) {
  const enzymeX = data.enzymes[0];
  const enzymeY = data.enzymes[1];
  renderScanFixedSummary("scan2dSummaryLine", data.enzymes, data.params || readCommonParams());

  const dim2Count = data.enzymes.length >= 3 ? data.values_list[2].length : 1;
  const dim3Count = data.enzymes.length >= 4 ? data.values_list[3].length : 1;
  const facetCount = dim2Count * dim3Count;

  if (data.enzymes.length <= 2 || facetCount > 16) {
    // With 2 enzymes there's only ever one chart anyway; with more than 16
    // slices, small multiples would be too crowded to read -- both cases
    // fall back to the older "single slice + dropdown" mode.
    els.scan2dSliceControls.hidden = data.enzymes.length <= 2;
    renderScan2dSliceControls(data);
    renderSingleSliceHeatmaps(data, state.scan2dSliceIndices || []);
  } else {
    // With 3-4 enzymes and a manageable slice count: lay out every slice's small chart at once, no switching needed.
    els.scan2dSliceControls.hidden = true;
    els.scan2dSliceControls.innerHTML = "";
    renderSmallMultipleHeatmaps("scan2dFinalHeatmap", data, "zFinal",
      [[0, cssVar("--bg-panel")], [1, cssVar("--accent")]], "mM",
      t("scan2dFinalTitle", enzymeX, enzymeY));
    renderSmallMultipleHeatmaps("scan2dFluxHeatmap", data, "zFlux",
      [[0, cssVar("--bg-panel")], [1, cssVar("--leaf")]], "mmol/gDW/h",
      t("scan2dFluxTitle", enzymeX, enzymeY));
  }

  renderDataTable("scan2d", data);
}

/* Single-slice mode (used for 2 enzymes, or as a fallback when there are too many slices). */
function renderSingleSliceHeatmaps(data, sliceIndices) {
  const enzymeX = data.enzymes[0];
  const enzymeY = data.enzymes[1];
  const { xVals, yVals, zFinal, zFlux } = extractScan2dSlice(data, sliceIndices);
  const xLabels = xVals.map((v) => `${v}`);
  const yLabels = yVals.map((v) => `${v}`);

  Plotly.newPlot("scan2dFinalHeatmap", [{
    type: "heatmap",
    x: xLabels,
    y: yLabels,
    z: zFinal,
    colorscale: [[0, cssVar("--bg-panel")], [1, cssVar("--accent")]],
    colorbar: { title: "mM", titleside: "right" },
  }], {
    ...plotlyLayoutBase(),
    title: { text: t("scan2dFinalTitle", enzymeX, enzymeY), font: { size: 13, color: chartTitleColor() } },
    xaxis: { ...plotlyLayoutBase().xaxis, title: `${enzymeX} Etot (µM)` },
    yaxis: { ...plotlyLayoutBase().yaxis, title: `${enzymeY} Etot (µM)` },
  }, { displayModeBar: false, responsive: true });

  Plotly.newPlot("scan2dFluxHeatmap", [{
    type: "heatmap",
    x: xLabels,
    y: yLabels,
    z: zFlux,
    colorscale: [[0, cssVar("--bg-panel")], [1, cssVar("--leaf")]],
    colorbar: { title: "mmol/gDW/h", titleside: "right" },
  }], {
    ...plotlyLayoutBase(),
    title: { text: t("scan2dFluxTitle", enzymeX, enzymeY), font: { size: 13, color: chartTitleColor() } },
    xaxis: { ...plotlyLayoutBase().xaxis, title: `${enzymeX} Etot (µM)` },
    yaxis: { ...plotlyLayoutBase().yaxis, title: `${enzymeY} Etot (µM)` },
  }, { displayModeBar: false, responsive: true });
}

/* With 3-4 enzymes, lays out the heatmap for every combination of "the 3rd
   enzyme's value" (and "the 4th enzyme's value," if present) all at once
   using a Plotly subplot grid (small multiples), instead of a dropdown
   that switches one at a time. With 4 enzymes: columns = the 3rd enzyme's
   values, rows = the 4th enzyme's values. All the small charts share the
   same color-scale range (taken from the min/max across all slices), so
   color intensity can be compared directly between small charts, rather
   than being distorted by each one scaling independently. */
function renderSmallMultipleHeatmaps(containerId, data, zKey, colorscale, colorbarTitle, baseTitle) {
  const enzymes = data.enzymes;
  const dim2Enzyme = enzymes[2];
  const dim3Enzyme = enzymes.length === 4 ? enzymes[3] : null;
  const dim2Vals = data.values_list[2];
  const dim3Vals = dim3Enzyme ? data.values_list[3] : [null];

  const cols = dim2Vals.length;
  const rows = dim3Vals.length;

  const cells = [];
  let globalMin = Infinity;
  let globalMax = -Infinity;
  for (let r = 0; r < rows; r++) {
    for (let c = 0; c < cols; c++) {
      const sliceIndices = dim3Enzyme ? [c, r] : [c];
      const { xVals, yVals, zFinal, zFlux } = extractScan2dSlice(data, sliceIndices);
      const z = zKey === "zFinal" ? zFinal : zFlux;
      z.forEach((row) => row.forEach((v) => {
        if (v < globalMin) globalMin = v;
        if (v > globalMax) globalMax = v;
      }));
      cells.push({ r, c, xVals, yVals, z });
    }
  }

  const traces = [];
  const annotations = [];
  cells.forEach(({ r, c, xVals, yVals, z }) => {
    const idx = r * cols + c + 1; // Plotly subplot numbering starts at 1, row-major order
    const axisSuffix = idx === 1 ? "" : String(idx);
    traces.push({
      type: "heatmap",
      x: xVals.map((v) => `${v}`),
      y: yVals.map((v) => `${v}`),
      z,
      zmin: globalMin,
      zmax: globalMax,
      colorscale,
      showscale: false,
      xaxis: `x${axisSuffix}`,
      yaxis: `y${axisSuffix}`,
    });
    const dim2Label = `${dim2Enzyme}=${fmt(dim2Vals[c], 2)}`;
    const dim3Label = dim3Enzyme ? `, ${dim3Enzyme}=${fmt(dim3Vals[r], 2)}` : "";
    annotations.push({
      text: dim2Label + dim3Label,
      showarrow: false,
      xref: `x${axisSuffix} domain`,
      yref: `y${axisSuffix} domain`,
      x: 0.5, y: 1.16,
      font: { size: 10, color: chartTitleColor() },
    });
  });
  // The color scale is only shown on the last small chart and shared by all of them, instead of repeating a colorbar in every panel.
  if (traces.length) {
    traces[traces.length - 1].showscale = true;
    traces[traces.length - 1].colorbar = { title: colorbarTitle, titleside: "right" };
  }

  Plotly.newPlot(containerId, traces, {
    ...plotlyLayoutBase(),
    title: { text: baseTitle, font: { size: 13, color: chartTitleColor() } },
    grid: { rows, columns: cols, pattern: "independent", xgap: 0.10, ygap: 0.22 },
    annotations,
    height: Math.max(320, rows * 230 + 60),
    margin: { t: 70, b: 40, l: 50, r: 60 },
  }, { displayModeBar: false, responsive: true });
}

/* LHS/Sobol sampled results aren't a regular grid, so a heatmap doesn't
   work -- this uses a scatter plot instead: the X/Y axes are enzymes[0]/
   enzymes[1]'s actual sampled values, each point's color intensity
   represents the result value, and the 3rd/4th enzyme's value (if present)
   goes into the hover text rather than being used for slice filtering
   (the sample points are already sparse; filtering would leave even fewer
   visible points). */
function renderScan2dSampledCharts(data) {
  const enzymeX = data.enzymes[0];
  const enzymeY = data.enzymes[1];
  els.scan2dSliceControls.innerHTML = "";
  renderScanFixedSummary("scan2dSummaryLine", data.enzymes, data.params || readCommonParams());

  const xs = data.results.map((r) => r.etot[enzymeX]);
  const ys = data.results.map((r) => r.etot[enzymeY]);
  const finalVals = data.results.map((r) => r.final_limonene_mM);
  const fluxVals = data.results.map((r) => r.max_post_iptg_flux);
  const hoverText = data.results.map((r) => {
    const extra = data.enzymes.slice(2).map((e) => `${e}=${fmt(r.etot[e], 3)}`).join(", ");
    return `${enzymeX}=${fmt(r.etot[enzymeX], 3)}, ${enzymeY}=${fmt(r.etot[enzymeY], 3)}` +
      (extra ? `, ${extra}` : "");
  });
  const methodLabel = data.method === "lhs" ? t("scanStrategyLhs") : t("scanStrategySobol");

  Plotly.newPlot("scan2dFinalHeatmap", [{
    type: "scatter",
    mode: "markers",
    x: xs,
    y: ys,
    text: hoverText,
    hovertemplate: "%{text}<br>final limonene: %{marker.color:.3f} mM<extra></extra>",
    marker: {
      size: 11,
      color: finalVals,
      colorscale: [[0, cssVar("--bg-panel")], [1, cssVar("--accent")]],
      colorbar: { title: "mM", titleside: "right" },
      line: { width: 1, color: cssVar("--border") },
    },
  }], {
    ...plotlyLayoutBase(),
    title: { text: t("scanSampledResultTitle", methodLabel), font: { size: 13, color: chartTitleColor() } },
    xaxis: { ...plotlyLayoutBase().xaxis, title: `${enzymeX} Etot (µM)` },
    yaxis: { ...plotlyLayoutBase().yaxis, title: `${enzymeY} Etot (µM)` },
  }, { displayModeBar: false, responsive: true });

  Plotly.newPlot("scan2dFluxHeatmap", [{
    type: "scatter",
    mode: "markers",
    x: xs,
    y: ys,
    text: hoverText,
    hovertemplate: "%{text}<br>max flux: %{marker.color:.4f} mmol/gDW/h<extra></extra>",
    marker: {
      size: 11,
      color: fluxVals,
      colorscale: [[0, cssVar("--bg-panel")], [1, cssVar("--leaf")]],
      colorbar: { title: "mmol/gDW/h", titleside: "right" },
      line: { width: 1, color: cssVar("--border") },
    },
  }], {
    ...plotlyLayoutBase(),
    title: { text: t("maxFluxLegend") + ` (${methodLabel})`, font: { size: 13, color: chartTitleColor() } },
    xaxis: { ...plotlyLayoutBase().xaxis, title: `${enzymeX} Etot (µM)` },
    yaxis: { ...plotlyLayoutBase().yaxis, title: `${enzymeY} Etot (µM)` },
  }, { displayModeBar: false, responsive: true });

  renderDataTable("scan2d", data);
}

/* ---------------- Etot quick estimate (external surrogate model, not live COBRA) ---------------- */