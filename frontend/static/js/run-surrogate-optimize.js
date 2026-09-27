async function runSurrogate() {
  // Uses the shared "total enzyme concentration Etot" field (readEtot())
  // directly, rather than opening a separate set of input boxes -- so
  // every tab sees the same set of numbers, and changing it in one place
  // can't drift out of sync with another.
  const etot = readEtot();
  const [data, curve] = await Promise.all([
    postJson("/api/etot-surrogate", { etot }),
    postJson("/api/etot-surrogate-curve", { etot: { DXS: etot.DXS, IDI: etot.IDI, GPPS: etot.GPPS }, n_points: 150 }),
  ]);
  state.lastResults.surrogate = data;
  state.lastResults.surrogateCurve = curve;
  state.lastResults.surrogateLS = etot.LS;
  renderSurrogateResults(data, curve, etot.LS);
}

/* The surrogate model is fundamentally a piecewise-cubic equation that
   takes LS Etot as input and outputs a limonene prediction -- DXS/IDI/GPPS
   don't affect the actual output value (they only affect the domain
   check) -- so the most meaningful chart is the full LS response curve
   (log-log), not a single point. The currently queried LS value is
   highlighted with a marker, and "termination classification uncertain"
   ranges like 1000-1250 uM are marked with a translucent band, echoing the
   hatched-band idea from the reference chart. */
function renderSurrogateCurveChart(curve, currentLS) {
  const points = curve.points;
  const shapes = (curve.termination_brackets_uM || []).map(([lo, hi]) => ({
    type: "rect", xref: "x", yref: "paper",
    x0: lo, x1: hi, y0: 0, y1: 1,
    fillcolor: accentColors()[3], opacity: 0.15, line: { width: 0 },
  }));

  const currentPoint = points.reduce((closest, p) =>
    Math.abs(p.LS_uM - currentLS) < Math.abs(closest.LS_uM - currentLS) ? p : closest, points[0]);

  // The surrogate model's equation fundamentally only consumes LS
  // (DXS/IDI/GPPS only affect the domain check, not the actual output
  // value -- as equations.md itself states), so no matter what these three
  // enzymes are set to, the prediction always equals the current query
  // point's value. Three horizontal reference lines are overlaid to
  // directly visualize "these three enzymes have zero influence in this
  // surrogate model," echoing the V9 reference chart layout you've seen
  // (the LSxDXS / LSxIDI / LSxGPPS projections, where DXS/IDI/GPPS are all
  // the same flat line).
  const xMin = Math.min(...points.map((p) => p.LS_uM));
  const xMax = Math.max(...points.map((p) => p.LS_uM));
  const flatTraces = ["DXS", "IDI", "GPPS"].map((enzyme, i) => ({
    x: [xMin, xMax],
    y: [currentPoint.limonene_mM, currentPoint.limonene_mM],
    type: "scatter", mode: "lines",
    line: { color: accentColors()[(i + 1) % accentColors().length], width: 1.5, dash: "dot" },
    name: enzyme,
    hovertemplate: `${enzyme}: %{y:.4g} mM (${t("surrogateFlatLineNote")})<extra></extra>`,
  }));

  Plotly.newPlot("surrogateCurveChart", [
    {
      x: points.map((p) => p.LS_uM),
      y: points.map((p) => p.limonene_mM),
      type: "scatter", mode: "lines",
      line: { color: accentColors()[0], width: 2 },
      name: "LS",
      hovertemplate: "LS=%{x:.4g} µM<br>%{y:.4g} mM<extra></extra>",
    },
    ...flatTraces,
    {
      x: [currentPoint.LS_uM], y: [currentPoint.limonene_mM],
      type: "scatter", mode: "markers",
      marker: { color: accentColors()[3], size: 12, symbol: "diamond" },
      name: t("surrogateCurrentQueryLegend"),
      hovertemplate: "LS=%{x:.4g} µM<br>%{y:.4g} mM<extra></extra>",
    },
  ], {
    ...plotlyLayoutBase(),
    title: { text: t("surrogateCurveTitle"), font: { size: 13, color: chartTitleColor() } },
    // The 5 legend entries (LS/DXS/IDI/GPPS/current query point) share the
    // same problem as the T7 flux chart and the Etot scan trace chart: a
    // horizontal legend with long labels competes with the X-axis title
    // for the same vertical space and overlaps it.
    margin: { ...plotlyLayoutBase().margin, b: 90 },
    legend: { ...plotlyLayoutBase().legend, y: -0.35 },
    xaxis: { ...plotlyLayoutBase().xaxis, title: "Etot (µM)", type: "log" },
    yaxis: { ...plotlyLayoutBase().yaxis, title: "mM", type: "log" },
    shapes,
    showlegend: true,
  }, { displayModeBar: false, responsive: true });
}

function renderSurrogateResults(data, curve, currentLS) {
  document.getElementById("mSurrogateLimonene").textContent =
    `${fmt(data.limonene_mM, 4)} mM`;

  const terminationKey = {
    horizon: "surrogateTerminationHorizon",
    model_infeasible_before_25h: "surrogateTerminationInfeasible",
    boundary_uncertain: "surrogateTerminationUncertain",
  }[data.termination] || data.termination;
  document.getElementById("mSurrogateTermination").textContent = t(terminationKey);

  document.getElementById("mSurrogateEndpoint").textContent =
    data.endpoint_time_h != null ? `${fmt(data.endpoint_time_h, 2)} h` : "—";

  document.getElementById("surrogateSummaryLine").textContent =
    `${t("surrogateSummaryPrefix")}${fmt(data.induction_h, 2)} h`;

  if (curve) renderSurrogateCurveChart(curve, currentLS);

  // provenance_note is a backend-hardcoded sentence (identical to the blue
  // hint box above it), doesn't get translated on a language switch, and
  // putting it in the raw data table would just produce the same sentence
  // appearing once in Chinese and once translated, an inconsistent
  // duplicate -- it's already shown via the translated hint box, so the raw
  // table doesn't need to repeat it; filter it out directly.
  const rows = Object.entries(data)
    .filter(([k]) => k !== "provenance_note")
    .map(([k, v]) => [k, JSON.stringify(v)]);
  document.getElementById("surrogateDataTable").innerHTML =
    buildTableHtml([t("rawDataMetricCol"), t("rawDataValueCol")], rows);
}

/* ---------------- Auto-optimize Etot (Nelder-Mead) ---------------- */
async function runOptimize() {
  const maxIter = Number(document.getElementById("optimizeMaxIter").value) || 200;
  const params = readCommonParams();
  const data = await postJson("/api/optimize-etot", { params, max_iter: maxIter });
  state.lastResults.optimize = data;
  renderOptimizeResults(data);
}

function renderOptimizeResults(data) {
  document.getElementById("mOptimizeFlux").textContent =
    `${fmt(data.best_limonene_flux, 5)} mmol/gDW/h`;
  document.getElementById("mOptimizeEvals").textContent = String(data.n_evaluations);
  document.getElementById("mOptimizeConverged").textContent =
    data.converged ? t("optimizeConvergedYes") : t("optimizeConvergedNo");

  // Best Etot combination: one bar per enzyme, colored in sequence.
  const enzymes = ["DXS", "IDI", "GPPS", "LS"];
  const values = enzymes.map((e) => data.best_etot[e]);
  Plotly.newPlot("optimizeBestEtotChart", [{
    type: "bar",
    x: enzymes,
    y: values,
    text: values.map((v) => fmt(v, 3)),
    textposition: "outside",
    cliponaxis: false,
    marker: { color: categoricalColors(enzymes.length) },
  }], {
    ...plotlyLayoutBase(),
    title: { text: t("optimizeBestEtotChartTitle"), font: { size: 13, color: chartTitleColor() } },
    margin: { l: 55, r: 20, t: 60, b: 45 },
    xaxis: { ...plotlyLayoutBase().xaxis, type: "category" },
    yaxis: { ...plotlyLayoutBase().yaxis, title: "Etot (µM)" },
  }, { displayModeBar: false, responsive: true });

  // Convergence trace: the limonene flux at each evaluation, with a
  // "running best so far" cumulative-max line overlaid, making it easy to
  // see whether the algorithm has converged (the curve flattens) or is
  // still improving.
  const history = data.history || [];
  const evalIdx = history.map((_, i) => i + 1);
  const fluxSeries = history.map((h) => h.limonene_flux);
  let runningBest = -Infinity;
  const bestSeries = fluxSeries.map((v) => {
    runningBest = Math.max(runningBest, v);
    return runningBest;
  });
  Plotly.newPlot("optimizeConvergeChart", [
    { x: evalIdx, y: fluxSeries, name: t("limFluxLegend"), mode: "markers",
      marker: { color: accentColors()[2], size: 5, opacity: 0.6 } },
    { x: evalIdx, y: bestSeries, name: t("optimizeBestFluxTrace"),
      line: { color: accentColors()[0], width: 2.5 } },
  ], {
    ...plotlyLayoutBase(),
    title: { text: t("optimizeConvergeChartTitle"), font: { size: 13, color: chartTitleColor() } },
    xaxis: { ...plotlyLayoutBase().xaxis, title: "evaluation #" },
    yaxis: { ...plotlyLayoutBase().yaxis, title: "mmol/gDW/h" },
  }, { displayModeBar: false, responsive: true });

  renderDataTable("optimize", data);
}

document.getElementById("optimizeApplyBtn").addEventListener("click", () => {
  const data = state.lastResults.optimize;
  if (!data) { flashHint(t("exportNoData"), true); return; }
  for (const enzyme of ["DXS", "IDI", "GPPS", "LS"]) {
    document.getElementById(`etot-${enzyme}`).value = fmt(data.best_etot[enzyme], 4);
  }
  flashHint(t("optimizeApplied"), false);
});