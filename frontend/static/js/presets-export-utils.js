function getPresetsStore() {
  try {
    return JSON.parse(localStorage.getItem(PRESET_STORAGE_KEY) || "{}");
  } catch (err) {
    return {};
  }
}

function savePresetsStore(store) {
  localStorage.setItem(PRESET_STORAGE_KEY, JSON.stringify(store));
}

function refreshPresetSelect() {
  const store = getPresetsStore();
  const names = Object.keys(store).sort();
  els.presetSelect.innerHTML = `<option value="">${t("presetNone")}</option>` +
    names.map((n) => `<option value="${escapeHtml(n)}">${escapeHtml(n)}</option>`).join("");
}

function escapeHtml(s) {
  return String(s).replace(/[&<>"']/g, (c) => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]
  ));
}

function getCurrentParamSnapshot() {
  return {
    mode: state.mode,
    etot: readEtot(),
    iptgTime: document.getElementById("iptgTime").value,
    minGrowthFrac: document.getElementById("minGrowthFrac").value,
    totalTime: document.getElementById("totalTime").value,
    induced: document.getElementById("induced").value,
    scanEnzyme: document.getElementById("scanEnzyme").value,
    scanValues: document.getElementById("scanValues").value,
    scanMultiCount: document.getElementById("scanMultiCount").value,
    scanStrategy: document.getElementById("scanStrategy").value,
    scanNSamples: document.getElementById("scanNSamples").value,
    scanEnzymeX: document.getElementById("scanEnzymeX").value,
    scanValuesX: document.getElementById("scanValuesX").value,
    scanEnzymeY: document.getElementById("scanEnzymeY").value,
    scanValuesY: document.getElementById("scanValuesY").value,
    scanEnzymeZ: document.getElementById("scanEnzymeZ").value,
    scanValuesZ: document.getElementById("scanValuesZ").value,
    scanEnzymeW: document.getElementById("scanEnzymeW").value,
    scanValuesW: document.getElementById("scanValuesW").value,
    glycerolUptake: document.getElementById("glycerolUptake").value,
    oxygenUptake: document.getElementById("oxygenUptake").value,
    aaUptake: document.getElementById("aaUptake").value,
    atpmLb: document.getElementById("atpmLb").value,
    dxpsUb: document.getElementById("dxpsUb").value,
    metatUb: document.getElementById("metatUb").value,
    fppsUb: document.getElementById("fppsUb").value,
  };
}

function applyParamSnapshot(snap) {
  if (!snap) return;
  if (snap.etot) {
    for (const enzyme of ["DXS", "IDI", "GPPS", "LS"]) {
      if (snap.etot[enzyme] !== undefined) document.getElementById(`etot-${enzyme}`).value = snap.etot[enzyme];
    }
  }
  const setIf = (id, v) => { if (v !== undefined && v !== null && v !== "") document.getElementById(id).value = v; };
  setIf("iptgTime", snap.iptgTime);
  setIf("minGrowthFrac", snap.minGrowthFrac);
  setIf("totalTime", snap.totalTime);
  setIf("induced", snap.induced);
  setIf("scanEnzyme", snap.scanEnzyme);
  setIf("scanValues", snap.scanValues);
  setIf("scanMultiCount", snap.scanMultiCount);
  setIf("scanStrategy", snap.scanStrategy);
  setIf("scanNSamples", snap.scanNSamples);
  setIf("scanEnzymeX", snap.scanEnzymeX);
  setIf("scanValuesX", snap.scanValuesX);
  setIf("scanEnzymeY", snap.scanEnzymeY);
  setIf("scanValuesY", snap.scanValuesY);
  setIf("scanEnzymeZ", snap.scanEnzymeZ);
  setIf("scanValuesZ", snap.scanValuesZ);
  setIf("scanEnzymeW", snap.scanEnzymeW);
  setIf("scanValuesW", snap.scanValuesW);
  updateScanMultiFieldsVisibility();
  updateScanStrategyUI();
  updateScanEnzymeDisabledStates();
  setIf("glycerolUptake", snap.glycerolUptake);
  setIf("oxygenUptake", snap.oxygenUptake);
  setIf("aaUptake", snap.aaUptake);
  setIf("atpmLb", snap.atpmLb);
  setIf("dxpsUb", snap.dxpsUb);
  setIf("metatUb", snap.metatUb);
  setIf("fppsUb", snap.fppsUb);
  if (snap.mode) {
    const btn = els.modeTabs.querySelector(`.mode-tab[data-mode="${snap.mode}"]`);
    if (btn) btn.click();
  }
}

function bindPresetControls() {
  document.getElementById("presetSaveBtn").addEventListener("click", () => {
    const name = els.presetName.value.trim();
    if (!name) { flashHint(t("presetNameRequired"), true); return; }
    const store = getPresetsStore();
    store[name] = getCurrentParamSnapshot();
    savePresetsStore(store);
    refreshPresetSelect();
    els.presetSelect.value = name;
    flashHint(t("presetSavedMsg", name), false);
  });

  document.getElementById("presetLoadBtn").addEventListener("click", () => {
    const name = els.presetSelect.value;
    if (!name) { flashHint(t("presetSelectRequired"), true); return; }
    const store = getPresetsStore();
    if (store[name]) {
      applyParamSnapshot(store[name]);
      flashHint(t("presetLoadedMsg", name), false);
    }
  });

  document.getElementById("presetDeleteBtn").addEventListener("click", () => {
    const name = els.presetSelect.value;
    if (!name) { flashHint(t("presetSelectRequired"), true); return; }
    const store = getPresetsStore();
    delete store[name];
    savePresetsStore(store);
    refreshPresetSelect();
    flashHint(t("presetDeletedMsg", name), false);
  });

  document.getElementById("shareLinkBtn").addEventListener("click", () => {
    const url = buildShareUrl();
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(url)
        .then(() => flashHint(t("shareLinkCopied"), false))
        .catch(() => window.prompt(t("shareLinkCopied"), url));
    } else {
      window.prompt(t("shareLinkCopied"), url);
    }
  });
}

function bindCompareControls() {
  els.compareToggle.addEventListener("change", () => {
    state.compare = els.compareToggle.checked;
  });
  document.getElementById("compareClearBtn").addEventListener("click", () => {
    state.dfbaHistory = [];
    els.compareClearRow.hidden = true;
    flashHint(t("compareClearedMsg"), false);
  });
}

function flashHint(text, isError) {
  els.statusHint.className = isError ? "hint hint--error" : "hint";
  els.statusHint.textContent = text;
}

/* ================= Share link (URL parameters) ================= */
function buildShareUrl() {
  const snap = getCurrentParamSnapshot();
  const params = new URLSearchParams();
  params.set("mode", snap.mode);
  params.set("lang", state.lang);
  params.set("dxs", snap.etot.DXS);
  params.set("idi", snap.etot.IDI);
  params.set("gpps", snap.etot.GPPS);
  params.set("ls", snap.etot.LS);
  params.set("iptg", snap.iptgTime);
  params.set("mgf", snap.minGrowthFrac);
  params.set("tt", snap.totalTime);
  params.set("induced", snap.induced);
  params.set("se", snap.scanEnzyme);
  params.set("sv", snap.scanValues);
  params.set("smc", snap.scanMultiCount);
  params.set("sst", snap.scanStrategy);
  params.set("sns", snap.scanNSamples);
  params.set("sex", snap.scanEnzymeX);
  params.set("svx", snap.scanValuesX);
  params.set("sey", snap.scanEnzymeY);
  params.set("svy", snap.scanValuesY);
  params.set("sez", snap.scanEnzymeZ);
  params.set("svz", snap.scanValuesZ);
  params.set("sew", snap.scanEnzymeW);
  params.set("svw", snap.scanValuesW);
  params.set("glyc", snap.glycerolUptake);
  params.set("o2", snap.oxygenUptake);
  params.set("aa", snap.aaUptake);
  params.set("atpm", snap.atpmLb);
  params.set("dxps", snap.dxpsUb);
  params.set("metat", snap.metatUb);
  params.set("fpps", snap.fppsUb);
  const url = new URL(window.location.href);
  url.search = params.toString();
  return url.toString();
}

function applyUrlParams() {
  const params = new URLSearchParams(window.location.search);
  if (![...params.keys()].length) return;

  if (params.has("lang") && params.get("lang") !== state.lang) {
    state.lang = params.get("lang") === "en" ? "en" : "zh";
    applyStaticI18n();
  }
  applyParamSnapshot({
    mode: params.get("mode") || undefined,
    etot: {
      DXS: params.get("dxs") || undefined,
      IDI: params.get("idi") || undefined,
      GPPS: params.get("gpps") || undefined,
      LS: params.get("ls") || undefined,
    },
    iptgTime: params.get("iptg") || undefined,
    minGrowthFrac: params.get("mgf") || undefined,
    totalTime: params.get("tt") || undefined,
    induced: params.get("induced") || undefined,
    scanEnzyme: params.get("se") || undefined,
    scanValues: params.get("sv") || undefined,
    scanMultiCount: params.get("smc") || undefined,
    scanStrategy: params.get("sst") || undefined,
    scanNSamples: params.get("sns") || undefined,
    scanEnzymeX: params.get("sex") || undefined,
    scanValuesX: params.get("svx") || undefined,
    scanEnzymeY: params.get("sey") || undefined,
    scanValuesY: params.get("svy") || undefined,
    scanEnzymeZ: params.get("sez") || undefined,
    scanValuesZ: params.get("svz") || undefined,
    scanEnzymeW: params.get("sew") || undefined,
    scanValuesW: params.get("svw") || undefined,
    glycerolUptake: params.get("glyc") || undefined,
    oxygenUptake: params.get("o2") || undefined,
    aaUptake: params.get("aa") || undefined,
    atpmLb: params.get("atpm") || undefined,
    dxpsUb: params.get("dxps") || undefined,
    metatUb: params.get("metat") || undefined,
    fppsUb: params.get("fpps") || undefined,
  });
}

/* ================= Result export ================= */
function bindExportControls() {
  document.querySelectorAll("[data-export-csv]").forEach((btn) => {
    btn.addEventListener("click", () => exportCsv(btn.dataset.exportCsv));
  });
  document.querySelectorAll("[data-export-png]").forEach((btn) => {
    btn.addEventListener("click", () => exportPng(btn.dataset.exportPng));
  });
  document.getElementById("fbaReportCopyBtn").addEventListener("click", () => {
    const text = document.getElementById("fbaReportBlock").textContent;
    if (!text) { flashHint(t("exportNoData"), true); return; }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text)
        .then(() => flashHint(t("reportCopiedMsg"), false))
        .catch(() => window.prompt(t("reportCopiedMsg"), text));
    } else {
      window.prompt(t("reportCopiedMsg"), text);
    }
  });
}

function csvEscape(v) {
  const s = v === null || v === undefined ? "" : String(v);
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

function rowsToCsv(headers, rows) {
  return [headers, ...rows].map((r) => r.map(csvEscape).join(",")).join("\n");
}

function downloadCsv(filename, csvContent) {
  const blob = new Blob(["\uFEFF" + csvContent], { type: "text/csv;charset=utf-8;" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  URL.revokeObjectURL(url);
}

function exportCsv(mode) {
  const data = (mode === "scan" && state.scanDisplayMode === "all")
    ? state.lastResults.scanAll : state.lastResults[mode];
  if (!data) { flashHint(t("exportNoData"), true); return; }

  let csv, filename;
  if (mode === "fba") {
    const rows = Object.entries(data.enzyme_fluxes || {}).map(([k, v]) => [k, v]);
    rows.push(["mu_max_h-1", data.mu_max]);
    rows.push(["growth_rate_h-1", data.growth_rate]);
    rows.push(["limonene_flux_mmol_gDW_h", data.limonene_flux]);
    const cof = data.cofactors || {};
    rows.push(["atp_production_mmol_gDW_h", cof.atp_production_mmol_gDW_h]);
    rows.push(["nadph_production_mmol_gDW_h", cof.nadph_production_mmol_gDW_h]);
    const mp = data.medium_pathway_report || {};
    rows.push(["glycerol_bounds", mp.glycerol_bounds ? mp.glycerol_bounds.join("/") : ""]);
    rows.push(["oxygen_bounds", mp.oxygen_bounds ? mp.oxygen_bounds.join("/") : ""]);
    rows.push(["atpm_lb", mp.atpm_lb]);
    rows.push(["dxps_ub", mp.dxps_ub]);
    rows.push(["metat_ub", mp.metat_ub]);
    rows.push(["fpps_ub", mp.fpps_ub]);
    csv = rowsToCsv(["metric", "value"], rows);
    filename = "fba_result.csv";
  } else if (mode === "dfba") {
    const headers = ["time_h", "phase", "solve_mode", "limonene_mM", "biomass_gDW_L", "glycerol_mM", "growth_rate_h", "limonene_flux"];
    const rows = data.time_h.map((tv, i) => [
      tv, data.phase?.[i] ?? "", data.mode?.[i] ?? "",
      data.limonene_mM[i], data.biomass_gDW_L[i], data.glycerol_mM[i],
      data.growth_rate_h[i], data.limonene_flux[i],
    ]);
    csv = rowsToCsv(headers, rows);
    filename = "dfba_result.csv";
  } else if (mode === "scan") {
    if (data.enzymes) {
      // Four-enzyme overlay-comparison mode: flattened into one full table, with an extra column per row naming which enzyme it is.
      const headers = ["enzyme", "etot_uM", "etot_factor", "final_limonene_mM", "max_post_iptg_flux",
        "mean_post_induction_flux_mmol_gDW_h", "max_volumetric_rate_mM_h", "peak_capacity_utilization"];
      const rows = [];
      for (const enzyme of data.enzymes) {
        for (const r of data.allResults[enzyme]) {
          rows.push([enzyme, r.etot_uM, r.etot_factor, r.final_limonene_mM, r.max_post_iptg_flux,
            r.mean_post_induction_flux_mmol_gDW_h, r.max_volumetric_rate_mM_h, r.peak_capacity_utilization]);
        }
      }
      csv = rowsToCsv(headers, rows);
      filename = "scan_all_enzymes_result.csv";
    } else {
      const headers = ["etot_uM", "etot_factor", "final_limonene_mM", "max_post_iptg_flux",
        "mean_post_induction_flux_mmol_gDW_h", "max_volumetric_rate_mM_h", "peak_capacity_utilization"];
      const rows = data.results.map((r) => [
        r.etot_uM, r.etot_factor, r.final_limonene_mM, r.max_post_iptg_flux,
        r.mean_post_induction_flux_mmol_gDW_h, r.max_volumetric_rate_mM_h, r.peak_capacity_utilization,
      ]);
      csv = rowsToCsv(headers, rows);
      filename = `scan_${data.enzyme}_result.csv`;
    }
  } else if (mode === "scan2d") {
    const headers = [...data.enzymes.map((e) => `${e}_Etot_uM`), "final_limonene_mM", "max_post_iptg_flux"];
    let rows;
    if (data.method) {
      // Sampled results are already a flat list; convert straight to CSV rows.
      rows = data.results.map((r) => [
        ...data.enzymes.map((e) => r.etot[e]), r.final_limonene_mM, r.max_post_iptg_flux,
      ]);
    } else {
      // An exhaustive grid is a nested structure; the CSV export includes the full N-dimensional data (unlike the table, which only shows the current slice).
      rows = [];
      const walk = (idxPath, nodeFinal, nodeFlux) => {
        if (idxPath.length === data.enzymes.length) {
          const etotVals = idxPath.map((idx, dim) => data.values_list[dim][idx]);
          rows.push([...etotVals, nodeFinal, nodeFlux]);
          return;
        }
        for (let i = 0; i < nodeFinal.length; i++) {
          walk([...idxPath, i], nodeFinal[i], nodeFlux[i]);
        }
      };
      walk([], data.final_limonene_mM, data.max_post_iptg_flux);
    }
    csv = rowsToCsv(headers, rows);
    filename = `scan_multi_${data.enzymes.join("_")}_result.csv`;
  } else if (mode === "optimize") {
    const headers = ["eval", "DXS", "IDI", "GPPS", "LS", "limonene_flux_mmol_gDW_h"];
    const rows = (data.history || []).map((h, i) => [
      i + 1, h.etot.DXS, h.etot.IDI, h.etot.GPPS, h.etot.LS, h.limonene_flux,
    ]);
    csv = rowsToCsv(headers, rows);
    filename = "optimize_etot_result.csv";
  } else {
    return;
  }
  downloadCsv(filename, csv);
}

function exportPng(blockId) {
  const block = document.getElementById(blockId);
  const charts = Array.from(block.querySelectorAll(".chart")).filter((el) => el.data);
  if (!charts.length) { flashHint(t("exportNoData"), true); return; }
  charts.forEach((chartDiv, i) => {
    setTimeout(() => {
      Plotly.downloadImage(chartDiv, { format: "png", filename: chartDiv.id, width: 1000, height: 600 });
    }, i * 400);
  });
}

function fmt(v, digits) {
  if (v === null || v === undefined || Number.isNaN(v)) return "—";
  return Number(v).toFixed(digits);
}

// Bootstraps the whole page. Called here, at the end of the last-loaded
// script, so that every function init() depends on (across all the other
// js/*.js files) is already defined by the time it runs.
init();