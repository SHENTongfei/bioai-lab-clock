/* app.js — UI wiring: language switch, drag&drop, in-browser compute, ECharts render. */
(function () {
  "use strict";
  var C = window.Clock;
  var CHARTS = [];
  var state = { last: null };

  var $ = function (id) { return document.getElementById(id); };

  /* ---------- language ---------- */
  function markLang() {
    var code = getLang();
    document.querySelectorAll(".lang-btn").forEach(function (b) {
      b.classList.toggle("is-active", b.getAttribute("data-lang") === code);
    });
  }
  function setupLang() {
    document.querySelectorAll(".lang-btn").forEach(function (b) {
      b.addEventListener("click", function () {
        setLang(b.getAttribute("data-lang"));
        markLang();
        document.documentElement.lang = b.getAttribute("data-lang");
        if (state.last) render(state.last);
      });
    });
    markLang();
    window.__onLang = function () { /* re-render charts with new labels */ if (state.last) render(state.last); };
  }

  /* ---------- status ---------- */
  function showStatus(msg, isErr) {
    var s = $("status");
    s.hidden = !msg;
    s.classList.toggle("err", !!isErr);
    $("statusText").textContent = msg || "";
  }

  /* ---------- data loading ---------- */
  async function onDrop(files) {
    showStatus(t("computing"), false);
    try {
      var texts = [];
      var arr = Array.from(files);
      for (var i = 0; i < arr.length; i++) texts.push(await arr[i].text());
      var phenTxt = null, phenIdx = -1;
      for (var k = 0; k < texts.length; k++) {
        var ph = C.parsePhenotype(texts[k]);
        if (ph.ok && ph.ids.length > 3) { phenTxt = texts[k]; phenIdx = k; break; }
      }
      if (!phenTxt) { showStatus(t("errParse"), true); return; }
      var phenObj = C.parsePhenotype(phenTxt);
      var exTxt = null, exOk = null;
      for (var m = 0; m < texts.length; m++) {
        if (m === phenIdx) continue;
        var cand = C.parseExpression(texts[m], phenObj.ids);
        if (cand.ok) { exTxt = texts[m]; exOk = cand; break; }
      }
      if (!exTxt) { showStatus(t("errParse"), true); return; }
      run(exTxt, phenTxt, exOk, phenObj);
    } catch (e) {
      console.error(e);
      showStatus(t("errUnknown"), true);
    }
  }

  function setupDrop() {
    var dz = $("drop");
    var inp = $("fileInput");
    ["dragenter", "dragover"].forEach(function (ev) {
      dz.addEventListener(ev, function (e) { e.preventDefault(); dz.classList.add("drag"); });
    });
    ["dragleave", "drop"].forEach(function (ev) {
      dz.addEventListener(ev, function (e) { e.preventDefault(); dz.classList.remove("drag"); });
    });
    dz.addEventListener("drop", function (e) {
      if (e.dataTransfer && e.dataTransfer.files.length) onDrop(e.dataTransfer.files);
    });
    dz.addEventListener("click", function (e) {
      if (e.target.closest("button") || e.target.closest("a")) return;
      inp.click();
    });
    inp.addEventListener("change", function () {
      if (inp.files.length) onDrop(inp.files);
      inp.value = "";
    });
    // download demo (both files)
    $("btnDownload").addEventListener("click", function (e) {
      e.preventDefault();
      dl("data/demo_expression.csv");
      dl("data/demo_phenotype.csv");
    });
  }
  function dl(href) {
    var a = document.createElement("a");
    a.href = href; a.download = "";
    document.body.appendChild(a); a.click(); a.remove();
  }

  async function loadDemo() {
    showStatus(t("computing"), false);
    try {
      var res = await Promise.all([
        fetch("data/demo_expression.csv").then(function (r) { return r.text(); }),
        fetch("data/demo_phenotype.csv").then(function (r) { return r.text(); })
      ]);
      run(res[0], res[1]);
    } catch (e) {
      console.error(e);
      showStatus(t("errUnknown"), true);
    }
  }

  /* ---------- compute + render ---------- */
  function run(exprText, phenText) {
    var ph = C.parsePhenotype(phenText);
    if (!ph.ok) { showStatus(t("errParse"), true); return; }
    var ex = C.parseExpression(exprText, ph.ids);
    if (!ex.ok) { showStatus(t("errParse"), true); return; }
    var r = C.buildClock(ex, ph);
    if (!r.ok) { showStatus(r.error === "Fewer than 8 expression samples matched the phenotype." ? t("errMatch") : t("errUnknown"), true); return; }
    showStatus(null);
    state.last = r;
    render(r);
    $("results").hidden = false;
  }

  function bandResidL(r) {
    var L = r.perSample.filter(function (s) { return s.band === "L"; });
    if (!L.length) return NaN;
    var s = 0; L.forEach(function (x) { s += x.pred - x.age; });
    return s / L.length; // negative => reads younger
  }

  function render(r) {
    // KPIs
    $("kpiMae").textContent = r.maeEns.toFixed(2);
    $("kpiMaeSub").textContent = "± yr";
    $("kpiR2").textContent = r.r2.toFixed(3);
    $("kpiSpearman").textContent = r.rho.toFixed(3);
    $("kpiSamples").textContent = String(r.n);
    $("kpiSex").textContent = (r.sexValue >= 0 ? "+" : "") + r.sexValue.toFixed(2);
    var lr = bandResidL(r);
    $("kpiGate").textContent = isFinite(lr) ? "−" + Math.abs(lr).toFixed(1) + " yr" : "–";

    // gate card
    var hasGate = isFinite(r.gateAuc) && r.gateN > 0;
    $("gateAuc").textContent = hasGate ? r.gateAuc.toFixed(2) : "–";
    var note = t("gateNote", { v: isFinite(lr) ? Math.abs(lr).toFixed(1) : "–" });
    document.querySelectorAll('[data-i18n="gateNote"]').forEach(function (el) { el.textContent = note; });
    var gateLine = "n=" + r.gateN + " (90+)  ·  " + r.p + " genes  ·  λ=" + r.lambda;
    if (hasGate && r.gateAuc < 0.5) gateLine += "\n" + t("gateMirror", { m: (1 - r.gateAuc).toFixed(2) });
    $("gateN").textContent = gateLine;
    $("gateN").style.whiteSpace = "pre-line";

    // charts
    renderScatter(r);
    renderViews(r);
    renderBands(r);
    renderLevers(r);
    renderTable(r);
  }

  function baseTheme() {
    return {
      textStyle: { color: "#9aa4bf", fontFamily: "Inter,system-ui,sans-serif" },
      grid: { left: 8, right: 16, top: 24, bottom: 8, containLabel: true },
      tooltip: {
        backgroundColor: "rgba(13,17,32,.96)", borderColor: "rgba(255,255,255,.12)",
        textStyle: { color: "#e8ecf6" }
      }
    };
  }
  function newChart(id) {
    var el = $(id);
    if (!el) return null;
    var old = echarts.getInstanceByDom(el);
    if (old) { old.dispose(); el.innerHTML = ""; }
    var c = echarts.init(el, null, { renderer: "canvas" });
    CHARTS.push(c);
    return c;
  }
  var BANDCOL = { Y: "#c2e0ee", M: "#776b97", L: "#c68dc0" };

  function renderScatter(r) {
    var c = newChart("chartScatter"); if (!c) return;
    var th = baseTheme();
    var series = ["Y", "M", "L"].map(function (b) {
      return {
        name: t("band" + b), type: "scatter", symbolSize: 8,
        itemStyle: { color: BANDCOL[b], opacity: .85 },
        data: r.scatter.filter(function (p) { return p.band === b; }).map(function (p) { return [p.x, p.y]; })
      };
    });
    var opts = th;
    opts.legend = { data: [t("bandY"), t("bandM"), t("bandL")], textStyle: { color: "#9aa4bf" }, top: 0 };
    opts.xAxis = { name: t("scatterX"), nameLocation: "middle", nameGap: 26, type: "value", splitLine: { lineStyle: { color: "rgba(255,255,255,.06)" } }, axisLine: { lineStyle: { color: "rgba(255,255,255,.15)" } } };
    opts.yAxis = { name: t("scatterY"), type: "value", splitLine: { lineStyle: { color: "rgba(255,255,255,.06)" } }, axisLine: { lineStyle: { color: "rgba(255,255,255,.15)" } } };
    // identity line
    var mn = Math.min(r.ageMin, Math.min.apply(null, r.scatter.map(function (p) { return p.y; })));
    var mx = Math.max(r.ageMax, Math.max.apply(null, r.scatter.map(function (p) { return p.y; })));
    series.push({ name: "", type: "line", showSymbol: false, lineStyle: { color: "rgba(255,255,255,.35)", type: "dashed" }, data: [[mn, mn], [mx, mx]], silent: true });
    opts.series = series;
    c.setOption(opts, true);
  }
  function renderViews(r) {
    var c = newChart("chartViews"); if (!c) return;
    var th = baseTheme();
    var labels = [t("viewAll"), t("view40"), t("viewBig"), t("viewEns")];
    var vals = r.views.map(function (v) { return +v.mae.toFixed(2); });
    var colors = vals.map(function (_, i) { return i === 3 ? "#c68dc0" : "#776b97"; });
    th.grid = { left: 8, right: 12, top: 24, bottom: 8, containLabel: true };
    th.xAxis = { type: "category", data: labels, axisLabel: { color: "#9aa4bf", interval: 0, fontSize: 11 }, axisLine: { lineStyle: { color: "rgba(255,255,255,.15)" } }, axisTick: { show: false } };
    th.yAxis = { type: "value", min: 0, splitLine: { lineStyle: { color: "rgba(255,255,255,.06)" } } };
    th.series = [{ type: "bar", data: vals.map(function (v, i) { return { value: v, itemStyle: { color: colors[i], borderRadius: [6, 6, 0, 0] } }; }), barWidth: "52%", label: { show: true, position: "top", color: "#e8ecf6", fontFamily: "monospace", formatter: function (p) { return p.value; } } }];
    th.tooltip = { trigger: "axis" };
    c.setOption(th, true);
  }
  function renderBands(r) {
    var c = newChart("chartBands"); if (!c) return;
    var th = baseTheme();
    var keys = ["Y", "M", "L"];
    var vals = keys.map(function (b) { return isFinite(r.bandMae[b]) ? +r.bandMae[b].toFixed(2) : 0; });
    th.xAxis = { type: "category", data: [t("bandY"), t("bandM"), t("bandL")].map(function (l, i) { return l + (r.bandN[keys[i]] ? " (n=" + r.bandN[keys[i]] + ")" : ""); }), axisLabel: { color: "#9aa4bf", fontSize: 11 }, axisLine: { lineStyle: { color: "rgba(255,255,255,.15)" } }, axisTick: { show: false } };
    th.yAxis = { type: "value", min: 0, splitLine: { lineStyle: { color: "rgba(255,255,255,.06)" } } };
    th.series = [{ type: "bar", data: vals.map(function (v, i) { return { value: v, itemStyle: { color: BANDCOL[keys[i]], borderRadius: [6, 6, 0, 0] } }; }), barWidth: "46%", label: { show: true, position: "top", color: "#e8ecf6", fontFamily: "monospace" } }];
    th.tooltip = { trigger: "axis" };
    c.setOption(th, true);
  }
  function renderLevers(r) {
    var c = newChart("chartLevers"); if (!c) return;
    var th = baseTheme();
    var lev = r.levers.slice(0, 12).reverse();
    th.grid = { left: 8, right: 30, top: 8, bottom: 8, containLabel: true };
    th.xAxis = { type: "value", min: 0, splitLine: { lineStyle: { color: "rgba(255,255,255,.06)" } } };
    th.yAxis = { type: "category", data: lev.map(function (l) { return l.gene; }), axisLabel: { color: "#e8ecf6", fontFamily: "monospace", fontSize: 11 }, axisLine: { lineStyle: { color: "rgba(255,255,255,.15)" } }, axisTick: { show: false } };
    th.series = [{ type: "bar", data: lev.map(function (l) { return { value: +Math.abs(l.w).toFixed(3), itemStyle: { color: l.w < 0 ? "#776b97" : "#c68dc0", borderRadius: l.w < 0 ? [0, 6, 6, 0] : [6, 0, 0, 6] } }; }), barWidth: "58%", label: { show: true, position: "right", color: "#9aa4bf", fontFamily: "monospace", formatter: function (p) { return p.value; } } }];
    th.tooltip = { trigger: "axis", valueFormatter: function (v) { return v; } };
    c.setOption(th, true);
  }
  function renderTable(r) {
    var tb = $("sampleTable").querySelector("tbody");
    tb.innerHTML = "";
    var rows = r.perSample.slice(0, 120);
    var frag = document.createDocumentFragment();
    rows.forEach(function (s) {
      var tr = document.createElement("tr");
      var resid = s.pred - s.age;
      var tr0 = document.createElement("td"); tr0.textContent = s.id;
      var tr1 = document.createElement("td"); tr1.textContent = s.age;
      var tr2 = document.createElement("td"); tr2.textContent = s.pred.toFixed(1);
      var tr3 = document.createElement("td"); tr3.textContent = (resid >= 0 ? "+" : "") + resid.toFixed(1);
      var tr4 = document.createElement("td");
      var sp = document.createElement("span"); sp.className = "band band-" + s.band; sp.textContent = s.band;
      tr4.appendChild(sp);
      tr.appendChild(tr0); tr.appendChild(tr1); tr.appendChild(tr2); tr.appendChild(tr3); tr.appendChild(tr4);
      frag.appendChild(tr);
    });
    tb.appendChild(frag);
  }

  /* ---------- init ---------- */
  function init() {
    setupLang();
    setupDrop();
    $("btnDemo").addEventListener("click", loadDemo);
    window.addEventListener("resize", function () { CHARTS.forEach(function (c) { if (c) c.resize(); }); });
    __applyDom();
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
