/* app.js — UI wiring: language switch, mobile drawer, drag&drop, in-browser compute, charts. */
(function () {
  "use strict";
  var C = window.Clock;
  var CHARTS = [];
  var state = { last: null, wReady: null };

  function tick(ms) { return new Promise(function (res) { setTimeout(res, ms); }); }
  function showBusy(msg) {
    var b = $("busy");
    if (!msg) { b.hidden = true; return; }
    $("busyText").textContent = msg;
    b.hidden = false;
  }
  function loadWeights() {
    // Lite (6.5 KB) lands almost instantly and renders results; the full clock
    // (~1.2 MB) downloads in the background and upgrades the page when ready.
    state.wLiteReady = fetch("data/weights-lite.json")
      .then(function (r) { if (!r.ok) throw new Error("http " + r.status); return r.json(); })
      .then(function (j) { state.lite = j; });
    state.wFullReady = fetch("data/weights.json?v=20261014")
      .then(function (r) { if (!r.ok) throw new Error("http " + r.status); return r.json(); })
      .then(function (j) { state.full = j; });
    state.wLiteReady.catch(function (e) { console.error("weights-lite:", e); });
    state.wFullReady.catch(function (e) { console.error("weights:", e); });
  }
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
      });
    });
    markLang();
    window.__onLang = function () { if (state.last) render(state.last); };
  }

  /* ---------- mobile drawer ---------- */
  function setupDrawer() {
    var menuBtn = $("menuBtn"), overlay = $("overlay");
    if (!menuBtn) return;
    function open() { document.body.classList.add("menu-open"); overlay.hidden = false; menuBtn.setAttribute("aria-expanded", "true"); }
    function close() { document.body.classList.remove("menu-open"); overlay.hidden = true; menuBtn.setAttribute("aria-expanded", "false"); }
    menuBtn.addEventListener("click", open);
    $("drawerClose").addEventListener("click", close);
    overlay.addEventListener("click", close);
    document.addEventListener("keydown", function (e) { if (e.key === "Escape") close(); });
    document.querySelectorAll(".drawer-nav a").forEach(function (a) {
      a.addEventListener("click", function () { close(); });
    });
    window.addEventListener("resize", function () { if (window.innerWidth > 600) close(); });
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
    showStatus(null);
    showBusy(t("busyRead"));
    try {
      var arr = Array.from(files);
      var texts = [];
      for (var i = 0; i < arr.length; i++) {
        var f = arr[i];
        var nm = (f.name || "").toLowerCase();
        if (nm.slice(-4) === ".zip") {
          var packed = await C.unzip(await f.arrayBuffer());
          Object.keys(packed).forEach(function (k) {
            var ext = k.slice(-4).toLowerCase();
            if (ext === ".csv" || ext === ".tsv") texts.push(packed[k]);
          });
        } else {
          texts.push(await f.text());
        }
      }
      var phenTxt = null, phenIdx = -1;
      for (var k = 0; k < texts.length; k++) {
        var ph = C.parsePhenotype(texts[k]);
        if (ph.ok && ph.ids.length > 3) { phenTxt = texts[k]; phenIdx = k; break; }
      }
      if (!phenTxt) { showStatus(t("errParse"), true); return; }
      var phenObj = C.parsePhenotype(phenTxt);
      var exTxt = null;
      for (var m = 0; m < texts.length; m++) {
        if (m === phenIdx) continue;
        var cand = C.parseExpression(texts[m], phenObj.ids);
        if (cand.ok) { exTxt = texts[m]; break; }
      }
      if (!exTxt) { showStatus(t("errParse"), true); return; }
      run(exTxt, phenTxt);
    } catch (e) {
      console.error(e);
      var isZip = e && /zip/.test(String(e.message || e));
      showStatus(t(isZip ? "errZip" : "errUnknown"), true);
    }
  }

  function setupDrop() {
    var dz = $("drop"), inp = $("fileInput");
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
  }

  async function loadDemo() {
    showStatus(null);
    showBusy(t("busyRead"));
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
  async function run(exprText, phenText) {
    showBusy(t("busyRead")); await tick(30);
    var ph = C.parsePhenotype(phenText);
    if (!ph.ok) { showBusy(null); showStatus(t("errParse"), true); return; }
    var ex = C.parseExpression(exprText, ph.ids);
    if (!ex.ok) { showBusy(null); showStatus(t("errParse"), true); return; }
    showBusy(t("busyLite"));
    var wj = null;
    try {
      if (!state.full) await state.wLiteReady;
      wj = state.lite || state.full;
      if (!wj) { await state.wFullReady; wj = state.full; }
      if (!wj) throw new Error("no weights");
      C.loadWeights(wj);
    } catch (e) { showBusy(null); showStatus(t("errWeights"), true); return; }
    showBusy(t("busyRun")); await tick(30);
    var r = C.predict(ex, ph);
    if (!r.ok) { showBusy(null); showStatus(t("errMatch"), true); return; }
    showBusy(null);
    showStatus(null);
    state.last = r;
    state.lastInput = { ex: ex, ph: ph };
    render(r);
    $("results").hidden = false;
    // silent upgrade: full clock arrived -> recompute in place
    state.wFullReady.then(function () {
      if (!state.full || !state.lastInput) return;
      C.loadWeights(state.full);
      var r2 = C.predict(state.lastInput.ex, state.lastInput.ph);
      if (r2.ok && state.last && state.last.meta && state.last.meta.kind === "lite") {
        state.last = r2;
        render(r2);
      }
    }).catch(function () {});
  }

  function bandResidL(r) {
    var L = r.perSample.filter(function (s) { return s.band === "L"; });
    if (!L.length) return NaN;
    var s = 0; L.forEach(function (x) { s += x.pred - x.age; });
    return s / L.length;
  }

  function render(r) {
    var lr = bandResidL(r);
    $("kpiMae").textContent = r.maeEns.toFixed(1);
    $("kpiR2").textContent = r.r2.toFixed(2);
    $("kpiSamples").textContent = String(r.n);
    $("kpiL90").textContent = isFinite(lr) ? Math.abs(lr).toFixed(1) : "–";
    var line = t("smallLine", { rho: r.rho.toFixed(2), prior: fmtPrior(r.meta && r.meta.priorValue) });
    if (r.meta && r.meta.kind === "lite") line += "  ·  " + t("modeLite");
    $("smallLine").textContent = line;
    var mn = t("matchedNote", { m: r.matched, p: r.p });
    var note = $("matchNote");
    note.textContent = mn;
    note.classList.toggle("warn", r.matched < r.p);

    // gate card (plain words: separation score, not "gate AUC")
    var hasGate = isFinite(r.gateAuc) && r.gateN > 0;
    var s = hasGate ? (r.gateAuc >= 0.5 ? r.gateAuc : 1 - r.gateAuc) : NaN;
    $("gateBig").textContent = hasGate ? s.toFixed(2) : "–";
    $("gateScore").textContent = t("gateScore", { s: hasGate ? s.toFixed(2) : "–" });
    $("gateNote").textContent = t("gateNote", { v: isFinite(lr) ? Math.abs(lr).toFixed(1) : "–" });
    $("gateMeta").textContent = t("gateMeta", { n: r.gateN, p: r.p });

    renderScatter(r);
    renderViews(r);
    renderBands(r);
    renderLevers(r);
    renderTable(r);
  }

  /* ---------- charts ---------- */
  function baseTheme() {
    return {
      textStyle: { color: "#a3add2", fontFamily: "Inter,system-ui,sans-serif" },
      tooltip: {
        backgroundColor: "rgba(16,22,48,.96)", borderColor: "rgba(255,255,255,.14)",
        textStyle: { color: "#eef1fb" }
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
  var BANDCOL = { Y: "#9cc7de", M: "#8f83b8", L: "#e29fd4" };

  function renderScatter(r) {
    var c = newChart("chartScatter"); if (!c) return;
    var o = baseTheme();
    var series = ["Y", "M", "L"].map(function (b) {
      return {
        name: t("band" + b), type: "scatter", symbolSize: 8,
        itemStyle: { color: BANDCOL[b], opacity: .85 },
        data: r.scatter.filter(function (p) { return p.band === b; }).map(function (p) { return [p.x, p.y]; })
      };
    });
    o.legend = { data: [t("bandY"), t("bandM"), t("bandL")], textStyle: { color: "#a3add2" }, top: 0 };
    o.grid = { left: 8, right: 16, top: 34, bottom: 8, containLabel: true };
    o.xAxis = { name: t("scatterX"), nameLocation: "middle", nameGap: 28, type: "value", splitLine: { lineStyle: { color: "rgba(255,255,255,.06)" } }, axisLine: { lineStyle: { color: "rgba(255,255,255,.16)" } } };
    o.yAxis = { name: t("scatterY"), type: "value", splitLine: { lineStyle: { color: "rgba(255,255,255,.06)" } }, axisLine: { lineStyle: { color: "rgba(255,255,255,.16)" } } };
    var mn = Math.min(r.ageMin, Math.min.apply(null, r.scatter.map(function (p) { return p.y; })));
    var mx = Math.max(r.ageMax, Math.max.apply(null, r.scatter.map(function (p) { return p.y; })));
    series.push({ name: "", type: "line", showSymbol: false, lineStyle: { color: "rgba(255,255,255,.35)", type: "dashed" }, data: [[mn, mn], [mx, mx]], silent: true });
    o.series = series;
    c.setOption(o, true);
  }
  function fmtPrior(v) { return (typeof v === "number" && isFinite(v)) ? (v >= 0 ? "+" : "") + v.toFixed(2) : "–"; }
  function renderViews(r) {
    var c = newChart("chartViews"); if (!c) return;
    var o = baseTheme();
    var labels = [t("viewDeep"), t("view40"), t("viewAll"), t("viewEns")];
    var vals = r.views.map(function (v) { return +v.mae.toFixed(2); });
    var barCol = ["#776b97", "#9cc7de", "#8f83b8", "#e29fd4"];
    o.grid = { left: 8, right: 12, top: 24, bottom: 8, containLabel: true };
    o.xAxis = { type: "category", data: labels, axisLabel: { color: "#a3add2", interval: 0, fontSize: 11 }, axisLine: { lineStyle: { color: "rgba(255,255,255,.16)" } }, axisTick: { show: false } };
    o.yAxis = { type: "value", min: 0, splitLine: { lineStyle: { color: "rgba(255,255,255,.06)" } } };
    o.series = [{ type: "bar", data: vals.map(function (v, i) { return { value: v, itemStyle: { color: barCol[i], borderRadius: [6, 6, 0, 0] } }; }), barWidth: "52%", label: { show: true, position: "top", color: "#eef1fb", fontFamily: "monospace", formatter: function (p) { return p.value; } } }];
    var mm = r.meta && r.meta.memberMae;
    if (mm && $("svgMaeDeep")) {
      $("svgMaeDeep").textContent = "MAE " + mm.deep;
      $("svgMae40").textContent = "MAE " + mm.en40;
      $("svgMaeAll").textContent = "MAE " + mm.enAll;
    }
    o.tooltip = { trigger: "axis" };
    c.setOption(o, true);
  }
  function renderBands(r) {
    var c = newChart("chartBands"); if (!c) return;
    var o = baseTheme();
    var keys = ["Y", "M", "L"];
    var vals = keys.map(function (b) { return isFinite(r.bandMae[b]) ? +r.bandMae[b].toFixed(2) : 0; });
    o.grid = { left: 8, right: 12, top: 24, bottom: 8, containLabel: true };
    o.xAxis = { type: "category", data: keys.map(function (b, i) { return t("band" + b) + (r.bandN[b] ? " (n=" + r.bandN[b] + ")" : ""); }), axisLabel: { color: "#a3add2", fontSize: 11 }, axisLine: { lineStyle: { color: "rgba(255,255,255,.16)" } }, axisTick: { show: false } };
    o.yAxis = { type: "value", min: 0, splitLine: { lineStyle: { color: "rgba(255,255,255,.06)" } } };
    o.series = [{ type: "bar", data: vals.map(function (v, i) { return { value: v, itemStyle: { color: BANDCOL[keys[i]], borderRadius: [6, 6, 0, 0] } }; }), barWidth: "46%", label: { show: true, position: "top", color: "#eef1fb", fontFamily: "monospace" } }];
    o.tooltip = { trigger: "axis" };
    c.setOption(o, true);
  }
  function renderLevers(r) {
    var c = newChart("chartLevers"); if (!c) return;
    var o = baseTheme();
    var lev = r.levers.slice(0, 12).reverse();
    o.grid = { left: 8, right: 34, top: 8, bottom: 8, containLabel: true };
    o.xAxis = { type: "value", min: 0, splitLine: { lineStyle: { color: "rgba(255,255,255,.06)" } } };
    o.yAxis = { type: "category", data: lev.map(function (l) { return l.gene; }), axisLabel: { color: "#eef1fb", fontFamily: "monospace", fontSize: 11 }, axisLine: { lineStyle: { color: "rgba(255,255,255,.16)" } }, axisTick: { show: false } };
    o.series = [{ type: "bar", data: lev.map(function (l) { return { value: +Math.abs(l.w).toFixed(3), itemStyle: { color: l.w < 0 ? "#8f83b8" : "#e29fd4", borderRadius: l.w < 0 ? [0, 6, 6, 0] : [6, 0, 0, 6] } }; }), barWidth: "58%", label: { show: true, position: "right", color: "#a3add2", fontFamily: "monospace", formatter: function (p) { return p.value; } } }];
    o.tooltip = { trigger: "axis" };
    c.setOption(o, true);
  }
  function renderTable(r) {
    var tb = $("sampleTable").querySelector("tbody");
    tb.innerHTML = "";
    var frag = document.createDocumentFragment();
    r.perSample.slice(0, 120).forEach(function (s) {
      var resid = s.pred - s.age;
      var tr = document.createElement("tr");
      var c0 = document.createElement("td"); c0.textContent = s.id;
      var c1 = document.createElement("td"); c1.textContent = s.age;
      var c2 = document.createElement("td"); c2.textContent = s.pred.toFixed(1);
      var c3 = document.createElement("td"); c3.textContent = (resid >= 0 ? "+" : "") + resid.toFixed(1);
      var c4 = document.createElement("td");
      var sp = document.createElement("span"); sp.className = "band band-" + s.band; sp.textContent = s.band;
      c4.appendChild(sp);
      tr.appendChild(c0); tr.appendChild(c1); tr.appendChild(c2); tr.appendChild(c3); tr.appendChild(c4);
      frag.appendChild(tr);
    });
    tb.appendChild(frag);
  }

  /* ---------- init ---------- */
  function init() {
    setupLang();
    setupDrawer();
    setupDrop();
    loadWeights();
    $("btnDemo").addEventListener("click", loadDemo);
    window.addEventListener("resize", function () { CHARTS.forEach(function (c) { if (c) c.resize(); }); });
    __applyDom();
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init);
  else init();
})();
