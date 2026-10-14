/* engine.js — in-browser aging-clock computation (LoCRANE linear-core demo).
 * Pure JS, zero dependencies. UMD: works in <script> (window.Clock) and in Node
 * (module.exports) so the numerics can be unit-tested before shipping.
 *
 * Honest scope: this is the *linear core* of the clock (ridge regression +
 * elastic-net-style feature budgets), computed 5-fold cross-validated in the
 * browser. It is NOT the full deep ensemble (that runs via code/run_demo.py).
 * The UI labels it exactly that way.
 */
(function (root, factory) {
  if (typeof module === "object" && module.exports) module.exports = factory();
  else root.Clock = factory();
})(typeof self !== "undefined" ? self : this, function () {
  "use strict";

  /* ---------------- CSV parsing ---------------- */
  // Split a CSV text file into rows. Supports comma and tab separators, quoted fields.
  function parseTable(text) {
    text = String(text).replace(/\r\n?/g, "\n");
    var lines = text.split("\n").filter(function (ln) { return ln.length > 0; });
    if (!lines.length) return [];
    // detect separator from the first line
    var first = lines[0];
    var tab = first.split("\t").length > first.split(",").length;
    var sep = tab ? "\t" : ",";
    var rows = [];
    for (var i = 0; i < lines.length; i++) {
      rows.push(splitLine(lines[i], sep));
    }
    return rows;
  }
  function splitLine(line, sep) {
    var out = [], cur = "", inQ = false;
    for (var i = 0; i < line.length; i++) {
      var c = line[i];
      if (inQ) {
        if (c === '"') {
          if (line[i + 1] === '"') { cur += '"'; i++; } else inQ = false;
        } else cur += c;
      } else {
        if (c === '"') inQ = true;
        else if (c === sep) { out.push(cur); cur = ""; }
        else cur += c;
      }
    }
    out.push(cur);
    return out.map(function (s) { return s.trim(); });
  }

  // Parse phenotype: needs sample_id + age (+ optional sex). Tolerant of header names.
  function parsePhenotype(text) {
    var rows = parseTable(text);
    if (rows.length < 2) return { ok: false, error: "Phenotype file is empty." };
    var head = rows[0].map(function (h) { return h.toLowerCase().replace(/[^a-z0-9]/g, ""); });
    var idIx = head.indexOf("sampleid"), ageIx = head.indexOf("age"), sexIx = head.indexOf("sex");
    if (ageIx < 0) ageIx = head.indexOf("chronage");
    if (idIx < 0) idIx = head.indexOf("id");
    if (idIx < 0) idIx = 0;
    var ids = [], age = [], sex = [];
    for (var r = 1; r < rows.length; r++) {
      var row = rows[r];
      var a = parseFloat(row[ageIx]);
      if (!isFinite(a)) continue;
      ids.push(row[idIx]);
      age.push(a);
      var s = (row[sexIx] || "").toUpperCase();
      sex.push(s === "F" ? 1 : (s === "M" ? 0 : 0.5)); // 0.5 = unknown/neutral
    }
    if (age.length < 8) return { ok: false, error: "Need at least 8 samples with a numeric age (got " + age.length + ")." };
    return { ok: true, ids: ids, age: age, sex: sex };
  }

  // Parse an expression matrix. Returns genes[], samples[], mat[genes][samples]
  // and orientation used. Auto-orient: the axis whose labels match the phenotype
  // sample ids is treated as samples.
  function parseExpression(text, phenIds) {
    var rows = parseTable(text);
    if (rows.length < 2) return { ok: false, error: "Expression file is empty." };
    var header = rows[0];
    var hasSampleHeader = header.length > 1 && header[0] === "" &&
      header.slice(1).some(function (c) { return c !== ""; });
    // rows = genes x samples (first col = gene, header = sample ids)
    var matRS = [];
    for (var r = 1; r < rows.length; r++) matRS.push(rows[r].slice(1).map(Number));
    var samplesRS = hasSampleHeader ? header.slice(1) :
      rows[0].slice(1).map(function (c, i) { return c || ("s" + i); });
    var genesRS = rows.slice(1).map(function (row) { return row[0]; });
    // match = how many sample ids line up with phenotype ids (optional; used only to orient)
    var set = {};
    if (phenIds && phenIds.length) phenIds.forEach(function (x) { set[x] = 1; });
    var matchRS = phenIds ? samplesRS.filter(function (s) { return set[s]; }).length : 999;
    // transposed: samples x genes (first col = sample id, first row = gene ids)
    var samplesCS = rows.slice(1).map(function (row) { return row[0]; });
    var genesCS = header.slice(1);
    var matchCS = phenIds ? samplesCS.filter(function (s) { return set[s]; }).length : 0;
    var useCS = matchCS > matchRS;
    var genes = useCS ? genesCS : genesRS;
    var samples = useCS ? samplesCS : samplesRS;
    var mat = useCS
      ? genesCS.map(function (_, gi) { return rows.slice(1).map(function (row) { return Number(row[gi + 1]); }); })
      : matRS;
    // drop non-finite gene rows / sample cols
    var goodCols = samples.map(function (_, ci) {
      return mat.every(function (row) { return row[ci] === undefined || isFinite(row[ci]); });
    });
    mat = mat.map(function (row) { return row.filter(function (_, ci) { return goodCols[ci]; }); });
    samples = samples.filter(function (_, ci) { return goodCols[ci]; });
    mat = mat.filter(function (row) { return row.every(function (v) { return isFinite(v); }); });
    genes = mat.map(function (_, gi) { return genes[gi]; });
    if (genes.length < 10 || samples.length < 8)
      return { ok: false, error: "Matrix too small after cleaning (" + genes.length + " genes × " + samples.length + " samples)." };
    return { ok: true, genes: genes, samples: samples, mat: mat, orientation: useCS ? "samples×genes" : "genes×samples" };
  }

  /* ---------------- linear algebra ---------------- */
  function chol(A) {
    var n = A.length, L = Array.from({ length: n }, function () { return new Array(n).fill(0); });
    for (var i = 0; i < n; i++) {
      for (var j = 0; j <= i; j++) {
        var s = 0;
        for (var k = 0; k < j; k++) s += L[i][k] * L[j][k];
        if (i === j) {
          var v = A[i][i] - s;
          if (v <= 1e-12) v = 1e-12;
          L[i][j] = Math.sqrt(v);
        } else {
          L[i][j] = (A[i][j] - s) / L[j][j];
        }
      }
    }
    return L;
  }
  function solveL(L, b) { // forward substitution: L x = b
    var n = b.length, x = new Array(n);
    for (var i = 0; i < n; i++) {
      var s = b[i];
      for (var k = 0; k < i; k++) s -= L[i][k] * x[k];
      x[i] = s / L[i][i];
    }
    return x;
  }
  function solveLT(L, b) { // back substitution: Lᵀ x = b  (L is lower-triangular)
    var n = b.length, x = new Array(n);
    for (var i = n - 1; i >= 0; i--) {
      var s = b[i];
      for (var k = i + 1; k < n; k++) s -= L[k][i] * x[k]; // Lᵀ[i][k] = L[k][i]
      x[i] = s / L[i][i];
    }
    return x;
  }
  function transpose(X) { // X = rows × cols → cols × rows
    var c = X[0].length, R = X.map(function (r) { return r.length; });
    var T = Array.from({ length: c }, function () { return new Array(R[0] === 0 ? 0 : 0); });
    for (var i = 0; i < X.length; i++) for (var j = 0; j < c; j++) T[j][i] = X[i][j];
    return T;
  }
  function matvec(A, x) {
    var n = A.length, m = A[0].length, y = new Array(m);
    for (var i = 0; i < m; i++) {
      var s = 0;
      for (var j = 0; j < n; j++) s += A[j][i] * x[j]; // A is n×m stored row-major
      y[i] = s;
    }
    return y;
  }
  // X: n samples × p features (rows=samples). y: n targets. lambda: ridge penalty.
  // Returns coefficients (p) + intercept for the centered problem.
  function ridgeFit(X, y, lambda) {
    var n = X.length, p = X[0].length;
    var mx = new Array(p).fill(0), my = 0;
    for (var i = 0; i < n; i++) { my += y[i]; for (var j = 0; j < p; j++) mx[j] += X[i][j]; }
    mx = mx.map(function (v) { return v / n; }); my /= n;
    var C = Array.from({ length: p }, function () { return new Array(p).fill(0); });
    var b = new Array(p).fill(0);
    for (var a = 0; a < n; a++) {
      var xa = X[a], ya = y[a] - my;
      for (var r = 0; r < p; r++) {
        var xr = xa[r] - mx[r];
        b[r] += xr * ya;
        for (var c = r; c < p; c++) C[r][c] += xr * (xa[c] - mx[c]);
      }
    }
    for (var rr = 0; rr < p; rr++) {
      C[rr][rr] += lambda;
      for (var cc = 0; cc < rr; cc++) C[rr][cc] = C[cc][rr];
    }
    var L = chol(C);
    var z = solveL(L, b);          // L z = b
    var beta = solveLT(L, z);      // Lᵀ beta = z  =>  (L Lᵀ) beta = b
    var intercept = my - mx.reduce(function (s, v, j) { return s + v * beta[j]; }, 0);
    return { coef: beta, intercept: intercept, meanX: mx, meanY: my };
  }
  function predict(fit, X) {
    var out = new Array(X.length);
    for (var i = 0; i < X.length; i++) {
      var s = fit.intercept;
      for (var j = 0; j < X[i].length; j++) s += X[i][j] * fit.coef[j];
      out[i] = s;
    }
    return out;
  }

  /* ---------------- metrics ---------------- */
  function mae(pred, act) { var s = 0; for (var i = 0; i < pred.length; i++) s += Math.abs(pred[i] - act[i]); return s / pred.length; }
  function rmse(pred, act) {
    var s = 0; for (var i = 0; i < pred.length; i++) { var d = pred[i] - act[i]; s += d * d; }
    return Math.sqrt(s / pred.length);
  }
  function r2(pred, act) {
    var m = 0; for (var i = 0; i < act.length; i++) m += act[i]; m /= act.length;
    var ss = 0, se = 0;
    for (i = 0; i < act.length; i++) { ss += (act[i] - m) * (act[i] - m); se += (act[i] - pred[i]) * (act[i] - pred[i]); }
    return 1 - se / ss;
  }
  function spearman(a, b) {
    function ranks(x) {
      var idx = x.map(function (v, i) { return [v, i]; }).sort(function (p, q) { return p[0] - q[0]; });
      var r = new Array(x.length);
      for (var i = 0; i < idx.length;) {
        var j = i; while (j < idx.length && idx[j][0] === idx[i][0]) j++;
        var avg = (i + j - 1) / 2 + 1;
        for (var k = i; k < j; k++) r[idx[k][1]] = avg;
        i = j;
      }
      return r;
    }
    var ra = ranks(a), rb = ranks(b), n = a.length;
    var ma = 0, mb = 0; for (var i = 0; i < n; i++) { ma += ra[i]; mb += rb[i]; } ma /= n; mb /= n;
    var num = 0, da = 0, db = 0;
    for (i = 0; i < n; i++) { var a2 = ra[i] - ma, b2 = rb[i] - mb; num += a2 * b2; da += a2 * a2; db += b2 * b2; }
    return num / Math.sqrt(da * db || 1);
  }
  function rocAuc(score, label) {
    var pos = score.filter(function (s, i) { return label[i] === 1; });
    var neg = score.filter(function (s, i) { return label[i] === 0; });
    var nP = pos.length, nN = neg.length;
    if (!nP || !nN) return NaN;
    var all = score.map(function (s, i) { return [s, label[i]]; }).sort(function (a, b) { return b[0] - a[0]; });
    var tp = 0, fp = 0, auc = 0;
    for (var i = 0; i < all.length; i++) {
      if (all[i][1] === 1) { tp++; auc += fp; } else fp++;
    }
    return auc / (nP * nN);
  }

  /* ---------------- the clock ---------------- */
  function buildClock(expr, phen, opts) {
    opts = opts || {};
    var p = expr.genes.length;
    var idAge = {}, idSex = {};
    phen.ids.forEach(function (id, i) { idAge[id] = phen.age[i]; idSex[id] = phen.sex[i]; });
    var keep = [], keepAge = [], keepSex = [], keepBand = [];
    expr.samples.forEach(function (id) {
      if (idAge[id] === undefined) return;
      var a = idAge[id];
      keep.push(expr.samples.indexOf(id));
      keepAge.push(a);
      keepSex.push(idSex[id]);
      keepBand.push(a >= 90 ? "L" : (a <= 50 ? "Y" : "M"));
    });
    if (keep.length < 8) return { ok: false, error: "Fewer than 8 expression samples matched the phenotype." };

    // X rows = samples (n), cols = genes (p); standardized per gene across samples
    var n = keep.length;
    var raw = keep.map(function (si) { return expr.mat.map(function (row) { return row[si]; }); });
    var colMean = new Array(p).fill(0);
    raw.forEach(function (row) { for (var j = 0; j < p; j++) colMean[j] += row[j]; });
    colMean = colMean.map(function (v) { return v / n; });
    var colVar = new Array(p).fill(0);
    raw.forEach(function (row) { for (var j = 0; j < p; j++) colVar[j] += (row[j] - colMean[j]) * (row[j] - colMean[j]); });
    var colSd = colVar.map(function (v, j) { return Math.sqrt(v / n) || 1; });
    var Xg = raw.map(function (row) { return row.map(function (v, j) { return (v - colMean[j]) / colSd[j]; }); });
    var y = keepAge;
    var useSex = opts.useSex !== false;

    // feature budgets by univariate |corr| with age (elastic-net-style top-K)
    var cor = new Array(p).fill(0);
    for (var j = 0; j < p; j++) {
      var sxy = 0;
      for (var i = 0; i < n; i++) sxy += Xg[i][j] * y[i];
      cor[j] = sxy / n;
    }
    var order = cor.map(function (v, j) { return [Math.abs(v), j]; }).sort(function (a, b) { return b[0] - a[0]; });
    var top40 = order.slice(0, Math.min(40, p)).map(function (o) { return o[1]; }).sort(function (a, b) { return a - b; });
    var bigK = Math.min(150, p);
    var topBig = order.slice(0, bigK).map(function (o) { return o[1]; }).sort(function (a, b) { return a - b; });

    function foldOf(i) { return i % 5; }

    // 5-fold CV of a ridge view. Design (D) = n x cols; `cols` = gene indices (or null = all).
    function cvRidge(cols, lambda) {
      var idx = cols || Array.from({ length: D[0].length }, function (_, k) { return k; });
      var Dc = D.map(function (row) { return idx.map(function (k) { return row[k]; }); });
      var pred = new Array(n).fill(0);
      for (var f = 0; f < 5; f++) {
        var tr = [], te = [];
        for (var i = 0; i < n; i++) (foldOf(i) === f ? te : tr).push(i);
        var fit = ridgeFit(tr.map(function (i) { return Dc[i]; }), tr.map(function (i) { return y[i]; }), lambda);
        te.forEach(function (i) { pred[i] = predict(fit, [Dc[i]])[0]; });
      }
      return pred;
    }
    function ensembleMa(Dm, lambda) {
      var a = cvRidge(null, lambda), b = cvRidge(top40, lambda), c = cvRidge(topBig, lambda);
      var ens = a.map(function (v, i) { return (v + b[i] + c[i]) / 3; });
      return { ens: ens, mae: mae(ens, y), a: a, b: b, c: c };
    }

    // design with the sex covariate (if used)
    var D = useSex
      ? Xg.map(function (row, i) { var r = row.slice(); r.push(keepSex[i]); return r; })
      : Xg;

    // locked lambda grid: pick the value with the lowest ensemble CV MAE
    var lamGrid = opts.lamGrid || [10, 100, 1000];
    var best = null;
    lamGrid.forEach(function (lam) {
      var e = ensembleMa(D, lam);
      if (!best || e.mae < best.mae) best = { lam: lam, mae: e.mae, e: e };
    });
    var lambda = best.lam;
    var ens = best.e;
    var predAll = ens.a, pred40 = ens.b, predBig = ens.c, predEns = ens.ens;

    var maeAll = mae(predAll, y), mae40 = mae(pred40, y), maeBig = mae(predBig, y), maeEns = best.mae;
    var r2Ens = r2(predEns, y), rhoEns = spearman(predEns, y);

    function bandMae(b) {
      var ys = [], ps = [];
      keepBand.forEach(function (bb, i) { if (bb === b) { ys.push(y[i]); ps.push(predEns[i]); } });
      return ys.length ? mae(ps, ys) : NaN;
    }
    var bandMaeOut = { Y: bandMae("Y"), M: bandMae("M"), L: bandMae("L") };
    var bandN = { Y: 0, M: 0, L: 0 };
    keepBand.forEach(function (b) { bandN[b]++; });

    // sex-prior value: ensemble MAE without the sex covariate vs with it
    var sexValue = 0;
    if (useSex) {
      var eNo = ensembleMa(Xg, lambda);
      sexValue = eNo.mae - maeEns; // >0 => the sex covariate lowers MAE
    }

    // gate: dedicated 5-fold ridge classifier for 90+ (L) vs the rest -> AUC
    var lab = keepBand.map(function (b) { return b === "L" ? 1 : 0; });
    var hasPos = lab.some(function (v) { return v === 1; });
    var gateAuc = NaN;
    if (hasPos && lab.some(function (v) { return v === 0; })) {
      var idxAll = Array.from({ length: D[0].length }, function (_, k) { return k; });
      var Dc = D.map(function (row) { return idxAll.map(function (k) { return row[k]; }); });
      var gpred = new Array(n).fill(0);
      for (var f = 0; f < 5; f++) {
        var tr = [], te = [];
        for (var i = 0; i < n; i++) (foldOf(i) === f ? te : tr).push(i);
        var fit = ridgeFit(tr.map(function (i) { return Dc[i]; }), tr.map(function (i) { return lab[i]; }), lambda);
        te.forEach(function (i) { gpred[i] = predict(fit, [Dc[i]])[0]; });
      }
      gateAuc = rocAuc(gpred, lab);
    }

    // levers: top |standardized coef| of the all-data, all-gene ridge
    var fitAll = ridgeFit(Xg, y, lambda);
    var levers = fitAll.coef.map(function (c, jj) { return [Math.abs(c), jj]; })
      .sort(function (a, b) { return b[0] - a[0]; })
      .slice(0, 15)
      .map(function (o) { return { gene: expr.genes[o[1]], w: fitAll.coef[o[1]] }; })
      .sort(function (a, b) { return Math.abs(b.w) - Math.abs(a.w); });

    return {
      ok: true,
      n: n, p: p, lambda: lambda,
      ageMin: Math.min.apply(null, y), ageMax: Math.max.apply(null, y),
      sexValue: sexValue,
      views: [
        { key: "all", labelAll: true, mae: maeAll },
        { key: "40", mae: mae40, k: top40.length },
        { key: "big", mae: maeBig, k: bigK },
        { key: "ens", mae: maeEns }
      ],
      maeEns: maeEns, r2: r2Ens, rho: rhoEns,
      bandMae: bandMaeOut, bandN: bandN,
      gateAuc: gateAuc, gateN: lab.filter(function (v) { return v === 1; }).length,
      levers: levers,
      scatter: keepBand.map(function (b, i) { return { x: y[i], y: predEns[i], band: b }; }),
      perSample: keep.map(function (id, i) { return { id: id, age: y[i], pred: predEns[i], band: keepBand[i] }; })
    };
  }

  return {
    parsePhenotype: parsePhenotype,
    parseExpression: parseExpression,
    buildClock: buildClock,
    _metrics: { mae: mae, rmse: rmse, r2: r2, spearman: spearman, rocAuc: rocAuc }
  };
});
