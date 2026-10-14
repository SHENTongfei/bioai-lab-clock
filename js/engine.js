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

  /* ---------------- pretrained clock (weights + pure-JS inference) ---------------- */
  // The clock is trained once (tools/export_weights.py, locked demo protocol) and
  // shipped as data/weights.json. Inference below is a faithful pure-JS port of the
  // PyTorch forward pass: standardized inputs -> ridge anchor -> residual deep
  // member (FiLM sex conditioning, band head) x 3 seeds -> EN-40 / EN-all ->
  // equal-weight blend. Nothing is learned in the browser.
  var W = null;
  function loadWeights(obj) { W = obj; return true; }
  function hasWeights() { return !!W; }

  function dot(a, b) {
    var s = 0;
    for (var i = 0; i < a.length; i++) s += a[i] * b[i];
    return s;
  }

  // one deep seed: returns {corr, band[3]} for a standardized row x and sex code
  function deepForward(x, sexv, sd, p) {
    var h = new Array(128), u, j;
    for (u = 0; u < 128; u++) {
      var s0 = sd.fc1B[u] + dot(sd.fc1W[u], x);
      h[u] = s0 > 0 ? s0 : 0;
    }
    // FiLM: sex (scalar) -> Linear(1,16) -> ReLU -> Linear(16,256) (gamma 128 | beta 128)
    var f1 = new Array(16), a;
    for (a = 0; a < 16; a++) {
      var fv = sd.film1B[a] + sexv * sd.film1W[a][0];
      f1[a] = fv > 0 ? fv : 0;
    }
    var gb = new Array(256);
    for (var b2 = 0; b2 < 256; b2++) {
      var sv = sd.film2B[b2] + dot(sd.film2W[b2], f1);
      gb[b2] = Math.tanh(sv);
    }
    for (u = 0; u < 128; u++) h[u] = h[u] * (1 + gb[u]) + gb[128 + u];
    // residual block
    var t1 = new Array(128);
    for (u = 0; u < 128; u++) {
      var s3 = sd.blk1B[u] + dot(sd.blk1W[u], h);
      t1[u] = s3 > 0 ? s3 : 0;
    }
    for (u = 0; u < 128; u++) {
      var o2 = sd.blk2B[u] + dot(sd.blk2W[u], t1);
      h[u] += o2 > 0 ? o2 : 0;  // the block ends with a ReLU on its own output (training semantics)
    }
    var corr = sd.headB[0] + dot(sd.headW[0], h);
    var bl = [0, 0, 0];
    for (var c2 = 0; c2 < 3; c2++) bl[c2] = sd.bandB[c2] + dot(sd.bandW[c2], h);
    return { corr: corr, band: bl };
  }

  function softmax3(z) {
    var m = Math.max(z[0], z[1], z[2]);
    var e = z.map(function (v) { return Math.exp(v - m); });
    var t = e[0] + e[1] + e[2];
    return [e[0] / t, e[1] / t, e[2] / t];
  }

  function predict(expr, phen) {
    if (!W) return { ok: false, error: "weights not loaded" };
    var p = W.genes.length;
    var idAge = {}, idSex = {};
    phen.ids.forEach(function (id, i) { idAge[id] = phen.age[i]; idSex[id] = phen.sex[i]; });
    var ids = [], age = [], sex = [], band = [];
    expr.samples.forEach(function (id) {
      if (idAge[id] === undefined) return;
      var a = idAge[id];
      ids.push(id); age.push(a); sex.push(idSex[id]);
      band.push(a >= 90 ? "L" : (a <= 50 ? "Y" : "M"));
    });
    if (ids.length < 1) return { ok: false, error: "none of the expression samples matched the sample list" };
    var n = ids.length;

    // align the user's matrix to the trained gene order (missing gene -> z 0 = mean)
    var colOf = {};
    expr.genes.forEach(function (g, j) { if (colOf[g] === undefined) colOf[g] = j; });
    var matched = 0, j;
    for (j = 0; j < p; j++) if (colOf[W.genes[j]] !== undefined) matched++;
    var Xs = new Array(n);
    for (var i = 0; i < n; i++) {
      var si = expr.samples.indexOf(ids[i]);
      var row = new Array(p);
      for (j = 0; j < p; j++) {
        var c = colOf[W.genes[j]];
        if (c === undefined) { row[j] = 0; continue; }
        var v = (expr.mat[j][si] - W.geneMean[j]) / W.geneSd[j];
        row[j] = isFinite(v) ? v : 0;
      }
      Xs[i] = row;
    }

    // ridge anchor
    var ap = Xs.map(function (x) { return W.anchor.intercept + dot(W.anchor.coef, x); });

    // deep member x seeds (anchor + residSd * tanh(corr)), band logits averaged
    var deepP = new Array(n), pL = new Array(n);
    for (i = 0; i < n; i++) {
      var acc = 0, bl = [0, 0, 0];
      for (var k = 0; k < W.seeds.length; k++) {
        var f = deepForward(Xs[i], sex[i], W.seeds[k], p);
        acc += Math.tanh(f.corr);
        for (var c3 = 0; c3 < 3; c3++) bl[c3] += f.band[c3];
      }
      deepP[i] = ap[i] + W.residSd * (acc / W.seeds.length);
      var pr = softmax3(bl.map(function (v) { return v / W.seeds.length; }));
      pL[i] = pr[2];
    }

    // linear members
    var en40P = Xs.map(function (x) {
      var s = W.en40.intercept;
      for (var q = 0; q < W.en40.idx.length; q++) s += W.en40.coef[q] * x[W.en40.idx[q]];
      return s;
    });
    var enAllP = Xs.map(function (x) { return W.enAll.intercept + dot(W.enAll.coef, x); });
    var ensP = deepP.map(function (v, i2) { return (v + en40P[i2] + enAllP[i2]) / 3; });

    var maeD = mae(deepP, age), mae40 = mae(en40P, age), maeA = mae(enAllP, age), maeE = mae(ensP, age);

    function bandMae(b) {
      var ys = [], ps = [];
      band.forEach(function (bb, i2) { if (bb === b) { ys.push(age[i2]); ps.push(ensP[i2]); } });
      return ys.length ? mae(ps, ys) : NaN;
    }
    var bandN = { Y: 0, M: 0, L: 0 };
    band.forEach(function (b) { bandN[b]++; });

    var lab = band.map(function (b) { return b === "L" ? 1 : 0; });
    var hasPos = lab.indexOf(1) >= 0, hasNeg = lab.indexOf(0) >= 0;
    var gateAuc = (hasPos && hasNeg) ? rocAuc(pL, lab) : NaN;

    // frozen levers: strongest elastic-net weights (fixed by training)
    var src = W.en40.coef.some(function (v) { return v !== 0; }) ? W.en40 : W.enAll;
    var levers = src.coef.map(function (w, q) {
      var gi = (src === W.en40) ? W.en40.idx[q] : q;
      return { gene: W.genes[gi], w: w };
    }).filter(function (o) { return o.w !== 0; })
      .sort(function (a, b) { return Math.abs(b.w) - Math.abs(a.w); })
      .slice(0, 15);

    return {
      ok: true,
      mode: "weights",
      n: n, p: p, matched: matched,
      ageMin: Math.min.apply(null, age), ageMax: Math.max.apply(null, age),
      views: [
        { key: "deep", mae: maeD },
        { key: "en40", mae: mae40 },
        { key: "enAll", mae: maeA },
        { key: "ens", mae: maeE }
      ],
      maeEns: maeE, r2: r2(ensP, age), rho: spearman(ensP, age),
      bandMae: { Y: bandMae("Y"), M: bandMae("M"), L: bandMae("L") }, bandN: bandN,
      gateAuc: gateAuc, gateN: bandN.L,
      levers: levers,
      scatter: band.map(function (b, i2) { return { x: age[i2], y: ensP[i2], band: b }; }),
      perSample: ids.map(function (id, i2) {
        return { id: id, age: age[i2], pred: ensP[i2], band: band[i2], lProb: pL[i2] };
      }),
      meta: W.meta || {}
    };
  }

  /* ---------------- zip (demo bundle) ---------------- */
  // Minimal zip reader: STORE (0) and DEFLATE (8) entries, no zip64.
  function unzip(buffer) {
    return new Promise(function (resolve, reject) {
      try {
        var u8 = new Uint8Array(buffer);
        var dv = new DataView(buffer);
        var eocd = -1;
        for (var i = u8.length - 22; i >= Math.max(0, u8.length - 65558); i--) {
          if (u8[i] === 0x50 && u8[i + 1] === 0x4b && u8[i + 2] === 0x05 && u8[i + 3] === 0x06) { eocd = i; break; }
        }
        if (eocd < 0) { reject(new Error("not a zip")); return; }
        var count = dv.getUint16(eocd + 10, true);
        var ptr = dv.getUint32(eocd + 16, true);
        var jobs = [];
        for (var k = 0; k < count; k++) {
          var method = dv.getUint16(ptr + 10, true);
          var csize = dv.getUint32(ptr + 20, true);
          var nlen = dv.getUint16(ptr + 28, true);
          var elen = dv.getUint16(ptr + 30, true);
          var clen = dv.getUint16(ptr + 32, true);
          var lho = dv.getUint32(ptr + 42, true);
          var name = new TextDecoder().decode(u8.slice(ptr + 46, ptr + 46 + nlen));
          var nlen2 = dv.getUint16(lho + 26, true);
          var elen2 = dv.getUint16(lho + 28, true);
          var start = lho + 30 + nlen2 + elen2;
          jobs.push({ name: name, method: method, data: u8.slice(start, start + csize) });
          ptr += 46 + nlen + elen + clen;
        }
        Promise.all(jobs.map(function (j) {
          if (j.method === 0) return Promise.resolve(new TextDecoder().decode(j.data));
          if (j.method !== 8) return Promise.reject(new Error("unsupported zip compression"));
          var ds = new DecompressionStream("deflate-raw");
          var writer = ds.writable.getWriter();
          writer.write(j.data);
          writer.close();
          var reader = ds.readable.getReader();
          var chunks = [];
          function read() {
            return reader.read().then(function (r) {
              if (r.done) {
                var total = 0; chunks.forEach(function (c) { total += c.length; });
                var out = new Uint8Array(total); var off = 0;
                chunks.forEach(function (c) { out.set(c, off); off += c.length; });
                return new TextDecoder().decode(out);
              }
              chunks.push(r.value);
              return read();
            });
          }
          return read();
        })).then(function (texts) {
          var files = {};
          texts.forEach(function (t, idx) { files[jobs[idx].name] = t; });
          resolve(files);
        }, reject);
      } catch (e) { reject(e); }
    });
  }

  return {
    parsePhenotype: parsePhenotype,
    parseExpression: parseExpression,
    loadWeights: loadWeights,
    hasWeights: hasWeights,
    predict: predict,
    unzip: unzip,
    _metrics: { mae: mae, rmse: rmse, r2: r2, spearman: spearman, rocAuc: rocAuc }
  };
});
