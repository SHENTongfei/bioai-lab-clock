# BioAI Lab · LoCRANE Aging Clock (web demo)

Drag in a transcriptomic dataset, get a longevity-certified aging clock — computed
entirely in your browser. No upload, no server: your expression matrix and
phenotype never leave the page.

## What it does

- **Drop zone**: drag in two CSVs — an expression matrix (genes × samples or
  samples × genes, auto-oriented) and a phenotype (`sample_id, age, sex`).
  Or press **Load built-in demo** to use the shipped real-data slice
  (300 genes × 150 samples from the LoCRANE Chinese longevity cohort,
  50 verified 90+ samples).
- **In-browser computation** (5-fold cross-validated):
  - model views: ridge · all genes, EN-style top-40 / top-150 feature budgets,
    equal-weight ensemble (LoCRANE `eq3′` composition analog);
  - locked λ grid `{10, 100, 1000}` — the value with the lowest ensemble CV MAE
    is kept (the protocol lock of the paper, in miniature);
  - KPIs: ensemble MAE, R², Spearman ρ, per-band MAE (Y ≤50 / M 51–89 / L ≥90),
    measured sex-prior value;
  - **90+ longevity gate**: the L-band residual (the 90+ tail reads younger than
    its chronological age — the gate ≠ rate decoupling) plus a dedicated
    5-fold gate classifier AUC;
  - **top causal levers**: largest standardized weights of the all-gene ridge clock.
- **Trilingual UI**: English / 中文 (zh-Hans) / 繁體 (zh-Hant), one click, persisted.
- **Responsive**: desktop and mobile layouts.

## How the clock runs (trained weights, local inference)

The clock is **trained once** (in Python, `tools/export_weights.py`: the LoCRANE
demo recipe — cross-fitted ridge anchor, FiLM deep member x 3 seeds with band
auxiliary + monotonicity, EN-40 / EN-all elastic nets, equal-weight eq3' blend).
The trained weights ship with the page as `data/weights.json` (~4.7 MB, ~1.9 MB
over the wire). When you drop data in, the browser runs **inference only** — a
pure-JS forward pass, bit-faithful to the PyTorch model (verified to <1e-2 on
check vectors). No training happens in the browser, and no file is uploaded.

The page also reports frozen training-time facts (5-split locked CV MAE 6.24 y,
per-member MAE, prior value +0.53 y) next to what it measures on your own data.

The paper's full-cohort model (15 locked splits) runs via
[`code/run_demo.py`](https://github.com/SHENTongfei/LoCRANE/blob/main/code/run_demo.py)
in the [LoCRANE repository](https://github.com/SHENTongfei/LoCRANE).

## Run it

Static site — open `index.html` directly, or serve it:

```
python -m http.server 8000
# → http://localhost:8000
```

Deployed on GitHub Pages: <https://shentongfei.github.io/bioai-lab-clock/>

## Files

| path | role |
|---|---|
| `index.html` | single-page UI |
| `css/style.css` | dark tech / glass theme, responsive |
| `js/engine.js` | CSV/zip parsing + weights inference (UMD, node-testable) |
| `data/weights.json` | trained clock weights (exported by `tools/export_weights.py`) |
| `tools/export_weights.py` | one-shot training + export (production tool, not part of the page) |
| `js/i18n.js` | EN / zh-Hans / zh-Hant dictionaries |
| `js/app.js` | drop zone, language switch, ECharts rendering |
| `data/demo_expression.csv` | real 300-gene × 150-sample slice |
| `data/demo_phenotype.csv` | real `sample_id, age, sex` (50 samples ≥90) |
| `assets/` | BioAI Lab logo + favicon |
