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

## Scope (read before citing)

This is a **demo-scale linear core** of the LoCRANE framework, computed
locally in the browser. The full deep LoCRANE ensemble (residual network with
the four measured priors, the 15-split locked protocol, external transfer,
mining and survival probe) runs via [`code/run_demo.py`](https://github.com/SHENTongfei/LoCRANE/blob/main/code/run_demo.py)
in the [LoCRANE repository](https://github.com/SHENTongfei/LoCRANE). Numbers on
this page are not the paper's frozen results.

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
| `js/engine.js` | data parsing + in-browser clock (UMD, node-testable) |
| `js/i18n.js` | EN / zh-Hans / zh-Hant dictionaries |
| `js/app.js` | drop zone, language switch, ECharts rendering |
| `data/demo_expression.csv` | real 300-gene × 150-sample slice |
| `data/demo_phenotype.csv` | real `sample_id, age, sex` (50 samples ≥90) |
| `assets/` | BioAI Lab logo + favicon |
