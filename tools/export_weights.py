# -*- coding: utf-8 -*-
"""export_weights.py — train the LoCRANE demo-scale clock once, export weights.

Same recipe as LoCRANE_repo/code/run_demo.py (cross-fitted ridge anchor, FiLM deep
member x 3 seeds with band auxiliary + monotonicity, EN-40 / EN-all elastic nets,
equal-weight eq3' blend). Two things happen here:

  1. the locked 5-split CV is re-run to freeze honest quality numbers (MAE per
     member, prior value) that the web page displays as "trained" facts;
  2. the FINAL model is retrained on the full released slice and exported to
     data/weights.json — the browser only ever does inference with these weights
     (pure-JS forward pass in js/engine.js; data never leaves the device).

Also writes data/weights_check.json: per-sample intermediate values of the final
model so the Node test can verify the JS forward pass against PyTorch to 1e-3.

Run:  python tools/export_weights.py     (needs torch/sklearn/scipy/pandas)
"""
import gzip, json, os, sys
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.linear_model import ElasticNet, Ridge
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "data")
SRC = sys.argv[1] if len(sys.argv) > 1 else r"C:\Users\TS\Desktop\crane\LoCRANE_repo\data\demo"
LAM_MON, LAM_SP = 0.05, 0.05
SEEDS = (42, 7, 2024)
EPOCHS = 60
SPLITS = 5


class PriorsMLP(nn.Module):
    """Residual age network with the four LoCRANE priors (identical to run_demo)."""

    def __init__(self, d_in, use_film=True, use_band=True):
        super().__init__()
        self.use_film, self.use_band = use_film, use_band
        self.fc1 = nn.Linear(d_in, 128)
        self.block = nn.Sequential(nn.Linear(128, 128), nn.ReLU(), nn.Linear(128, 128), nn.ReLU())
        self.head = nn.Linear(128, 1)
        self.band = nn.Linear(128, 3)
        self.film = nn.Sequential(nn.Linear(1, 16), nn.ReLU(), nn.Linear(16, 256))

    def forward(self, x, sex):
        h = F.relu(self.fc1(x))
        if self.use_film:
            gb = self.film(sex.view(-1, 1))
            gamma, beta = gb[:, :128], gb[:, 128:]
            h = h * (1 + torch.tanh(gamma)) + torch.tanh(beta)
        h = h + self.block(h)
        return self.head(h).squeeze(-1), self.band(h)


def band_of(a):
    return np.where(a <= 50, 0, np.where(a < 90, 1, 2)).astype(np.int64)


def train_deep(Xtr, ytr, sex_tr, ap_tr, s, Xte, sex_te, ap_te,
               yL, XL, sexL, use_film, use_band, use_mono, seed, epochs=EPOCHS):
    """run_demo.train_deep, verbatim: net outputs a tanh-bounded correction on the
    anchor; final prediction = anchor + s * tanh(g(x))."""
    torch.manual_seed(seed); np.random.seed(seed)
    net = PriorsMLP(Xtr.shape[1], use_film, use_band)
    opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-3)
    Xt = torch.tensor(Xtr, dtype=torch.float32); yt = torch.tensor(ytr, dtype=torch.float32)
    st = torch.tensor(sex_tr, dtype=torch.float32)
    at = torch.tensor((ytr - ap_tr) / s, dtype=torch.float32)
    bt = torch.tensor(band_of(ytr), dtype=torch.int64)
    XLt = torch.tensor(XL, dtype=torch.float32); sL = torch.tensor(sexL, dtype=torch.float32)
    bL = torch.tensor(band_of(yL), dtype=torch.int64)
    Xe = torch.tensor(Xte, dtype=torch.float32)
    n = len(yt)
    for ep in range(epochs):
        net.train()
        for i in torch.split(torch.randperm(n), 64):
            corr, band_logits = net(Xt[i], st[i])
            g = torch.tanh(corr)
            pred = at[i] + g
            loss = F.smooth_l1_loss(pred, torch.zeros_like(pred))
            if use_band:
                loss = loss + LAM_SP * F.cross_entropy(band_logits, bt[i])
            if use_mono and i.numel() > 1:
                d = pred.view(-1, 1) - pred.view(1, -1)
                sign = (yt[i].view(-1, 1) - yt[i].view(1, -1)) < 0
                loss = loss + LAM_MON * (F.softplus(d) * sign.float()).mean()
            opt.zero_grad(); loss.backward(); opt.step()
        if use_band and len(yL):
            net.train()
            _, bl = net(XLt, sL)
            loss = 0.3 * F.cross_entropy(bl, bL)
            opt.zero_grad(); loss.backward(); opt.step()
    net.eval()
    with torch.no_grad():
        corr, _ = net(Xe, torch.tensor(sex_te, dtype=torch.float32))
        return ap_te + s * torch.tanh(corr).numpy()


def load_slice():
    with gzip.open(os.path.join(SRC, "expression_demo.csv.gz"), "rt") as f:
        E = pd.read_csv(f, sep="\t", index_col=0)
    P = pd.read_csv(os.path.join(SRC, "phenotype_demo.csv")).set_index("sample_id")
    common = [s for s in P.index if s in E.columns]
    E, P = E[common], P.loc[common]
    genes = list(E.index)
    X = np.nan_to_num(E.T.values, nan=0.0).astype(np.float32)
    y = P["age"].to_numpy(float)
    sex = (P["sex"].astype(str).str.upper().str[0] == "F").to_numpy().astype(np.float32)
    return genes, X, y, sex, list(P.index)


def frozen_cv(X, y, sex, is_L):
    """The locked 5-split CV of run_demo, with per-member errors kept."""
    rng = np.random.default_rng(2026)
    SPL = [rng.permutation(np.where(~is_L)[0]) for _ in range(SPLITS)]
    err = {"eq3p": [], "deep": [], "en40": [], "enAll": [], "nopriors": [], "ridge": []}
    for r, perm in enumerate(SPL):
        k = max(1, len(perm) // 5)
        te, tr = perm[:k], perm[k:]
        trL = np.where(is_L)[0]
        sc = StandardScaler().fit(X[tr]); Xs = sc.transform(X).astype(np.float32)
        r0 = Ridge(alpha=10.0).fit(Xs[tr], y[tr])
        ap_kf = np.zeros(len(tr), dtype=np.float32)
        for a, b in KFold(5, shuffle=True, random_state=42).split(tr):
            r_ = Ridge(alpha=10.0).fit(Xs[tr[a]], y[tr[a]]); ap_kf[b] = r_.predict(Xs[tr[b]])
        s = float(np.std(y[tr] - ap_kf)) + 1e-6
        ap_te = r0.predict(Xs[te])
        ds = [train_deep(Xs[tr], y[tr], sex[tr], ap_kf, s, Xs[te], sex[te], ap_te,
                         y[trL], Xs[trL], sex[trL], True, True, True, sd) for sd in SEEDS]
        pd_ = np.mean(ds, axis=0)
        fc = np.array([abs(np.corrcoef(Xs[tr, j], y[tr])[0, 1]) if Xs[tr, j].std() > 0 else 0
                       for j in range(Xs.shape[1])])
        top40, top500 = np.argsort(-fc)[:40], np.argsort(-fc)[:500]
        e40 = ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=20000).fit(Xs[tr][:, top40], y[tr]).predict(Xs[te][:, top40])
        eA = ElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=20000).fit(Xs[tr][:, top500], y[tr]).predict(Xs[te][:, top500])
        eq = np.mean([pd_, e40, eA], axis=0)
        pn = np.mean([train_deep(Xs[tr], y[tr], sex[tr], ap_kf, s, Xs[te], sex[te], ap_te,
                                 y[trL], Xs[trL], sex[trL], False, False, False, sd) for sd in SEEDS], axis=0)
        rg = Ridge(alpha=10.0).fit(Xs[tr], y[tr]).predict(Xs[te])
        err["eq3p"].append(float(np.abs(eq - y[te]).mean()))
        err["deep"].append(float(np.abs(pd_ - y[te]).mean()))
        err["en40"].append(float(np.abs(e40 - y[te]).mean()))
        err["enAll"].append(float(np.abs(eA - y[te]).mean()))
        err["nopriors"].append(float(np.abs(pn - y[te]).mean()))
        err["ridge"].append(float(np.abs(rg - y[te]).mean()))
        print(f"[cv split {r}] eq3' {err['eq3p'][-1]:.2f} | deep {err['deep'][-1]:.2f} "
              f"| en40 {err['en40'][-1]:.2f} | enAll {err['enAll'][-1]:.2f} | priors-off {err['nopriors'][-1]:.2f}")
    return {k: float(np.mean(v)) for k, v in err.items()}


def net_to_list(net):
    """Export one seed's parameters with PyTorch [out,in] row layout (JS reads rows)."""
    def t(x, r=6):
        return np.round(x.detach().numpy(), r).tolist()
    return {
        "fc1W": t(net.fc1.weight), "fc1B": t(net.fc1.bias),
        "blk1W": t(net.block[0].weight), "blk1B": t(net.block[0].bias),
        "blk2W": t(net.block[2].weight), "blk2B": t(net.block[2].bias),
        "headW": t(net.head.weight), "headB": t(net.head.bias),
        "bandW": t(net.band.weight), "bandB": t(net.band.bias),
        "film1W": t(net.film[0].weight), "film1B": t(net.film[0].bias),
        "film2W": t(net.film[2].weight), "film2B": t(net.film[2].bias),
    }


def pack_arrays(obj):
    """Recursively pack long numeric lists as base64 little-endian float32.
    Cuts the download roughly in half and lets the browser skip JSON parsing."""
    import base64, struct
    if isinstance(obj, dict):
        return {k: pack_arrays(v) for k, v in obj.items()}
    if isinstance(obj, list) and len(obj) >= 16 and all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in obj):
        return {"__f32": base64.b64encode(struct.pack("<%df" % len(obj), *obj)).decode("ascii")}
    if isinstance(obj, list):
        return [pack_arrays(v) for v in obj]
    return obj


def main():
    genes, X, y, sex, ids = load_slice()
    n, p = X.shape
    is_L = y >= 90
    print(f"[data] n={n} (working-age {int((~is_L).sum())}, 90+ {int(is_L.sum())}) genes={p}")

    print("[1/3] frozen 5-split CV (honest numbers for the page)...")
    cv = frozen_cv(X, y, sex, is_L)

    print("[2/3] final model on the full slice...")
    scF = StandardScaler().fit(X)
    XsF = scF.transform(X).astype(np.float32)
    rF = Ridge(alpha=10.0).fit(XsF, y)
    ap_kf = np.zeros(n, dtype=np.float32)
    for a, b in KFold(5, shuffle=True, random_state=42).split(np.arange(n)):
        r_ = Ridge(alpha=10.0).fit(XsF[a], y[a]); ap_kf[b] = r_.predict(XsF[b])
    sF = float(np.std(y - ap_kf)) + 1e-6
    trL = np.where(is_L)[0]
    nets = []
    for sd in SEEDS:
        net = PriorsMLP(p, True, True)
        torch.manual_seed(sd); np.random.seed(sd)
        opt = torch.optim.Adam(net.parameters(), lr=1e-3, weight_decay=1e-3)
        Xt = torch.tensor(XsF, dtype=torch.float32); yt = torch.tensor(y, dtype=torch.float32)
        st = torch.tensor(sex, dtype=torch.float32)
        at = torch.tensor((y - ap_kf) / sF, dtype=torch.float32)
        bt = torch.tensor(band_of(y), dtype=torch.int64)
        XLt = torch.tensor(XsF[trL], dtype=torch.float32); sL = torch.tensor(sex[trL], dtype=torch.float32)
        bL = torch.tensor(band_of(y[trL]), dtype=torch.int64)
        for ep in range(EPOCHS):
            net.train()
            for i in torch.split(torch.randperm(n), 64):
                corr, bl = net(Xt[i], st[i])
                g = torch.tanh(corr)
                pred = at[i] + g
                loss = F.smooth_l1_loss(pred, torch.zeros_like(pred))
                loss = loss + LAM_SP * F.cross_entropy(bl, bt[i])
                d = pred.view(-1, 1) - pred.view(1, -1)
                sign = (yt[i].view(-1, 1) - yt[i].view(1, -1)) < 0
                loss = loss + LAM_MON * (F.softplus(d) * sign.float()).mean()
                opt.zero_grad(); loss.backward(); opt.step()
            net.train()
            _, blL = net(XLt, sL)
            loss = 0.3 * F.cross_entropy(blL, bL)
            opt.zero_grad(); loss.backward(); opt.step()
        net.eval()
        nets.append(net)
        print(f"  seed {sd} trained")

    # final-model inference on the full slice (check vectors)
    with torch.no_grad():
        Xt = torch.tensor(XsF, dtype=torch.float32); st = torch.tensor(sex, dtype=torch.float32)
        corrs, bands = [], []
        for net in nets:
            c_, b_ = net(Xt, st)
            corrs.append(c_.numpy()); bands.append(b_.numpy())
    corr_m = np.mean([np.tanh(c) for c in corrs], axis=0)   # tanh per seed, then average (run_demo semantics)
    band_m = np.mean(bands, axis=0)            # (n,3)
    deep_all = rF.predict(XsF) + sF * corr_m
    fcF = np.array([abs(np.corrcoef(XsF[:, j], y)[0, 1]) if XsF[:, j].std() > 0 else 0 for j in range(p)])
    top40F = np.argsort(-fcF)[:40]
    en40F = ElasticNet(alpha=0.1, l1_ratio=0.5, max_iter=20000).fit(XsF[:, top40F], y)
    enAllF = ElasticNet(alpha=0.01, l1_ratio=0.5, max_iter=20000).fit(XsF, y)
    e40_all = en40F.predict(XsF[:, top40F])
    eA_all = enAllF.predict(XsF)
    ens_all = np.mean([deep_all, e40_all, eA_all], axis=0)

    print("[3/3] export...")
    weights = {
        "version": 1,
        "genes": genes,
        "geneMean": np.round(scF.mean_, 5).tolist(),
        "geneSd": np.round(scF.scale_, 5).tolist(),
        "residSd": round(sF, 6),
        "anchor": {"coef": np.round(rF.coef_, 6).tolist(), "intercept": round(float(rF.intercept_), 6)},
        "seeds": [net_to_list(net) for net in nets],
        "en40": {"idx": top40F.tolist(), "coef": np.round(en40F.coef_, 6).tolist(),
                 "intercept": round(float(en40F.intercept_), 6)},
        "enAll": {"coef": np.round(enAllF.coef_, 6).tolist(), "intercept": round(float(enAllF.intercept_), 6)},
        "bandNames": ["Y", "M", "L"],
        "meta": {
            "protocol": "demo-scale locked 5-split CV; final model retrained on the full released slice",
            "frozenMae": round(cv["eq3p"], 3),
            "memberMae": {"deep": round(cv["deep"], 3), "en40": round(cv["en40"], 3), "enAll": round(cv["enAll"], 3)},
            "priorValue": round(cv["nopriors"] - cv["eq3p"], 3),
            "ridgeRef": round(cv["ridge"], 3),
            "n": int(n), "genes": int(p), "epochs": EPOCHS, "seeds": list(SEEDS),
            "source": "LoCRANE released slice (100 working-age + 50 verified 90+)",
        },
    }
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "weights.json"), "w") as f:
        json.dump(pack_arrays(weights), f, separators=(",", ":"))
    kb = os.path.getsize(os.path.join(OUT, "weights.json")) / 1024

    check = []
    for i in range(min(8, n)):
        check.append({
            "id": ids[i], "sex": float(sex[i]), "age": float(y[i]),
            "anchor": float(rF.predict(XsF[i:i + 1])[0]),
            "corr": [float(c[i]) for c in corrs],
            "deep": float(deep_all[i]), "en40": float(e40_all[i]), "enAll": float(eA_all[i]),
            "ens": float(ens_all[i]),
            "bandL": float(band_m[i][2]),
        })
    with open(os.path.join(OUT, "weights_check.json"), "w") as f:
        json.dump(check, f, indent=1)

    print(f"\nweights.json  {kb:.0f} KB   check vectors: {len(check)}")
    print(f"frozen eq3' MAE {cv['eq3p']:.3f} y | deep {cv['deep']:.3f} | en40 {cv['en40']:.3f} "
          f"| enAll {cv['enAll']:.3f} | prior value {cv['nopriors'] - cv['eq3p']:+.3f} y")
    print(f"full-slice in-sample ens MAE {np.abs(ens_all - y).mean():.3f} y (reference only)")


if __name__ == "__main__":
    main()
