# -*- coding: utf-8 -*-
"""diag_forward.py — dump per-layer intermediates of the final model for sample 0
(full precision) so the Node side can locate the first divergence."""
import gzip, json, os, sys
import numpy as np
import pandas as pd
import torch
from sklearn.preprocessing import StandardScaler

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = sys.argv[1] if len(sys.argv) > 1 else r"C:\Users\TS\Desktop\crane\LoCRANE_repo\data\demo"

with gzip.open(os.path.join(SRC, "expression_demo.csv.gz"), "rt") as f:
    E = pd.read_csv(f, sep="\t", index_col=0)
P = pd.read_csv(os.path.join(SRC, "phenotype_demo.csv")).set_index("sample_id")
common = [s for s in P.index if s in E.columns]
E, P = E[common], P.loc[common]
X = np.nan_to_num(E.T.values, nan=0.0).astype(np.float32)
sex = (P["sex"].astype(str).str.upper().str[0] == "F").to_numpy().astype(np.float32)

W = json.load(open(os.path.join(ROOT, "data", "weights.json")))
# rebuild the trained nets from the exported (rounded) weights to mirror the JS side
import torch.nn as nn
class PriorsMLP(nn.Module):
    def __init__(self, d_in):
        super().__init__()
        self.fc1 = nn.Linear(d_in, 128)
        self.block = nn.Sequential(nn.Linear(128, 128), nn.ReLU(), nn.Linear(128, 128), nn.ReLU())
        self.head = nn.Linear(128, 1)
        self.band = nn.Linear(128, 3)
        self.film = nn.Sequential(nn.Linear(1, 16), nn.ReLU(), nn.Linear(16, 256))
    def forward(self, x, sex):
        h = torch.relu(self.fc1(x))
        gb = self.film(sex.view(-1, 1))
        h = h * (1 + torch.tanh(gb[:, :128])) + torch.tanh(gb[:, 128:])
        h = h + self.block(h)
        return self.head(h).squeeze(-1), self.band(h)
def setw(lin, Wm, b):
    lin.weight.data = torch.tensor(Wm, dtype=torch.float32)
    lin.bias.data = torch.tensor(b, dtype=torch.float32)

nets = []
for sd in W["seeds"]:
    net = PriorsMLP(len(W["genes"]))
    setw(net.fc1, sd["fc1W"], sd["fc1B"])
    setw(net.block[0], sd["blk1W"], sd["blk1B"])
    setw(net.block[2], sd["blk2W"], sd["blk2B"])
    setw(net.head, sd["headW"], sd["headB"])
    setw(net.band, sd["bandW"], sd["bandB"])
    setw(net.film[0], sd["film1W"], sd["film1B"])
    setw(net.film[2], sd["film2W"], sd["film2B"])
    net.eval()
    nets.append(net)

# numpy forward with the same rounded params (mirrors JS exactly)
def np_forward(x, sexv, sd):
    h = np.maximum(sd["fc1W"] @ x + np.array(sd["fc1B"]), 0.0)
    f1 = np.array(sd["film1B"]) + sexv * np.array(sd["film1W"]).ravel()
    gb = np.tanh(np.array(sd["film2W"]) @ f1 + np.array(sd["film2B"]))
    h = h * (1 + gb[:128]) + gb[128:]
    t1 = np.maximum(np.array(sd["blk1W"]) @ h + np.array(sd["blk1B"]), 0.0)
    h = h + (np.array(sd["blk2W"]) @ t1 + np.array(sd["blk2B"]))
    corr = float(np.array(sd["headW"]).ravel() @ h + sd["headB"][0])
    bl = np.array(sd["bandW"]) @ h + np.array(sd["bandB"])
    return h, corr, bl

i = 0
x_raw = X[i]
x_std = (x_raw - np.array(W["geneMean"])) / np.array(W["geneSd"])
ap = W["anchor"]["intercept"] + float(np.array(W["anchor"]["coef"]) @ x_std)
tanh_mean = np.mean([np.tanh(np_forward(x_std, sex[i], sd)[1]) for sd in W["seeds"]])
deep = ap + W["residSd"] * tanh_mean

layers = {"x_std": x_std[:8].tolist(), "anchor": ap}
with torch.no_grad():
    xt = torch.tensor(x_std, dtype=torch.float32).view(1, -1)
    stt = torch.tensor([float(sex[i])])
    for k, net in enumerate(nets):
        h = torch.relu(net.fc1(xt))
        layers[f"h_fc1_{k}"] = h.numpy().ravel().tolist()
        gb = net.film(stt.view(-1, 1)).view(-1)
        h = h * (1 + torch.tanh(gb[:128])) + torch.tanh(gb[128:])
        layers[f"h_film_{k}"] = h.numpy().ravel().tolist()
        t1_t = torch.relu(net.block[0](h))
        layers[f"t1_{k}"] = t1_t.numpy().ravel().tolist()
        h_res_t = (h + net.block(h)).numpy().ravel()
        layers[f"h_res_{k}"] = h_res_t.tolist()
        sd = W["seeds"][k]
        t1_np = np.maximum(np.array(sd["blk1W"]) @ h.numpy().ravel() + np.array(sd["blk1B"]), 0.0)
        h_res_np = h.numpy().ravel() + (np.array(sd["blk2W"]) @ t1_np + np.array(sd["blk2B"]))
        layers[f"h_res_np_{k}"] = h_res_np.tolist()
        print(f"  seed{k}: torch-vs-numpy h_res maxdiff {np.abs(h_res_t - h_res_np).max():.2e}")
        w2 = torch.tensor(np.array(sd["blk2W"]), dtype=torch.float32)
        b2 = torch.tensor(np.array(sd["blk2B"]), dtype=torch.float32)
        print(f"    blk2.weight vs JSON maxdiff {(net.block[2].weight.data - w2).abs().max():.2e}")
        print(f"    block[2](t1_t) vs JSON-blk2(t1_t) maxdiff {(net.block[2](t1_t).view(-1) - (w2 @ t1_t.view(-1) + b2)).abs().max():.2e}")
        print(f"    t1_t vs t1_np maxdiff {(t1_t.view(-1) - torch.tensor(t1_np)).abs().max():.2e}")
        print(f"    h(film) vs h_film dump maxdiff {(h.view(-1) - torch.tensor(layers[f'h_film_{k}'])).abs().max():.2e}")
        print(f"    block structure: {net.block}")
        h_manual = h + net.block[2](torch.relu(net.block[0](h)))
        print(f"    h+block(h) vs h+block[2](relu(block[0](h))) maxdiff {(h_res_t - h_manual.view(-1)).abs().max():.2e}")
        print(f"    h_res dump vs h+block(h) maxdiff {(torch.tensor(h_res_t) - (h + net.block(h)).view(-1)).abs().max():.2e}")
        c_ = net.head(h)
        b_ = net.band(h)
        layers[f"corr_{k}"] = float(c_)
        layers[f"bandL_{k}"] = float(b_.numpy().ravel()[2])

# torch full-precision reference (nets rebuilt from rounded weights, same as JS)
with torch.no_grad():
    xt = torch.tensor(x_std, dtype=torch.float32).view(1, -1)
    st = torch.tensor([float(sex[i])])
    out = []
    for net in nets:
        c_, b_ = net(xt, st)
        out.append((float(c_), b_.numpy().ravel()))
layers["corr_torch"] = [c for c, _ in out]
layers["bandL_torch"] = [float(b[2]) for _, b in out]
json.dump(layers, open(os.path.join(ROOT, "data", "weights_diag.json"), "w"), indent=1)
print("diag written; sample0 anchor:", ap, " deep:", deep)
