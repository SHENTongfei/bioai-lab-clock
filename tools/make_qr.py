# -*- coding: utf-8 -*-
"""make_qr.py — tech-styled QR for the clock site, verified decodable, saved to desktop."""
import qrcode, math, os
import numpy as np
from PIL import Image, ImageDraw, ImageFont, ImageFilter

URL = "https://shentongfei.github.io/bioai-lab-clock/"
OUT = r"C:\Users\TS\Desktop\LoCRANE-clock-qr.png"
LOGO = r"C:\Users\TS\Desktop\bioai-lab-clock\assets\logo-256.png"

# ---- matrix (H error correction, quiet zone handled manually) ----
qr = qrcode.QRCode(error_correction=qrcode.constants.ERROR_CORRECT_H, box_size=1, border=0)
qr.add_data(URL)
qr.make(fit=True)
M = np.array(qr.get_matrix(), dtype=bool)   # True = dark module
N = M.shape[0]
print("matrix", N, "x", N)

# ---- canvas ----
SS = 2                       # supersample
CELL = 24 * SS               # px per module
QUIET = 4 * CELL             # quiet zone
TOP = 150 * SS               # header band
BOTTOM = 190 * SS            # footer band
W = N * CELL + 2 * QUIET
H = TOP + N * CELL + 2 * QUIET + BOTTOM
BG = (247, 246, 252)
EDGE = (10, 14, 28)
img = Image.new("RGB", (W, H), BG)
d = ImageDraw.Draw(img)

# soft frame + glow
pad = 26 * SS
d.rounded_rectangle([pad, pad, W - pad, H - pad], radius=40 * SS,
                    fill=(16, 22, 48), outline=(120, 110, 170), width=3 * SS)
inner = Image.new("RGB", (W, H), (16, 22, 48))
# subtle corner glows on the dark frame
glow = Image.new("L", (W, H), 0)
gd = ImageDraw.Draw(glow)
for (cx, cy, r, a) in [(90*SS, 90*SS, 260*SS, 60), (W-90*SS, H-200*SS, 300*SS, 55)]:
    gd.ellipse([cx-r, cy-r, cx+r, cy+r], fill=a)
glow = glow.filter(ImageFilter.GaussianBlur(60 * SS))
frame_mask = Image.new("L", (W, H), 0)
fd = ImageDraw.Draw(frame_mask)
fd.rounded_rectangle([pad, pad, W - pad, H - pad], radius=40 * SS, fill=255)
frame = Image.composite(inner, img, frame_mask)
frame_col = np.array(frame, dtype=np.int16)
glow_col = np.array(Image.new("RGB", (W, H), (110, 90, 160)), dtype=np.int16)
alpha = np.array(glow, dtype=np.int16)[..., None] / 255.0
mixed = (frame_col * (1 - alpha * 0.55) + glow_col * alpha * 0.55).astype(np.uint8)
img.paste(Image.fromarray(mixed), (0, 0))
d = ImageDraw.Draw(img)

# header text (on dark frame)
def font(sz, bold=True):
    for name in (["msyhbd.ttc"] if bold else ["msyh.ttc"]):
        try: return ImageFont.truetype("C:/Windows/Fonts/" + name, sz)
        except Exception: pass
    return ImageFont.load_default()
f_head = font(56 * SS // 2)
d.text((W // 2, pad + 66 * SS), "BioAI Lab · LoCRANE", font=f_head, fill=(238, 241, 251), anchor="mm")
f_sub = font(30 * SS // 2, bold=False)
d.text((W // 2, pad + 118 * SS), "扫一扫 · 打开老化时钟", font=f_sub, fill=(163, 173, 210), anchor="mm")

# ---- code area geometry ----
ox, oy = QUIET + (W - 2 * QUIET - N * CELL) // 2, TOP + QUIET
c2 = CELL / 2

def grad_color(u, v):   # u,v in 0..1 → deep indigo → plum
    t = (u * 0.55 + v * 0.45)
    c1 = (61, 63, 107); c2p = (138, 79, 134)
    return tuple(int(c1[i] + (c2p[i] - c1[i]) * t) for i in range(3))

def rounded_ring(dr, x0, y0, x1, y1, r, width, fill):
    dr.rounded_rectangle([x0, y0, x1, y1], radius=r, outline=fill, width=width)

# dark inner panel for the code (keeps scanner polarity: dark dots on light panel)
panel0 = (ox - QUIET // 2, oy - QUIET // 2)
panel1 = (ox + N * CELL + QUIET // 2, oy + N * CELL + QUIET // 2)
d.rounded_rectangle([panel0[0] - 8*SS, panel0[1] - 8*SS, panel1[0] + 8*SS, panel1[1] + 8*SS],
                    radius=24 * SS, fill=(250, 249, 253))

finder = [(0, 0), (N - 7, 0), (0, N - 7)]
fcol = (45, 47, 84)   # deep indigo, verified scannable
dot_r = CELL * 0.36
for r in range(N):
    for c in range(N):
        if not M[r, c]:
            continue
        in_finder = (r < 7 and c < 7) or (r < 7 and c >= N - 7) or (r >= N - 7 and c < 7)
        if in_finder:
            continue
        cx = ox + c * CELL + c2
        cy = oy + r * CELL + c2
        col = grad_color(c / N, r / N)
        d.ellipse([cx - dot_r, cy - dot_r, cx + dot_r, cy + dot_r], fill=col)

for (fr, fc) in finder:
    x0, y0 = ox + fc * CELL, oy + fr * CELL
    x1, y1 = x0 + 7 * CELL, y0 + 7 * CELL
    d.rounded_rectangle([x0, y0, x1, y1], radius=int(0.6 * CELL), fill=fcol)
    d.rounded_rectangle([x0 + CELL, y0 + CELL, x1 - CELL, y1 - CELL], radius=int(0.3 * CELL), fill=(250, 249, 253))
    d.rounded_rectangle([x0 + 2 * CELL, y0 + 2 * CELL, x1 - 2 * CELL, y1 - 2 * CELL], radius=int(0.5 * CELL), fill=fcol)

# center logo on white rounded card (covers ~9 modules, H-level EC covers it)
L = 9 * CELL
lcx, lcy = ox + N * CELL / 2, oy + N * CELL / 2
d.rounded_rectangle([lcx - L/2 - 6*SS, lcy - L/2 - 6*SS, lcx + L/2 + 6*SS, lcy + L/2 + 6*SS],
                    radius=18 * SS, fill=(255, 255, 255), outline=(210, 200, 235), width=2 * SS)
logo = Image.open(LOGO).convert("RGBA").resize((int(L * 0.92), int(L * 0.92)), Image.LANCZOS)
img.paste(logo, (int(lcx - logo.width / 2), int(lcy - logo.height / 2)), logo)

# footer (on dark frame)
f_foot = font(40 * SS // 2)
d.text((W // 2, H - BOTTOM + 62 * SS), "LoCRANE 老化时钟 · Aging Clock", font=f_foot, fill=(238, 241, 251), anchor="mm")
f_url = font(26 * SS // 2, bold=False)
d.text((W // 2, H - BOTTOM + 116 * SS), "shentongfei.github.io/bioai-lab-clock", font=f_url, fill=(150, 160, 200), anchor="mm")

# ---- downscale (AA) and save ----
final = img.resize((W // SS, H // SS), Image.LANCZOS)
final.save(OUT, "PNG")
print("saved", OUT, final.size)

# ---- decode verification (crop to code panel, invert-free standard polarity) ----
import cv2
a = np.array(final)
det = cv2.QRCodeDetector()
data, _, _ = det.detectAndDecode(a)
if not data:
    # retry on the code region only
    scale = final.width / W
    box = (int((ox - QUIET) * scale), int((oy - QUIET + TOP * 0) * scale),
           int((ox + N * CELL + QUIET) * scale), int((oy + N * CELL + QUIET) * scale))
    data, _, _ = det.detectAndDecode(a[box[1]:box[3], box[0]:box[2]])
print("DECODED:", repr(data))
assert data == URL, "decode mismatch!"
print("VERIFIED: scans back to", URL)
