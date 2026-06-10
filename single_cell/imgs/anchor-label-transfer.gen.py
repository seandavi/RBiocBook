#!/usr/bin/env python3
"""Original 'perspective-plane point cloud' schematic of label transfer, in two
steps: (1) anchors carry labels from the scRNA reference onto the gray scATAC
query, recoloring it (hatched); (2) the two assays merge into one plane.
RNA cells = solid fill; ATAC cells = hatched fill. No title/caption (done in text)."""
import random

cells = []
def add(s): cells.append(s)

def ellipse(cid, x, y, style, d=16):
    add(f'<mxCell id="{cid}" style="{style}" vertex="1" parent="1">'
        f'<mxGeometry x="{x:.0f}" y="{y:.0f}" width="{d}" height="{d}" as="geometry"/></mxCell>')

def edge(cid, src, tgt, style):
    add(f'<mxCell id="{cid}" style="{style}" edge="1" parent="1" source="{src}" target="{tgt}">'
        f'<mxGeometry relative="1" as="geometry"/></mxCell>')

def line(cid, x1, y1, x2, y2, style):
    add(f'<mxCell id="{cid}" style="{style}" edge="1" parent="1"><mxGeometry relative="1" as="geometry">'
        f'<mxPoint x="{x1:.0f}" y="{y1:.0f}" as="sourcePoint"/><mxPoint x="{x2:.0f}" y="{y2:.0f}" as="targetPoint"/>'
        f'</mxGeometry></mxCell>')

def text(cid, x, y, w, h, val, fs=12, bold=False, italic=False, color="#000000"):
    fst = (1 if bold else 0) + (2 if italic else 0)
    add(f'<mxCell id="{cid}" value="{val}" style="text;html=1;align=center;fontSize={fs};fontStyle={fst};fontColor={color};" '
        f'vertex="1" parent="1"><mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>')

types = {  # cell type -> (fill, stroke)
    "T":    ("#36c5d6", "#1b8a96"),
    "B":    ("#e85a9b", "#b13571"),
    "Mono": ("#f08a3c", "#c2651f"),
    "NK":   ("#6dc25a", "#4a9139"),
}
GRAY = ("#d9d9d9", "#9a9a9a")
W, H, SKEW = 240, 120, 40

def place(ox, oy, lx, ly):
    return ox + lx + SKEW * (H - ly) / H, oy + ly

def gen_layout(centers, npts, seed):
    rnd = random.Random(seed)
    out = []
    for t, (cx, cy) in centers.items():
        for i in range(npts[t]):
            out.append((t, i, cx + rnd.uniform(-20, 20), cy + rnd.uniform(-18, 18)))
    return out

def plane_border(prefix, ox, oy):
    tl = place(ox, oy, 0, 0); tr = place(ox, oy, W, 0)
    br = place(ox, oy, W, H); bl = place(ox, oy, 0, H)
    st = "endArrow=none;html=1;strokeColor=#555555;strokeWidth=1.5;"
    line(prefix+"_t", *tl, *tr, st); line(prefix+"_r", *tr, *br, st)
    line(prefix+"_b", *br, *bl, st); line(prefix+"_l", *bl, *tl, st)

def style_for(mode, t):
    fill, stroke = types[t]
    if mode == "gray":
        fill, stroke = GRAY
        return f"ellipse;html=1;fillColor={fill};strokeColor={stroke};"
    if mode == "hashed":
        return f"ellipse;html=1;sketch=1;fillStyle=hachure;fillColor={fill};strokeColor={stroke};"
    return f"ellipse;html=1;fillColor={fill};strokeColor={stroke};"

def draw_plane(prefix, ox, oy, layout, mode):
    plane_border(prefix, ox, oy)
    ids = {}
    for (t, i, lx, ly) in layout:
        x, y = place(ox, oy, lx, ly)
        cid = f"{prefix}_{t}_{i}"
        ellipse(cid, x - 8, y - 8, style_for(mode, t))
        ids[(t, i)] = cid
    return ids

# RNA and ATAC are organized differently (different centers, counts, jitter)
rna_centers  = {"T": (48, 40), "B": (150, 30), "Mono": (105, 92), "NK": (200, 58)}
atac_centers = {"T": (62, 54), "B": (136, 44), "Mono": (122, 80), "NK": (186, 72)}
mrg_centers  = {"T": (54, 46), "B": (150, 36), "Mono": (112, 88), "NK": (197, 62)}
npts_rna  = {"T": 7, "B": 5, "Mono": 6, "NK": 8}
npts_atac = {"T": 6, "B": 6, "Mono": 7, "NK": 7}

lay_rna  = gen_layout(rna_centers, npts_rna, 1)
lay_atac = gen_layout(atac_centers, npts_atac, 2)
lay_mR   = gen_layout(mrg_centers, npts_rna, 3)
lay_mA   = gen_layout(mrg_centers, npts_atac, 4)

# --- Scene 1: RNA (solid) over ATAC (gray), with anchors ---
s1r = draw_plane("s1r", 40, 60, lay_rna, "solid")
s1q = draw_plane("s1q", 40, 250, lay_atac, "gray")
ast = "endArrow=none;html=1;dashed=1;strokeColor=#888888;"
for t in types:
    for i in range(min(2, npts_rna[t], npts_atac[t])):
        edge(f"an_{t}_{i}", s1r[(t, i)], s1q[(t, i)], ast)

# --- Scene 2: RNA (solid) over ATAC (now hatched-colored) ---
draw_plane("s2r", 470, 60, lay_rna, "solid")
draw_plane("s2q", 470, 250, lay_atac, "hashed")

# --- Scene 3: merged single plane (RNA solid + ATAC hatched together) ---
draw_plane("s3R", 900, 150, lay_mR, "solid")
draw_plane("s3A", 900, 150, lay_mA, "hashed")

# --- step arrows ---
astyle = "endArrow=classic;html=1;strokeColor=#cc2a82;strokeWidth=3;fillColor=#cc2a82;"
line("arr1", 360, 235, 452, 235, astyle)
line("arr2", 800, 235, 884, 235, astyle)

# --- minimal labels (no title / caption) ---
text("l1a", 60, 38, 240, 18, "scRNA-seq reference", 12, bold=True)
text("l1b", 60, 372, 240, 18, "scATAC-seq query", 12, bold=True)
text("lanch", 150, 205, 80, 30, "anchors", 11, italic=True)
text("larr1", 330, 210, 150, 18, "1. transfer labels", 11, italic=True)
text("l2b", 490, 372, 240, 18, "scATAC-seq query (labeled)", 12, bold=True)
text("larr2", 790, 210, 110, 18, "2. merge", 11, italic=True)
text("l3", 905, 128, 240, 18, "integrated", 12, bold=True)

doc = ('<mxfile host="rbiocbook"><diagram name="anchors" id="anchors">'
       '<mxGraphModel dx="1200" dy="640" grid="0" gridSize="10" guides="1" tooltips="1" connect="1" arrows="1" '
       'fold="1" page="1" pageScale="1" pageWidth="1200" pageHeight="410" math="0" shadow="0"><root>'
       '<mxCell id="0"/><mxCell id="1" parent="0"/>' + "".join(cells) +
       '</root></mxGraphModel></diagram></mxfile>')
import sys
sys.stdout.write(doc)
