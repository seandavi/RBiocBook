#!/usr/bin/env python3
"""Generate an original 'perspective-plane point cloud' schematic of anchor-based
label transfer (scRNA reference -> scATAC query) as a drawio file."""
import random
random.seed(7)

cells = []
def add(s): cells.append(s)

def ellipse(cid, x, y, fill, stroke, d=18):
    add(f'<mxCell id="{cid}" style="ellipse;html=1;fillColor={fill};strokeColor={stroke};" '
        f'vertex="1" parent="1"><mxGeometry x="{x:.0f}" y="{y:.0f}" width="{d}" height="{d}" as="geometry"/></mxCell>')

def edge(cid, src, tgt, style):
    add(f'<mxCell id="{cid}" style="{style}" edge="1" parent="1" source="{src}" target="{tgt}">'
        f'<mxGeometry relative="1" as="geometry"/></mxCell>')

def line(cid, x1, y1, x2, y2, style):
    add(f'<mxCell id="{cid}" style="{style}" edge="1" parent="1"><mxGeometry relative="1" as="geometry">'
        f'<mxPoint x="{x1:.0f}" y="{y1:.0f}" as="sourcePoint"/><mxPoint x="{x2:.0f}" y="{y2:.0f}" as="targetPoint"/>'
        f'</mxGeometry></mxCell>')

def text(cid, x, y, w, h, val, fs=12, bold=False, italic=False, color="#000000", align="center"):
    fst = (1 if bold else 0) + (2 if italic else 0)
    add(f'<mxCell id="{cid}" value="{val}" style="text;html=1;align={align};fontSize={fs};fontStyle={fst};fontColor={color};" '
        f'vertex="1" parent="1"><mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>')

# cell types: name -> (fill, stroke), and a local cluster center within the plane content box
types = {
    "T":    ("#36c5d6", "#1b8a96"),
    "B":    ("#e85a9b", "#b13571"),
    "Mono": ("#f08a3c", "#c2651f"),
    "NK":   ("#6dc25a", "#4a9139"),
}
centers = {"T": (55, 48), "B": (165, 36), "Mono": (120, 112), "NK": (245, 72)}
npts = {"T": 7, "B": 5, "Mono": 7, "NK": 8}
GRAY = ("#d9d9d9", "#9a9a9a")

W, H, SKEW = 300, 150, 52

def place(ox, oy, lx, ly):
    """map local (lx,ly) onto the skewed plane -> absolute (x,y)."""
    return ox + lx + SKEW * (H - ly) / H, oy + ly

def plane_border(prefix, ox, oy):
    tl = place(ox, oy, 0, 0); tr = place(ox, oy, W, 0)
    br = place(ox, oy, W, H); bl = place(ox, oy, 0, H)
    st = "endArrow=none;html=1;strokeColor=#444444;strokeWidth=1.5;"
    line(prefix+"t", *tl, *tr, st); line(prefix+"r", *tr, *br, st)
    line(prefix+"b", *br, *bl, st); line(prefix+"l", *bl, *tl, st)

# generate a fixed cluster layout (local coords + which point) shared by all planes
layout = []  # (type, idx, lx, ly)
for t, (cx, cy) in centers.items():
    for i in range(npts[t]):
        lx = cx + random.uniform(-22, 22)
        ly = cy + random.uniform(-20, 20)
        layout.append((t, i, lx, ly))

def draw_plane(prefix, ox, oy, mode):
    """mode: 'ref' (colored), 'unlabeled' (gray), 'pred' (colored = transferred)."""
    plane_border(prefix, ox, oy)
    ids = {}
    for (t, i, lx, ly) in layout:
        x, y = place(ox, oy, lx, ly)
        cid = f"{prefix}_{t}_{i}"
        if mode == "unlabeled":
            fill, stroke = GRAY
        else:
            fill, stroke = types[t]
        ellipse(cid, x - 9, y - 9, fill, stroke)
        ids[(t, i)] = cid
    return ids

# ---- planes ----
REF_OX, REF_OY = 70, 70
QRY_OX, QRY_OY = 70, 320
RES_OX, RES_OY = 660, 195

ref_ids = draw_plane("ref", REF_OX, REF_OY, "ref")
qry_ids = draw_plane("qry", QRY_OX, QRY_OY, "unlabeled")
res_ids = draw_plane("res", RES_OX, RES_OY, "pred")

# ---- anchor lines: connect a couple of points per type between ref and query ----
anchor_style = "endArrow=none;html=1;dashed=1;strokeColor=#888888;"
k = 0
for t in types:
    for i in range(min(2, npts[t])):
        edge(f"anch_{t}_{i}", ref_ids[(t, i)], qry_ids[(t, i)], anchor_style)
        k += 1

# ---- arrow query -> result (the transfer) ----
add('<mxCell id="arrow" style="endArrow=classic;html=1;strokeColor=#cc2a82;strokeWidth=3;fillColor=#cc2a82;" '
    'edge="1" parent="1"><mxGeometry relative="1" as="geometry">'
    '<mxPoint x="470" y="300" as="sourcePoint"/><mxPoint x="640" y="300" as="targetPoint"/></mxGeometry></mxCell>')

# ---- labels ----
text("title", 120, 12, 760, 26, "Anchor-based label transfer: borrowing labels across datasets", 16, bold=True)
text("lref", 70, 44, 350, 20, "scRNA-seq reference (cell types known)", 12, bold=True)
text("lqry", 70, 296, 350, 20, "scATAC-seq query (cell types unknown)", 12, bold=True)
text("lres", 640, 168, 360, 20, "scATAC-seq query (labels transferred)", 12, bold=True)
text("lanch", 430, 220, 150, 48, "anchors\n(mutual nearest neighbors)", 11, italic=True)
text("foot", 70, 520, 940, 22,
     "Mutual-nearest-neighbor cells act as anchors; the query inherits the cell-type label of the reference cells it anchors to.",
     12, italic=True)

doc = (
'<mxfile host="rbiocbook"><diagram name="anchors" id="anchors">'
'<mxGraphModel dx="1100" dy="640" grid="0" gridSize="10" guides="1" tooltips="1" connect="1" arrows="1" '
'fold="1" page="1" pageScale="1" pageWidth="1040" pageHeight="560" math="0" shadow="0"><root>'
'<mxCell id="0"/><mxCell id="1" parent="0"/>'
+ "".join(cells) +
'</root></mxGraphModel></diagram></mxfile>'
)
import sys
sys.stdout.write(doc)
