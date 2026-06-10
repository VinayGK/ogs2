#!/usr/bin/env python3
"""
Split the axisymmetric column boundary into named face meshes that OGS RM
boundary conditions can reference. Lab oedometer dimensions:
  radius  x in [0, 0.0175] m   (35 mm diameter)
  height  y in [0, 0.0125] m   (12.5 mm)
Faces: left (x=0, axis), right (x=R, ring wall), bottom (y=0, hydration),
       top (y=H, piston).
Writes col_*_left/right/bottom/top.vtu next to the domain mesh.

Pure-meshio extraction of edge cells (line elements) lying on each face, with
bulk_node_ids preserved so OGS can map BC nodes back to the domain.
"""
import numpy as np, meshio, os, sys

D = os.path.dirname(os.path.abspath(__file__))
dom = os.path.join(D, "col_r175_h125.vtu")
R, H = 0.0175, 0.0125
tol = 1e-9

m = meshio.read(dom)
pts = m.points
# OGS needs a bulk_node_ids array on the domain for BC mapping
nnodes = len(pts)
bulk_ids = np.arange(nnodes, dtype=np.uint64)
# write domain with bulk_node_ids (idempotent; adds the array OGS expects)
m.point_data["bulk_node_ids"] = bulk_ids
# also bulk_element_ids on cells
ncells = sum(len(cb.data) for cb in m.cells)
# find quad block
quad = None
for cb in m.cells:
    if cb.type == "quad":
        quad = cb.data
m.cell_data["bulk_element_ids"] = [np.arange(len(quad), dtype=np.uint64)]
meshio.write(dom, m, file_format="vtu", binary=False)

def face(mask_fn, name):
    # collect nodes on the face
    on = np.where(mask_fn(pts))[0]
    onset = set(on.tolist())
    # build line cells from quad edges fully on the face
    edges = []
    for q in quad:
        # quad node order: 0-1-2-3 around; edges (0,1)(1,2)(2,3)(3,0)
        for a, b in ((0,1),(1,2),(2,3),(3,0)):
            na, nb = int(q[a]), int(q[b])
            if na in onset and nb in onset:
                edges.append((na, nb))
    if not edges:
        print(f"  WARN: no edges for face {name}")
        return
    used = sorted({n for e in edges for n in e})
    remap = {old:i for i,old in enumerate(used)}
    fpts = pts[used]
    lines = np.array([[remap[a], remap[b]] for a,b in edges], dtype=int)
    fm = meshio.Mesh(
        fpts,
        [("line", lines)],
        point_data={"bulk_node_ids": np.array(used, dtype=np.uint64)},
    )
    out = os.path.join(D, f"col_r175_h125_{name}.vtu")
    meshio.write(out, fm, file_format="vtu", binary=False)
    print(f"  {name}: {len(used)} nodes, {len(lines)} line cells -> {os.path.basename(out)}")

print("splitting column boundary (lab dims R=17.5mm, H=12.5mm):")
face(lambda p: np.abs(p[:,0]-0.0) < tol, "left")    # axis
face(lambda p: np.abs(p[:,0]-R)   < tol, "right")   # ring wall
face(lambda p: np.abs(p[:,1]-0.0) < tol, "bottom")  # hydration face
face(lambda p: np.abs(p[:,1]-H)   < tol, "top")     # piston
print("done.")
