#!/usr/bin/env python3
"""
MGR23 two-material axisymmetric column mesh (FEBEX, verified geometry D4.1:3407-3411).
Cell: inner diameter 10.0 cm -> axisym radius r in [0, 0.05] m; sample length 10 cm ->
height y in [0, 0.10] m. Pellet layer 0-5 cm (BOTTOM, hydration face), block 5-10 cm (top).
MaterialIDs: 0 = pellet (bottom), 1 = block (top).
Hydration through the bottom (pellet end) per D4.1:3419 'water intake through the bottom
surface'; assembly was overturned so pellets sit at the wetting face (D4.1:3417).

Mesh: structured quad, nx=10 (radial), ny=40 (axial, 2.5 mm cells -> resolves the front
and the 5cm material interface lands on an element boundary). Boundary face meshes
(left/right/top/bottom) written with bulk_node_ids for OGS BCs.
"""
import numpy as np, meshio, os

D = os.path.dirname(os.path.abspath(__file__))
R, H = 0.05, 0.10
nx, ny = 10, 40
xs = np.linspace(0, R, nx + 1)
ys = np.linspace(0, H, ny + 1)

pts = np.array([[x, y, 0.0] for y in ys for x in xs], dtype=float)
def nid(i, j): return j * (nx + 1) + i   # i radial, j axial

quads, matids = [], []
for j in range(ny):
    yc = 0.5 * (ys[j] + ys[j + 1])
    mat = 0 if yc < 0.05 else 1          # 0 pellet (bottom <5cm), 1 block (top)
    for i in range(nx):
        quads.append([nid(i, j), nid(i + 1, j), nid(i + 1, j + 1), nid(i, j + 1)])
        matids.append(mat)
quads = np.array(quads); matids = np.array(matids, dtype=np.int32)

dom = meshio.Mesh(
    pts, [("quad", quads)],
    cell_data={"MaterialIDs": [matids],
               "bulk_element_ids": [np.arange(len(quads), dtype=np.uint64)]},
    point_data={"bulk_node_ids": np.arange(len(pts), dtype=np.uint64)},
)
meshio.write(f"{D}/mgr_col.vtu", dom, file_format="vtu", binary=False)
print(f"domain: {len(pts)} nodes, {len(quads)} quads, "
      f"{int((matids==0).sum())} pellet + {int((matids==1).sum())} block cells")

tol = 1e-9
def face(mask_fn, name):
    on = set(np.where(mask_fn(pts))[0].tolist())
    edges = []
    for q in quads:
        for a, b in ((0,1),(1,2),(2,3),(3,0)):
            if int(q[a]) in on and int(q[b]) in on:
                edges.append((int(q[a]), int(q[b])))
    used = sorted({n for e in edges for n in e})
    rm = {o:i for i,o in enumerate(used)}
    lines = np.array([[rm[a], rm[b]] for a,b in edges])
    fm = meshio.Mesh(pts[used], [("line", lines)],
                     point_data={"bulk_node_ids": np.array(used, dtype=np.uint64)})
    meshio.write(f"{D}/mgr_col_{name}.vtu", fm, file_format="vtu", binary=False)
    print(f"  {name}: {len(used)} nodes")

face(lambda p: np.abs(p[:,0]-0.0) < tol, "left")    # axis
face(lambda p: np.abs(p[:,0]-R)   < tol, "right")   # ring wall
face(lambda p: np.abs(p[:,1]-0.0) < tol, "bottom")  # hydration (pellet end)
face(lambda p: np.abs(p[:,1]-H)   < tol, "top")     # piston (block end)
print("done.")
