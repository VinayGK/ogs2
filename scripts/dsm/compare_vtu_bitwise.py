#!/usr/bin/env python3
"""Bitwise field-data comparison of two VTU output directories.

Regression tool for PI_OF_NL_EV_IMPLEMENTATION.md T-1 (off-mode bit-for-bit),
mirroring the 2026-06-09 strained-film verification: compares point
coordinates and EVERY point/cell data array of every *.vtu pair byte-for-byte
(numpy array_equal on the raw arrays); the embedded OGS_VERSION string is
intentionally NOT compared (two binaries, one input).

Usage: compare_vtu_bitwise.py DIR_REF DIR_NEW
Exit 0 = all VTUs bitwise-identical field data; 1 = any mismatch.
"""
import sys
from pathlib import Path

import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy


def load(path):
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(str(path))
    r.Update()
    return r.GetOutput()


def arrays(data):
    return {
        data.GetArrayName(i): vtk_to_numpy(data.GetArray(i))
        for i in range(data.GetNumberOfArrays())
    }


def compare(a_path, b_path):
    a, b = load(a_path), load(b_path)
    bad = []
    pa = vtk_to_numpy(a.GetPoints().GetData())
    pb = vtk_to_numpy(b.GetPoints().GetData())
    if not np.array_equal(pa, pb):
        bad.append("points")
    for kind, da, db in (
        ("point", arrays(a.GetPointData()), arrays(b.GetPointData())),
        ("cell", arrays(a.GetCellData()), arrays(b.GetCellData())),
    ):
        if set(da) != set(db):
            bad.append(f"{kind}-array-names {sorted(set(da) ^ set(db))}")
            continue
        for name in sorted(da):
            if not np.array_equal(da[name], db[name]):
                d = np.max(np.abs(da[name] - db[name]))
                bad.append(f"{kind}:{name} (max|diff|={d:g})")
    return bad


def main():
    ref, new = Path(sys.argv[1]), Path(sys.argv[2])
    ref_vtus = sorted(ref.glob("*.vtu"))
    new_vtus = sorted(new.glob("*.vtu"))
    if [p.name for p in ref_vtus] != [p.name for p in new_vtus]:
        print("FILE LIST MISMATCH")
        print(" ref:", [p.name for p in ref_vtus])
        print(" new:", [p.name for p in new_vtus])
        return 1
    rc = 0
    for r, n in zip(ref_vtus, new_vtus):
        bad = compare(r, n)
        if bad:
            rc = 1
            print(f"DIFF  {r.name}: {bad}")
        else:
            print(f"OK    {r.name}: bitwise-identical field data")
    print("RESULT:", "PASS (bitwise)" if rc == 0 else "FAIL")
    return rc


if __name__ == "__main__":
    sys.exit(main())
