# SPDX-License-Identifier: BSD-3-Clause
# MS33 Model III - technological-gap closure by a Python BC SWITCH on the
# bentonite outer radial boundary (free swelling -> rigid Dirichlet at gap).
#
# Requested by Vinay 2026-06-23: "if contact is not possible, run it with
# python BC switching from free swelling to dirichlet once gap is closed."
#
# Mechanism (mirrors the PROVEN Task-13 2023 getDirichletBCValue switch,
# ~/ogs-models/EBS/Task13/2023_August_submitted/prj-common/top-gap-bc.py,
# which actually ran and produced the submitted gap-closure figures - applied
# here to the RADIAL boundary instead of the axial top):
#   the bentonite cylinder (r = 23 mm since 2026-09-08, spec Fig3gapelementa.jpg
#   "46 mm" clay in the "50 mm" cell; r = 25 mm through gen 4, which put the
#   wall at r = 27 mm) swells radially into a 2 mm technological gap. The
#   switch tests the DISPLACEMENT u_r, not the radius, so the r = 23 mm mesh
#   latches at r = 25 mm (= the spec cell wall) with NO change to this script. While the outer-boundary radial displacement
#   u_r < GAP the boundary is LEFT FREE (return apply_bc=False) -> free
#   swelling, gap still open. Once u_r reaches GAP the node is CONSTRAINED to
#   u_r = GAP (return apply_bc=True, value=GAP) -> rigid container contact;
#   further swelling then builds confined swelling pressure against the wall.
#
# Once latched (u_r = GAP exactly), the test `u_r < GAP` is False so the node
# stays constrained - a stable one-way latch (swelling here is monotonic, no
# unloading), so no contact chatter. GAP = 2 mm is the EURAD-2 MS33 Model III
# spec geometry (PRJ header), NOT a free parameter.
try:
    import ogs.callbacks as OpenGeoSys
except ModuleNotFoundError:
    import OpenGeoSys

GAP = 0.002  # technological gap width [m] (2 mm), EURAD-2 MS33 Model III spec


class RadialGapSwitch(OpenGeoSys.BoundaryCondition):
    def getDirichletBCValue(self, t, coords, node_id, primary_vars):
        # RM primary_vars order = [pressure, u_x(=u_r radial), u_y(=u_z axial)].
        u_r = primary_vars[1]
        if u_r < GAP:
            return (False, 0.0)   # gap open -> free-swelling boundary
        return (True, GAP)        # gap closed -> rigid Dirichlet wall at +2 mm


gap_switch = RadialGapSwitch()
