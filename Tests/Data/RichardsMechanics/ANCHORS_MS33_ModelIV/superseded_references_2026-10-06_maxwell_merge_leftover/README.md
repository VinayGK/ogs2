# Superseded reference frames moved here by the maxwell conjugate merge, 2026-10-06

Files: ms33_modelIV_pellets_ts_16494_t_17280000.000000.vtu 

These are the IC-fix V0 reference frames (registered 2026-09-30 on gen5_icfix_conformant_2026-09-30,
d64c47cb30). On the submitted line they were superseded on 2026-10-02 (Mass-conserving A) and
2026-10-05 (K-fix set, E = 52 MPa); identical copies are in
`../superseded_references_2026-10-02_masscons_1a/`.

The merge of the submitted line into dsm_native_maxwell_conjugate (integration branch
integration/maxwell_conjugate_2026-10-06) left them beside the PRJ, because the maxwell line had
added them there. Beside the PRJ they match the vtkdiff regex of the deck's test_definition next
to the live reference, and the ctest fails (MEASURED 2026-10-06: Reference, III, IV and both VII
decks failed with "Error opening file ...ts_<old>..."; Model I passed). They are moved, not
deleted. HISTORICAL; not referenced by any test.
