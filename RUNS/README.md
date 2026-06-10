# RUNS/ — per-run snapshot bundles (the standard mechanism)

Standardized 2026-06-10 at Vinay's instruction; rule codified in CLAUDE.md §13.

Every completed simulation campaign gets ONE self-contained snapshot folder
here, committed and pushed, so that any result can be traced to the exact
code, input, and post-processing that produced it.

## Naming

```
RUNS/<YYYY-MM-DD>T<HHMM>_<branch>_<short-slug>/
```

Date+time = when the run (last sim of the campaign) finished, local time.
Branch = the git branch of the binary that produced the results.

## Mandatory contents

- `README.md` — starts with the RUN BUNDLE header block:
  bundle name, finish time, **branch @ commit**, **binary path**, where the
  full raw outputs live (untracked, local), and pointers to report / beamer /
  tex artifacts (or `TODO` if not yet written). Then: what was run, why, the
  headline results, findings, caveats. The README must allow a cold reader to
  re-run and re-derive every figure.
- The PRJ file(s) **as run** (with any variant edits applied).
- All post-processing / comparison scripts.
- All figures (PNG/PDF) and extracted tables (CSV) presented anywhere.
- Run logs (or their tails if huge).
- `final_state/` — the final-time VTU per variant (NOT the full series).
- Report / beamer / paper-passage PDFs+tex if they exist for this run; else a
  `TODO` line in the README header.

## Size policy

Full VTU time series stay OUT of git (local, path recorded in the README).
Meshes, final states, figures, CSVs, logs, PRJs go in. Target < ~30 MB per
bundle; justify exceptions in the README.

## Provenance rules (inherited)

Bundles obey CLAUDE.md §5 (predicted ≠ verified in every README claim) and
§12 (DSM PRJs carry their provenance header). A bundle is a HISTORICAL
RECORD: never edit results retroactively — supersede with a new bundle and
cross-link.
