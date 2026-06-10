#!/usr/bin/env python3
"""Reusable provenance footer for MS33 inter-team overlay figures.

A single small grey line stamped at the bottom of every overlay so that any
exported figure is traceable back to the exact OGS branch + commit, the solver
binary, the simulation run/output dir its data came from, and the plot date.

Use it in any MS33 overlay plotting script: import the function, then call it on
the Figure object just before each savefig, passing the run id for that figure.

    from ms33_provenance_footer import add_provenance_footer
    ...
    add_provenance_footer(fig, run="conv3_out/final (conv3_diff_final)")
    fig.savefig(path, dpi=220)

By default `commit` auto-detects the current HEAD of the
dsm_native_pdisj_maxwell worktree, so figures produced from a *new* run stamp
whatever HEAD is checked out at plot time. Pass commit="..." explicitly only to
stamp figures that were produced on a known earlier binary (e.g. the existing
overlays were run on 36dfeaddcc, the pre-Option-B binary).

Footer text:
  OGS <branch> @<commit> · binary pdisj_maxwell_revref_20260605 · run: <run> · <date>
"""
import subprocess

# OGS worktree whose HEAD identifies the binary that produced a fresh run.
_WORKTREE = "/Users/vinaykumar/git/ogs-worktrees/dsm_native_pdisj_maxwell_wt"
# Fixed solver-binary tag (the compiled binary in use for the MS33 suite).
_BINARY = "pdisj_maxwell_revref_20260605"


def detect_commit(worktree=_WORKTREE):
    """Return the short (10-char) HEAD commit of the OGS worktree, or 'unknown'
    if git is unavailable / the path is not a repo. Never raises."""
    try:
        out = subprocess.run(
            ["git", "-C", worktree, "rev-parse", "--short=10", "HEAD"],
            capture_output=True, text=True, check=True,
        )
        return out.stdout.strip() or "unknown"
    except Exception:
        return "unknown"


def provenance_text(run, commit=None, branch="dsm_native_pdisj_maxwell",
                    date="2026-06-08"):
    """Build the single-line footer string (no rendering). commit=None -> HEAD."""
    if commit is None:
        commit = detect_commit()
    return (f"OGS {branch} @{commit} · binary {_BINARY} "
            f"· run: {run} · {date}")


def add_provenance_footer(fig, run, commit=None, branch="dsm_native_pdisj_maxwell",
                          date="2026-06-08", bottom=0.085, fontsize=7):
    """Stamp a small grey provenance line at the bottom of `fig`.

    fig    : matplotlib Figure to annotate (call once per figure, before savefig).
    run    : simulation output dir / id that produced this figure's data.
    commit : OGS short commit to stamp; None -> auto-detect current worktree HEAD.
    branch : OGS branch name (default the DSM Maxwell branch).
    date   : plot date stamp.
    bottom : fraction of the figure height reserved below the axes for the line,
             so it never overlaps the plotted content.

    Returns the footer text actually stamped (handy for logging/reporting).
    """
    text = provenance_text(run, commit=commit, branch=branch, date=date)
    # Reserve room beneath the axes so the footer never collides with them.
    # subplots_adjust must come before adding the text; guard against scripts
    # that already pulled the bottom in further than our default.
    try:
        cur = fig.subplotpars.bottom
    except Exception:
        cur = 0.0
    fig.subplots_adjust(bottom=max(cur, bottom))
    fig.text(0.5, 0.005, text, ha="center", va="bottom", fontsize=fontsize,
             color="0.4")
    return text


if __name__ == "__main__":
    # Self-test: show the auto-detected line and an explicitly-stamped one.
    print("auto-HEAD :", provenance_text("demo_run"))
    print("explicit  :", provenance_text("demo_run", commit="36dfeaddcc"))
