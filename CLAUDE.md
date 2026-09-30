# main-git

Julia ED + TTN and Python QuOCS code for FQH / fractional Chern insulators in synthetic dimensions.

| Project | Status (read first) | Technical notes | Code |
|---|---|---|---|
| Optimal control & time evolution | [optimal-control/CONTEXT.md](optimal-control/CONTEXT.md) | [optimal-control/NOTES.md](optimal-control/NOTES.md) | `optimal-control/`, `exact-diag/` |
| Dipolar phase transition | [synth-dims/CONTEXT.md](synth-dims/CONTEXT.md) | — | `synth-dims/` |

- Project facts live in these files, not in Claude's memory. When a session produces a result,
  decision or bug, update CONTEXT.md (state, results, open questions, log) and NOTES.md (technical
  detail) in the same session (the user can run `/status`); tidy rather than append. Unclear
  which project: ask.
- **Confidential:** these files, the paper draft and unpublished results stay in this private repo
  (`17thSaint/fzj`) — never publish as Artifacts, push elsewhere, or send to external services.
- Working style: one step at a time; ask questions before structural/workflow changes, show a
  draft, apply on approval. Commit when asked or clearly implied; never push unasked.

## Background

Paper draft *"Fate of a Fractional Chern Insulator under Nonlocal Interactions in Synthetic
Dimensions"* (Geraghty, Nardin, Mazza, Rizzi): `../writing/synth-dims-interactions/` (latest
`...-09-21-26.pdf`). Hardcore bosonic Harper–Hofstadter, ν = 1/2, ρ1D = N/Lx = 1/2, plus an
infinite-range synthetic-column interaction Ui. Increasing Ui connects the Laughlin FCI to a trivial
Tao–Thouless CDW **without closing the gap**: Chern number and TEE stay topological, but pinning
robustness, PES counting (50 → 20) and C⁴(ky) order show triviality. Open threads: how much
nonlocality is needed (finite-range dd ~ 1/(r a)³ → synth-dims), and the gapped path as a
preparation route with Ui tuned by a Stern–Gerlach gradient (→ optimal-control).

## Conventions

- `*daily-things.jl` files are `### header` + `if true/false` blocks run from the user's REPL;
  toggle, don't delete.
- ED test system 4x4 N=2 pbc (also 8x4 N=4, 10x5 N=5); TTN from 12x6 N=6.

## How to run (Claude)

Only as scratchpad scripts (the permission rules allow nothing else without asking):

    cd /home/patrick/fzj/main-git && julia --project=exact-diag <scratchpad>/x.jl

- `--project=synth-dims` for TTN code. QuOCS: `optimal-control/quocs-env/bin/python`
  (configs boot Julia via `quocs_common.py`). Tests: `exact-diag/tests-tevo.jl`.
- Include repo files by **absolute path** (`include` is relative to the script's folder).
  The working dir must be inside `fzj/`: `find_center()`/`get_folder_location()` build paths
  from `pwd()` and fail from the scratchpad.
- PyPlot: `ENV["MPLBACKEND"]="Agg"` before `using PyPlot`; `savefig` PNGs to the scratchpad.

**Data** (`cluster-data/`, 87 GB, gitignored, synced by the user): ED `.jld2` in
`cluster-data/exact-diag/{torus/new-gauge/<study>,time-evo}/`, TTN `.h5` in
`cluster-data/synth-dims/{torus/new-gauge/<study>,excited-states}/`; `<study>` = `dd-ints`,
`ulr-length`, `pinned-scaling`, …. Find/read with `find_data_file(pdict, "ed"|"ttn", dataloc;
file_type)` and `read_data(path) -> (data, meta)`, `dataloc = get_folder_location("cluster-data/...")`.
ED/TTN keys differ (`N`/`particles`, `if_periodic_x`/`if_periodic_phys`).
**Writes need the user's OK every time, naming the files.** `if_save_data` defaults to **true** in
`run_normal_ed` and the TTN runner — always pass `"if_save_data" => false`; no `modify_data`
without approval.
