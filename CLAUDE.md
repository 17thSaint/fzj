# main-git

Research code (Julia ED + TTN, Python QuOCS) for FQH / fractional Chern insulator physics in
synthetic dimensions. Two active projects, each with a living status file. **Read the relevant
one before working, and update it (Current state, results, open questions, Log) when a result,
decision, or bug comes out of the session.**

| Project | Context file | Code |
|---|---|---|
| Optimal control & time evolution (state preparation) | [optimal-control/CONTEXT.md](optimal-control/CONTEXT.md) | `optimal-control/`, `exact-diag/` |
| Phase transition vs. dipolar interactions | [synth-dims/CONTEXT.md](synth-dims/CONTEXT.md) | `synth-dims/` |

If it's unclear which project a task belongs to, ask.

**CONFIDENTIAL.** This file, both CONTEXT.md files, the paper draft and unpublished results are
private, evolving research. They are version-controlled only in this private repo
(`17thSaint/fzj`); never push them anywhere else, never publish them as Artifacts, and never send
them to any external service.

## Shared background

Both projects build on the paper draft *"Fate of a Fractional Chern Insulator under Nonlocal
Interactions in Synthetic Dimensions"* (Geraghty, Nardin, Mazza, Rizzi), latest at
`../writing/synth-dims-interactions/synth-dim-interactions-draft-09-21-26.pdf`.

- Model: hardcore bosonic Harper–Hofstadter on Lx × Ly (x physical, y synthetic), flux α,
  filling ν = 1/2, plus an infinite-range column interaction Ui (all particles in the same
  physical column interact), mimicking synthetic-dimension nonlocality. Fixed ρ1D = N/Lx = 1/2.
- Result: increasing Ui adiabatically connects the Laughlin FCI (Ui = 0) to a trivial
  Tao–Thouless-like CDW **without closing the many-body gap**. Chern number (C = 1/2) and TEE
  (γ ≈ 1/2) stay "topological", but pinning robustness is lost (Δ01 plateaus at ~δ), PES counting
  drops 50 → 20 (10×5, N=5, NA=3), and momentum-space four-point correlator C⁴(ky, 0) orders.
- At commensurate ρ1D = 1 the degeneracy breaks and the system becomes noninsulating (Zeng, Wang,
  Zhai 2015).
- Open threads from the paper: how much nonlocality (e.g. finite-range dd ~ 1/(r a)³ instead of
  infinite range) is needed to drive the transition → **synth-dims project**; the gapped path as a
  finite-time preparation route (CDW → FCI, Ui tuned by a Stern–Gerlach gradient that separates
  synthetic levels) → **optimal-control project**.

## Conventions

- Julia scratch files (`*daily-things.jl`) are sequences of `### header` + `if true/false` blocks;
  toggle, don't delete.
- Standard ED test system: 4x4 N=2 pbc (larger: 8x4 N=4, 10x5 N=5); TTN for 12x6 N=6 and up.

## How to run (for Claude)

Claude runs code only as scripts in its scratchpad (the permission rules allow nothing else
without asking); the user's REPL session and repo files stay the user's.

**Commands & environments**
- Julia 1.12 via juliaup. Use the Project.toml of the code's folder:
  `julia --project=/home/patrick/fzj/main-git/exact-diag <scratchpad>/x.jl` for ED, time evolution
  and optimal-control Julia code; `--project=.../synth-dims` for TTN code. (`other-funcs/`,
  `review-practice-codes/`, `j1j2/` have their own too.)
- Python/QuOCS: `/home/patrick/fzj/main-git/optimal-control/quocs-env/bin/python`; the configs boot
  Julia themselves through `quocs_common.py` (exact-diag project).
- Time-evolution tests: `exact-diag/tests-tevo.jl`.

**Path & plotting gotchas** (verified 2026-09-30)
- `include` resolves relative to the *script's* folder, so a scratchpad script must include repo
  files by absolute path: `include("/home/patrick/fzj/main-git/exact-diag/execute-ed.jl")`.
- `find_center()` / `get_folder_location()` / `include_other_files()` build paths from `pwd()`,
  which must contain a `fzj` folder. From the scratchpad they fail ("Not sure where the center
  is"). Run with the working directory inside the repo (e.g. `cd /home/patrick/fzj/main-git`
  before `julia ...`, or `cd(...)` at the top of the script).
- Plotting is PyPlot. Set `ENV["MPLBACKEND"]="Agg"` before `using PyPlot` and `savefig` PNGs into
  the scratchpad (not the repo's `local-figs/`), then look at them / tell the user the path.

**Data layout** (`cluster-data/`, 87 GB, gitignored, synced from clusters by the user)
- ED: `cluster-data/exact-diag/torus/new-gauge/{dd-ints,ulr-length,pinned-scaling,periodic-potential}/`,
  `.jld2`, named `ed-<key>-<value>-...jld2` from the run's parameters.
- TTN: `cluster-data/synth-dims/torus/new-gauge/{dd-ints,ulr-length,pinned-scaling,bonddim-scaling}/`
  and `cluster-data/synth-dims/excited-states/`, `.h5`, named `ttn-<key>-<value>-...h5`.
- Time evolution: `cluster-data/exact-diag/time-evo/`.
- Find/read with `find_data_file(pdict, "ed"|"ttn", dataloc; file_type="jld2")` and
  `read_data(path) -> (data, metadata)` (other-funcs/data-storage-funcs.jl); dataloc from
  `get_folder_location("cluster-data/...")`. ED and TTN parameter keys differ (e.g. `N` vs
  `particles`, `if_periodic_x` vs `if_periodic_phys`).
- **Writing data: ask the user first, every time**, naming the files. `if_save_data` defaults to
  **true** in both `run_normal_ed` and the TTN runner, so every Claude run must pass
  `"if_save_data" => false` unless the user approved a write; never call `modify_data` without
  approval. Temporary outputs go in the scratchpad.
