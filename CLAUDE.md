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
