# Dipolar phase transition — status

Update **Current state**, **Open questions** and **Log** as things change. Last updated: 2026-09-30

## Goal

The paper leaves open *"the degree of nonlocality required"* to drive the FCI → CDW transition.
Replace infinite-range Ui by a finite-range dipolar interaction along the synthetic direction,

    U(r) = intstren / (r a)^3      (scaling "dd", a = magnetic_spacing, intstren = 300)

and find the critical spacing l_c where the CDW sets in (small a → CDW, large a → Laughlin), and
whether it survives the thermodynamic limit. Fixed: ν = 1/2, ρ1D = 1/2, Lx × Lx/2 torus, α = 1/Ly.
Sizes: 6x3, 8x4, 10x5 (ED); 12x6, 14x7, 16x8 (TTN, H100).

## Layout

- `synth-dims/ed-plus-ttns.jl`: ED+TTN analysis; l_c scaling blocks at the end (Hill fit, CDW order
  parameter, boundary-interaction check).
- `synth-dims/long-range-ttn.jl`: TTN runner. `synth-dims/daily-things.jl`: TTN scratch.
- `exact-diag/daily-things.jl`: ED dd blocks from "look at energy spectrum vs magnetic spacing…"
  through "plot all cdw contrasts…" (data collection, l_c from derivatives).
- `exact-diag/observables.jl`: `get_manifold_occupancy` (→ `occs_eig1/2`).
- Data: `cluster-data/{exact-diag,synth-dims}/torus/new-gauge/dd-ints/` — ED 6x3 (42 files), 8x4 (65),
  10x5 (42); TTN 12x6 (12), 14x7 (10), 16x8 (11) as of 2026-09-30.

## Methods for l_c

1. CDW contrast max − min of the manifold density, normalized per size, Hill fit
   `1/(1 + (a/x50)^α)` → l_c = x50; small-a cutoffs `cutoff_magspacs = [1.0, 0.65, 0.8, 0.1, 0.1, 0.1]`.
2. Derivative of the contrast (ED; extra files at a + 0.0002).
3. **Current:** staggered order parameter |Σ occs e^{iπx}| / (Lx Ly), Hill fit, l_c / l_B with
   l_B = sqrt(Ly/2π), fit vs Ly to `A e^{-B Ly} + C`; C = thermodynamic-limit value.

## Current state (2026-09-30)

- **Result: l_c → size-independent C = 0.91 l_B** (sizes 6x3–16x8 without 14x7). The nonlocality
  threshold is set by an intrinsic FQH length (l_B), not system size — answers the paper's open
  question and rules out the boundary-interaction hypothesis (would make l_c size dependent).
- Analysis lives in the "check CDW order parameter vs magnetic spacing…" block, `ed-plus-ttns.jl`
  (committed a24a9b9).
- **14x7 TTN still running.** When it lands: add to the size list, refit, and recheck the
  hand-excluded points (14x7 a=2.1, 12x6 a=1.53).

## Open questions

- [ ] C stable with 14x7?
- [ ] Error bar on C (fit uncertainty; sensitivity to cutoffs and fit window).
- [ ] Contrast vs order-parameter vs derivative l_c consistent?
- [ ] C⁴ / PES vs a at TTN sizes (the paper's more robust diagnostics).

## Gotchas

- Manifold is 2-fold: densities from diagonalizing ⟨Ψi|n|Ψj⟩ (`occs_eig1`), not one eigenvector.
- `occs` layout: ED `[x,y]`, TTN `[y,x]`. TTN files without `"occs"` are skipped.
- The current block reads `magspac_vals`, `orderparam_vals`, `cols`, `lc_vals` from the REPL (its
  loading loop is commented out) — run that loop once first.

## Log

- 2026-05-22 — dd interaction added. 2026-07-13 — fidelity susceptibility along Ui.
- 2026-07-23/28 — "correct" dipolar term tried in TTN code, reverted to plain 1/r³.
- 2026-08-03/12 — derivative and Hill-fit l_c; 12x6 TTN dd runs.
- 2026-09-28 — CDW order parameter, l_c/l_B vs Ly with plateau fit.
- 2026-09-30 — C = 0.91 l_B, size independent; 14x7 still running.
