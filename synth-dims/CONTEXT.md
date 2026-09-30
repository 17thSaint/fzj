# Dipolar-interaction phase transition — project context

Living status file. Update the **Current state** and **Log** sections whenever something changes.
Last updated: 2026-09-30

## Goal

Follow-up to the paper draft (see root `CLAUDE.md`): the paper uses an **infinite-range** column
interaction Ui and finds an FCI → Tao–Thouless CDW transition without gap closing. Its conclusion leaves open
*"the degree of nonlocality required for the interactions to drive this type of transition"*. This
project answers that with a physical, finite-range **dipole–dipole (dd)** interaction along the
synthetic direction,

    U(r) = intstren / (r * a)^3        r = synthetic separation, a = magnetic_spacing

(`scaling = "dd"` / `scaling_type = "dd"`), intstren = 300 (the paper's large-Ui value), and asks:
**at what spacing a does the CDW set in, and does that critical spacing l_c survive the
thermodynamic limit?** Small a → strongly nonlocal → CDW; large a → interaction dies off → Laughlin.
Experimentally a is set by a Stern–Gerlach gradient separating the synthetic levels (links to the
optimal-control project).

Fixed throughout: hardcore bosons, ν = 1/2, ρ1D = N/Lx = 1/2, Lx × Lx/2 torus (pbc both ways),
α = 1/Ly, hopping_anisotropy = 1.0. Sizes: (6,3,3), (8,4,4), (10,5,5) by ED; (12,6,6), (14,7,7),
(16,8,8) by TTN (GPU, H100).

## Layout

| Where | What |
|---|---|
| `synth-dims/ed-plus-ttns.jl` | Combined ED+TTN analysis; **the l_c finite-size scaling blocks live at the end** (Hill fit, CDW order parameter, boundary-interaction check) |
| `synth-dims/daily-things.jl` | TTN scratch; last block "look at ttn occs under DD ints" (16x8) |
| `synth-dims/long-range-ttn.jl` | TTN model/runner (`scaling = "dd"`, `magnetic_spacing`), data → `cluster-data/synth-dims/torus/new-gauge/dd-ints/` (`.h5`) |
| `exact-diag/daily-things.jl` | ED side: dd spectrum vs spacing (~3584), data collection 8x4 (~3671, ~3812), l_c from derivatives 6x3/8x4/10x5 (~3856–4033), CDW contrast plots (~4034) |
| `exact-diag/execute-ed.jl` | `run_normal_ed` (`scaling_type = "dd"`), data → `cluster-data/exact-diag/torus/new-gauge/dd-ints/` (`.jld2`) |
| `exact-diag/observables.jl` | `get_manifold_occupancy` (diagonalizes the 2×2 manifold density matrix → `occs_eig1/2`), fidelity susceptibility |
| `synth-dims/local-paperstuff/` | Figures for the paper (gitignored) |
| `cluster-codes/` | Slurm/JSC/Jupiter submission and data-sync scripts |

Data on disk (2026-09-30): ED 6x3 (42 files), 8x4 (65), 10x5 (42); TTN 12x6 (12), 14x7 (10), 16x8 (11).

## Methods for locating l_c

1. **CDW contrast** `max(occs) - min(occs)` of the manifold-rotated density (as in the paper's
   Fig. 5), normalized to [0,1] per size, vs a (log axis), fit to a Hill function
   `1 / (1 + (a/x50)^alpha)`; l_c := x50. Small-a points below a per-size cutoff
   (`cutoff_magspacs = [1.0, 0.65, 0.8, 0.1, 0.1, 0.1]`) are dropped.
2. **Derivatives** (ED only): finite difference of the contrast at a ± 1e-4 (`shift_value`,
   extra files at spacing + 0.0002), l_c from the peak of the derivative.
3. **CDW order parameter** (current, uncommitted): `|sum occs(x,y) e^{iπx}| / (Lx Ly)`, same Hill
   fit, then l_c / l_B with l_B = sqrt(Ly / 2π), fit vs Ly to `A e^{-B Ly} + C`; the plateau C is
   the thermodynamic-limit estimate.
4. Checked earlier: energy spectrum vs spacing (gap stays open, like the Ui path), fidelity
   susceptibility along Ui.

## Current state (2026-09-30)

- Working in the **"check CDW order parameter vs magnetic spacing for various system sizes"**
  block at the end of `ed-plus-ttns.jl` (uncommitted). Moved from contrast to the staggered
  order parameter, and to l_c / l_B with an exponential-plateau fit in Ly.
- **Main result: l_c converges to a size-independent value of about one magnetic length.**
  Fitted plateau of l_c / l_B vs Ly (`A e^{-B Ly} + C`, sizes 6x3–16x8 without 14x7):
  **C = 0.91 l_B**.
- **Interpretation:** the amount of nonlocality needed to destroy the FCI is set by an intrinsic
  length of the FQH state (l_B), not by the system size. This answers the paper's open question
  and rules out the earlier hypothesis (block below it, now commented out) that the transition is
  controlled by the interaction strength at the cylinder boundary, which would have made l_c
  size dependent.
- **14x7 TTN data still running** — add it to the order-parameter block's size list and refit
  when it lands (two TTN points are hand-excluded: 14x7 a=2.1, 12x6 a=1.53; recheck whether the
  new 14x7 set still needs that).

## Open questions / next steps

- [ ] Add 14x7 when the runs finish; check C = 0.91 l_B is stable with it.
- [ ] Error bar on C (fit uncertainty, sensitivity to the per-size `cutoff_magspacs` and to the
      Hill-fit window).
- [ ] Is the contrast-based and order-parameter-based l_c the same? Cross-check against the
      derivative method on ED sizes.
- [ ] Beyond real-space order: momentum four-point correlator C⁴ / PES counting vs a (the
      paper's more robust diagnostics) at the TTN sizes.

## Conventions & gotchas

- Groundstate manifold is 2-fold; densities must come from diagonalizing the 2×2 matrix
  `<Ψi|n|Ψj>` (`occs_eig1`), not from a single eigenvector.
- Occupation array layout differs: ED `occs[x,y]`, TTN `occs[y,x]`.
- TTN files missing `"occs"` in metadata are skipped (printed "No occupancy data").
- The current order-parameter block reads `magspac_vals`, `orderparam_vals`, `cols`, `lc_vals`
  from the REPL: its data-loading loop is commented out, so it only runs after that loop has been
  run once in the session.
- ED data files are keyed by `magnetic_spacing` in the filename (10-digit precision); TTN by the
  `ttn-...-magnetic_spacing-X-...h5` pattern with `lr` = interaction range.

## Log

- 2026-05-22 — dd interaction added.
- 2026-07-13 — ED fidelity susceptibility along Ui.
- 2026-07-23/28 — tried the "correct" dipolar term (m = m' fix) in the TTN code, reverted to plain 1/r³.
- 2026-08-03 — first derivative-based l_c from manifold occupations.
- 2026-08-12 — Hill-function l_c fits; 12x6 TTN dd runs on H100.
- 2026-09-04/08 — dd + magnetic-gradient time evolution (→ optimal-control project).
- 2026-09-28 — analyzing dd critical spacing: CDW order parameter, l_c/l_B vs Ly with plateau fit.
- 2026-09-30 — Result: l_c plateaus at C = 0.91 l_B, size independent → nonlocality threshold is
  an intrinsic FQH length scale. 14x7 TTN still running.
