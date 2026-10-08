# Dipolar phase transition — status

Update **Current state**, **Open questions** and **Log** as things change. Last updated: 2026-10-08

## Goal

The paper leaves open *"the degree of nonlocality required"* to drive the FCI → CDW transition.
Replace infinite-range Ui by a finite-range interaction along the synthetic direction,

    U(r) = V / (r s)^p      (p = 3: scaling "dd", s = magnetic_spacing, V = intstren = 300 or 50)

and find the critical point (small s → CDW, large s → Laughlin) and whether it survives the
thermodynamic limit. Scaling variable Γ = V / s^p (nearest-neighbour synthetic coupling). Fixed: ν = 1/2.
Torus: ρ1D = 1/2, Lx × Lx/2, α = 1/Ly; sizes 6x3, 8x4, 10x5 (ED); 12x6, 14x7, 16x8 (TTN, H100).
Open y (OBC / cylinder, ED): 6×Ly N=3 with Ly = 4–7, 8x4 N=4, V = 50.

## Current state (2026-10-08)

- **p = 3 (dd) is done for now.** Γ_c grows strongly with Ly on both torus and OBC:
  - The earlier torus "l_c → 0.91 l_B, size independent" was wrong; against Γ, Γ_c diverges with size.
  - OBC 6×Ly N=3, V=50: block "cdw transition for 6xi N=3 scaling…" in `exact-diag/daily-things.jl`.
    S(π) and Binder vs Γ, min-max normalized, fitted with `1 - hill_model` (Γ_c = x50) from the
    pre-transition dip upward, s ≥ 0.2, 24 points per Ly:

    | | 6x4 | 6x5 | 6x6 | 6x7 |
    |---|---|---|---|---|
    | Γ_c S(π) | 13.4 ± 0.4 | 26.6 ± 0.5 | 78.6 ± 0.8 | 200 ± 4 |
    | Γ_c Binder | 11.3 ± 0.4 | 24.2 ± 0.4 | 75.1 ± 0.7 | 194 ± 3 |
    | s_c / l_B | 2.2 | 1.5 | 0.96 | 0.64 |

    Γ_c ratio per unit Ly ~2, 3, 2.6. Exponential e^{~1.0 Ly} fits better than a power law (∝ Ly^5),
    but there are only 4 points. Smooth crossover, not a jump. Errors are fit errors only. A
    free-baseline 4-parameter sigmoid gives Γ_c within ~12% (6x4) and ~6% (larger).
  - Interpretation (user): 1/r³ is short-ranged in a quasi-1D system.
- **Caveat, deferred:** at fixed Lx=6, N=3 the flux (Nφ = 6) is fixed. Growing Ly only refines
  the lattice (α = 1/(Ly−1) → 0) and adds synthetic states; it is not a thermodynamic limit. In
  experiment Ly (number of synthetic states) is fixed and Lx → ∞. Planned check: fine-grid 8x4
  and 10x4 N=5 (s ∈ [1.0, 2.6]) to see if Γ_c(Lx) at Ly=4 is converged. A proper Binder crossing
  is better done with TTN on open strips (ED N ≤ 6 too small, see Gotchas).

## Next: other power laws p (start here)

Goal: Γ_c(Ly) for longer-ranged p (e.g. p = 1, 2; maybe 0.5) on the same cheap OBC 6×Ly N=3,
V=50 series. Does the growth of Γ_c with Ly track p, and does any p > 0 give a finite Γ_c?
- **Hypothesis (Claude, untested):** dd acts within a synthetic column and penalizes two particles
  in the same column. For any p > 0 that penalty decays with their synthetic separation, which
  grows with Ly. Naively Γ_c ∝ Ly^p, so only p = 0 (the paper's Ui) would be finite. The p = 3
  data grows faster than Ly³, so the naive picture is incomplete.
- **Code:** the profile is `long_range_scaling` in `other-funcs/basic-2d-stuff.jl` (branch
  `"dd"`: `strengths[x+1] = V / (a x)^3`, index = synthetic distance, `strengths[1] = V` onsite,
  irrelevant for hardcore bosons). Plumbing to touch: `make_filename_dict` and
  `get_normal_model_params_ed` in `exact-diag/execute-ed.jl` (the kwargs go through
  `other_params_dict`, and `"dd"` sets the default `dataloc` to torus `dd-ints`), the filename
  dict in `exact-diag/time-evolution.jl`, and `synth-dims/long-range-ttn.jl` (`"dd"` branches) for
  TTN. `"magnetic_gradient"` in `long_range_scaling` also hard-codes ³.
- **Decide first (ask the user):** a new scaling type (e.g. `"power"` + exponent key) or a power
  key on `"dd"` that is written to filenames only when ≠ 3, so the existing dd files keep
  matching. Also where the data goes (e.g. `cluster-data/exact-diag/obc/power-ints/`).
- Runner: OBC `if false` block in `exact-diag/execute-ed.jl` (loops `xis`; saves states plus
  `column_cdw_statistics`). Analysis: copy the block above with Γ = V/s^p. The s window for the
  jump has to be found again for each p (coarse grid first, then fine).

## Methods for the critical point

1. CDW contrast max − min of the manifold density, Hill fit `1/(1 + (s/x50)^α)` → l_c = x50.
2. Derivative of the contrast (ED; extra files at s + 0.0002).
3. Staggered order parameter |Σ occs e^{iπx}| / (Lx Ly), Hill fit, l_c / l_B (l_B = sqrt(Ly/2π))
   fit vs Ly to `A e^{-B Ly} + C` → superseded.
4. Same order parameter vs Γ, rising Hill fit, Γ_c vs Ly (torus, V = 300, 50).
5. **Current (open y):** from column correlations averaged over the two lowest states:
   S(π) = ⟨|Σ_x e^{iπx} n_x|²⟩/N and the Binder cumulant of the staggered column magnetization.
   On the cylinder ⟨n_x⟩ is uniform, so 3./4. don't work there.

## Layout

- `exact-diag/daily-things.jl`: ED torus dd blocks; at the end the OBC/cylinder blocks ("look at
  CDW transition for Gamma with cylinder or open boundary conditions", "cdw transition for 6xi N=3
  scaling…" = current Γ_c fit and plots).
- `exact-diag/execute-ed.jl`: `if false` runner block for the OBC/cylinder sweeps.
- `exact-diag/observables.jl`: `get_manifold_occupancy` (→ `occs_eig1/2`); `column_cdw_statistics`
  (⟨n_x n_x'⟩, P(k on even columns)) and `cdw_order_parameters` (S(π), correlation ratio, Binder).
- `synth-dims/ed-plus-ttns.jl`: ED+TTN torus analysis; scaling blocks at the end (s_c/l_B plateau,
  Γ with V = 300, 50, boundary-interaction check), all `#= =#`.
- `synth-dims/dd-paper-plots.jl`: paper figures (so far: shifted Zeeman minima + 1/(s r)³ profile).
- `fpeps/execute-fpeps.jl`: fPEPS (QuantumNaturalfPEPS) for dd on OBC, 8x4 V=50 s=0.5, bond dim 2;
  demo only (10 iterations), goal is larger open systems.
- `synth-dims/long-range-ttn.jl`: TTN runner. `synth-dims/daily-things.jl`: TTN scratch.
- Data (dd-ints, p = 3), as of 2026-10-08:
  - ED torus `cluster-data/exact-diag/torus/new-gauge/dd-ints/`: V=300: 6x3 (42), 8x4 (65),
    10x5 (42); V=50: 6x3 (21), 8x4 (9), 10x5 (1 file holding the whole sweep).
  - TTN torus `cluster-data/synth-dims/torus/new-gauge/dd-ints/`: V=300, 12x6, 14x7, 16x8 (12 each).
  - ED `cluster-data/exact-diag/obc/dd-ints/`: V=50, 6×{4,5,6,7} N=3 (26 each: coarse s = 0.1–20
    plus a fine grid across the jump), 8x4 (11, coarse).
  - ED `cluster-data/exact-diag/cylinder/dd-ints/`: V=50, 6×{4,5,6} N=3 (15 each), 6x7 (11), 8x4 (22).

## Open questions

- [ ] Γ_c(Ly) for p = 1, 2, …: does the growth track p; is any p > 0 finite?
- [ ] Γ_c converged in Lx at fixed Ly (8x4, 10x4)? Proper Binder crossing in Lx with TTN?
- [ ] What is the pre-transition dip in S(π)/Binder (strong in 6x4)?
- [ ] V = 50 TTN data at 12x6+ for the torus?
- [ ] Does fPEPS reach useful open sizes / accuracy?
- [ ] C⁴ / PES vs s at TTN sizes (the paper's more robust diagnostics).

## Gotchas

- Manifold is 2-fold: torus densities from diagonalizing ⟨Ψi|n|Ψj⟩ (`occs_eig1`); open-y order
  parameters average the column statistics of states 1 and 2 (basis independent).
- Binder at fixed N has a disordered baseline 2/(3N) (0.222 for N=3, 0.167 for N=4); CDW plateau
  for perfect TT is 2/3 (S(π)/N → 1). With N ≤ 6 the offset distorts crossings. Growing Ly at fixed
  N doesn't lengthen the ordering direction, so there's no Binder crossing in Ly.
- The OBC CDW plateau is far from perfect TT and drops with Ly (S(π)/N 0.71, 0.53, 0.44, 0.40).
- OBC 6x6, 6x7 dip at s ≤ 0.15 and cylinder 8x4 is unconverged for s < 0.25 (|H| > 1e4,
  Lanczos); the analysis uses s ≥ 0.2 / 0.25. For p ≠ 3 the large-|H| limit moves (V/s^p).
- Open y: α = N / (ν Lx' (Ly − 1)), so 6×Ly N=3 files carry α = 1/(Ly−1); on OBC α drifts with Lx
  at fixed Ly (0.40, 0.38, 0.37 for Lx = 6, 8, 10 at Ly=4).
- `occs` layout: ED `[x,y]`, TTN `[y,x]`. TTN files without `"occs"` are skipped.
- ED torus 10x5 V=50 is one file with `all_magspacs` / `all_occs_eig1`, not one file per s.
- OBC 8x4 s = 0.1, 0.17 have no `column_corr*` in metadata.
- The "look at CDW transition for Gamma…" block computes `column_corr*` if missing and **writes
  them to the file metadata** with `modify_data` — ask first.
- Hand-excluded TTN points: 14x7 s = 20.0, 12x6 s = 1.53.

## Log

- 2026-05-22 — dd interaction added. 2026-07-13 — fidelity susceptibility along Ui.
- 2026-07-23/28 — "correct" dipolar term tried in TTN code, reverted to plain 1/r³.
- 2026-08-03/12 — derivative and Hill-fit l_c; 12x6 TTN dd runs.
- 2026-09-28/30 — CDW order parameter, l_c/l_B vs Ly, apparent C = 0.91 l_B plateau (later wrong).
- 2026-10-02/05 — 14x7 TTN in, error bars on l_c, paper-style figure; `dd-paper-plots.jl` started.
- 2026-10-05/06 — Γ = V/s³ axis with V = 300, 50: Γ_c diverges with size (1/r³ short-ranged in quasi-1D).
- 2026-10-07 — OBC/cylinder ED dd data; S(π), correlation ratio, Binder from column correlations.
- 2026-10-08 — OBC fine Γ grid: Γ_c(Ly=4..7) ≈ 13, 27, 79, 200 (S(π)), ~e^{Ly}; fit added to the
  block. Lx check deferred; next: other power laws p.
