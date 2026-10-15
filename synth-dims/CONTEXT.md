# Dipolar phase transition — status

Update **Current state**, **Open questions** and **Log** as things change. Last updated: 2026-10-15

## Goal

The paper leaves open *"the degree of nonlocality required"* to drive the FCI → CDW transition.
Replace infinite-range Ui by a finite-range interaction along the synthetic direction,

    U(r) = V / (r s)^p      (p = 3: scaling "dd", s = magnetic_spacing, V = intstren = 300 or 50)

and find the critical point (small s → CDW, large s → Laughlin) and whether it survives the
thermodynamic limit. Scaling variable Γ = V / s^p (nearest-neighbour synthetic coupling). Fixed: ν = 1/2.
Torus: ρ1D = 1/2, Lx × Lx/2, α = 1/Ly; sizes 6x3, 8x4, 10x5 (ED); 12x6, 14x7, 16x8 (TTN, H100).
Open y (OBC / cylinder, ED), V = 50: 6×Ly N=3 and 8×Ly N=4 with Ly = 4–7, 10x4 N=5.
Other exponents: scaling `"powerlaw"` with key `"exponent"` = p (p = 3 equals `"dd"`).

## Current state (2026-10-14)

**Entanglement / central charge (2026-10-15, side thread):** vertical cuts (A = all sites with
x ≤ ℓ, all y) on the cylinder (x periodic, y open), Calabrese–Cardy fit
S = (c/3) log[(Lx/π) sin(πℓ/Lx)] + const over ℓ ≤ Lx/2 (S(ℓ) = S(Lx − ℓ) for the pure GS).
Flat Ui = 0, Ly = 4: 6x4 N=3 c = 3.2 ± 0.1, 8x4 N=4 c = 2.7 ± 0.07 (3 resp. 4 independent
points). Not a CFT value: the fit is dominated by the area-law offset at Lx ≤ 8; the 6x4 GS is
gapped (E2 − E1 = 0.10; 8x4: 0.016). Block "try Calabrese-Cardy scaling cut size…" at the end
of `exact-diag/daily-things.jl` (reads states, saves `ee_cuts`, fits c with LsqFit).

**Bottom line:** OBC ground-state data for p = 0 (flat), 0.5, …, 3 at 6×Ly N=3 and 8×Ly N=4
(Ly = 4–7), plus 10x4 N=5. All p go from the same Laughlin state to the same large-Γ state;
below Γ_c, p ≳ 1.5 takes a detour (precursor), no new phase seen. Γ_c grows with Ly for every p
(including p = 0, i.e. partly geometry) at both Lx, and faster for larger p. At Ly = 4, Γ_c does
not grow with Lx (6, 8, 10). Hints that the exponent changes the physics: small p is far from a
single power law in Ly (curving down), p ≥ 2 is close to a clean power law, and b(p) − b(0) keeps
growing past p = 1 with no saturation. But only 4 Ly, the GS CDW signal collapses with Ly at fixed
N (S(π)/N rise ≈ 0.4, 0.15, 0.06, 0.035 for Ly = 4–7, at both Lx), and Hill vs midpoint Γ_c
disagree at Ly ≥ 6. User (2026-10-14): these are good limited-size studies, larger sizes are
needed to resolve the p dependence (see Next).

- **8×Ly N=4 (cluster sweep, 2026-10-14):** block "cdw transition for 8xi N=4 scaling…" in
  `exact-diag/daily-things.jl` (copy of the 6×Ly block, reads the `all_*` metadata, GS only).
  GS Hill Γ_c from S(π), 6×Ly in brackets; b ± its standard error from the residuals of the
  log-log line (Ly = 4–7):

    | Γ_c GS, Hill | 8x4 | 8x5 | 8x6 | 8x7 | b |
    |---|---|---|---|---|---|
    | p = 0 | 1.95 (1.94) | 3.33 (2.78) | 3.83 (3.41) | 4.42 (3.98) | 1.43 ± 0.30 (1.28) |
    | p = 0.5 | 2.63 (2.75) | 5.04 (4.23) | 6.13 (5.35) | 7.36 (6.43) | 1.80 ± 0.33 (1.51) |
    | p = 1 | 3.62 (4.0) | 7.63 (6.31) | 9.91 (8.4) | 12.6 (10.7) | 2.19 ± 0.34 (1.74) |
    | p = 1.5 | 5.37 (5.86) | 11.5 (9.29) | 16.0 (13.4) | 22.1 (19.0) | 2.48 ± 0.27 (2.09) |
    | p = 2 | 8.01 (8.6) | 16.9 (13.8) | 26.0 (22.0) | 39.0 (32.6) | 2.80 ± 0.16 (2.39) |
    | p = 2.5 | 12.0 (12.8) | 25.1 (20.4) | 42.7 (35.2) | 70.2 (56.9) | 3.14 ± 0.06 (2.68) |
    | p = 3 | 17.2 (18.9) | 37.0 (29.7) | 70.5 (56.3) | 127 (103) | 3.57 ± 0.06 (3.05) |

  - Ly ≥ 5: 8×Ly Γ_c is 10–25% above 6×Ly; still grows with Ly for every p, so not an N = 3
    artifact (Lx = 8 also changes α and Nφ).
  - Ly = 4, Lx = 6, 8, 10, Γ_c S(π) (Binder): p = 0: 1.94, 1.95, 1.73 (1.55, 1.70, 1.56);
    p = 1: 4.0, 3.62, 3.34; p = 2: 8.6, 8.01, 8.14; p = 3: 18.9, 17.2, 17.7 (17.4, 16.6, 16.3).
    Flat or slightly falling with Lx.
  - Ly = 4 transitions are steeper than the Hill shape (the fit misses the upper tail), as at 6x4.
- **Single power law Γ_c ∝ Ly^b is a poor fit (local exponents, 2026-10-14):** at 8×Ly the 4→5
  step has slope 2.4–3.4 for every p, then the slope drops, so Ly = 4 dominates b. The b error
  bars show it: ±0.3 for p ≤ 1.5, ±0.06 for p ≥ 2.5. b over Ly = 5–7 only (S(π)): Hill 0.84,
  1.12, 1.49, 1.94, 2.47, 3.05 for p = 0–2.5 (b < 1 for small p); midpoint 1.26, 1.69, 2.09, 2.52,
  3.00, 3.55 (never < 1). A local slope carries ±0.2–0.3, so "b crosses 1" isn't resolved. Both
  give b(p) − b(0) still growing past p = 1 (midpoint ≈ 0.85 p): no saturation or kink at p = 1,
  where Σ r^−p starts converging in 1D. The naive picture (pair spread over Ly feels
  Γ Ly^−min(p,1)) predicts saturation; the data don't show it, which fits the pairs separating in
  y (precursor). Hill vs midpoint differ mainly in the baseline (dip minimum vs Γ → 0 value).
- **6×Ly N=3 GS (2026-10-09):** Hill Γ_c from S(π) (Binder in brackets):

    | Γ_c GS, Hill | 6x4 | 6x5 | 6x6 | 6x7 | b |
    |---|---|---|---|---|---|
    | p = 0 | 1.94 (1.55) | 2.78 (2.47) | 3.41 (3.27) | 3.98 (3.88) | 1.28 (1.65) |
    | p = 0.5 | 2.75 (2.32) | 4.23 (3.78) | 5.35 (5.14) | 6.43 (6.29) | 1.51 (1.79) |
    | p = 1 | 4.0 (3.49) | 6.31 (5.68) | 8.4 (8.1) | 10.7 (10.4) | 1.74 (1.96) |
    | p = 1.5 | 5.86 (5.26) | 9.29 (8.37) | 13.4 (12.9) | 19.0 (18.6) | 2.09 (2.26) |
    | p = 2 | 8.6 (7.86) | 13.8 (12.4) | 22.0 (21.2) | 32.6 (31.9) | 2.39 (2.53) |
    | p = 2.5 | 12.8 (11.8) | 20.4 (18.4) | 35.2 (34.0) | 56.9 (56.0) | 2.68 (2.82) |
    | p = 3 | 18.9 (17.4) | 29.7 (26.9) | 56.3 (54.5) | 103 (102) | 3.05 (3.18) |

  - Midpoint b (GS, S(π)): 1.58, 1.94, 2.24, 2.50, 2.84, 3.15, 3.43 for p = 0–3. Local b: Hill
    p = 0 1.6 → 1.1 → 1.0 (concave), p = 3 2.0 → 3.5 → 3.9 (convex); midpoint convex for all p.
  - Energies are sorted in all 464 saved OBC 6×Ly N=3 files, so state 1 is the GS.
- **GS only, not the 2-state average (`if_gs_only`):** on open y, state 2 is an excited state for
  Ly ≥ 5 (E2 − E1 ≈ 0.03–0.09 at Γ → ∞, checked at Lx = 6 only), so the average describes no
  physical state and borrows CDW signal from it: Γ_c inflated (6x7 p = 3: 200 vs 103 GS).
  At Ly = 4 the pair is near-degenerate only at large Γ: the plateau agrees but Γ_c doesn't
  (p = 0: 1.94 GS vs 1.40 averaged at 6x4, 1.95 vs 1.57 at 8x4), since below Γ_c state 2 isn't the
  CDW partner and 6x4 has true GS level crossings. S(π)/Binder come from ⟨n_x n_x'⟩, so they see
  the order in a cat-like GS. The 6×Ly block's toggle is currently `false` (user's setting), the
  8×Ly block's `true`.
- **2-state-average results (2026-10-08, superseded, numbers in git history):** p = 3 6×Ly Γ_c
  13.4, 26.6, 78.6, 200 (S(π)); b ≈ 1.1 p + 1.6; midpoint Γ_c 10–20% above Hill.
- **Same physics for all p? (2026-10-08, scratchpad, dense no-save paths Γ = 0.03–3000, p = 0–3,
  6×Ly):**
  - Endpoints identical: S(π) at Γ → 0 (Laughlin) and Γ → ∞ (CDW) agree across p to 3–4 digits;
    above Γ_c raw S(π), Binder, E3−E2 collapse onto one curve against Γ/Γ_c.
  - Below Γ_c they differ: larger p gives a pre-transition dip in S(π)/Binder below the Laughlin
    value, which grows with Ly, while the same-column pairs move apart in y (mean separation up).
    For p ≲ 1 there is (almost) no dip.
  - Ground state: smooth path for Ly ≥ 5 (neighbour fidelity ≥ 0.993), one broad GS
    fidelity-susceptibility bump at Γ/Γ_c ≈ 0.3–1, height weakly growing with p; E3−E1 dips but
    stays open. Best GS overlap of the p = 3 path with the p = 0 path: 0.95, 0.94, 0.91
    (Ly = 5–7), at Γ/Γ_c ≈ 0.15–0.3; p = 0.5: ≥ 0.993.
  - Ly = 4 has true GS level crossings (near Γ → 0 and near Γ_c).
  - GS-only, p = 1: Γ_c(1)/Γ_c(0) = 2.06, 2.27, 2.46, 2.69 (∝ Ly^0.5).
  - Reading (Claude): same Laughlin → CDW crossover for all p, but for p ≳ 1.5 the system first
    avoids the interaction by separating same-column pairs in y before the CDW forms.
- **Caveat:** at fixed N growing Ly keeps the flux fixed (6×Ly: Nφ = 6, α = 1.2/(Ly−1) → 0), so
  it refines the lattice and adds synthetic states; it is not a thermodynamic limit. In
  experiment Ly (number of synthetic states) is fixed and Lx → ∞.

## Next (start here)

**Finite-size scaling of the fitted c (2026-10-15, next).** Follow the comment above the open
block: put the per-size cut scans together (one loop over sizes) and scale the fitted c with
Lx, for Laughlin (Ui = 0) and the CDW (large Ui), to see where c lands (expect c → 0 in the
gapped CDW; Laughlin on a thin cylinder may also be gapped, or show c = 1 from the two edges).
- Data now: only Ui = 0, 6x4 N=3 and 8x4 N=4 in `cylinder/flat-ints`. CDW points and 10x4 N=5
  need new ED runs (saving = write, ask; the user makes folders).
- 10x4 needs a better EE route: dense ρ_A at ℓ = 5 is 21700² (7.5 GB). ρ_A is block-diagonal in
  N_A, and the Schmidt values come from SVDs of the blocks Ψ_{N_A} (C(20,k) × C(20,5−k) at
  10x4), so a block-SVD version of `compute_reduced_density_matrix` makes every cut cheap.
- With 3–4 independent ℓ per size the fit is crude; TTN/MPS would give longer Lx.

**Larger sizes to resolve the p dependence (decision open, 2026-10-14).** User's idea: fPEPS
(`fpeps/execute-fpeps.jl`) for more system sizes. Discussed first:
- Which size family? b(p) so far is at fixed N, Ly growing. Bigger sizes in that family don't
  need fPEPS: ED reaches 8x8–8x9 N=4 (0.6–1M states) and 6×Ly N=3 to Ly ≈ 10 (C(60,3) ≈ 34k,
  minutes no-save), but the GS signal keeps shrinking there. New methods are needed for fixed
  density with N growing: (a) fixed Ly, Lx growing (experiment-like; signal holds with Lx at
  Ly = 4; gives a real Binder crossing in Lx, removing the Hill/midpoint ambiguity), or (b) both
  growing.
- Claude's suggestion: MPS/TTN on open strips for (a). The interaction acts only within a
  synthetic column, so with x as the long direction it's local in x for every p (small MPO);
  `"powerlaw"` TTN code exists (`synth-dims/long-range-ttn.jl`). fPEPS for (b), after a benchmark
  against ED.
- fPEPS risks: chiral topological states are hard for finite-D PEPS while the CDW is near
  product, so low D biases Γ_c down, possibly in a Ly/p-dependent way (needs D convergence of
  Γ_c); the demo uses `siteinds("Boson",…)` with no U(1) or chemical potential, so N isn't fixed
  unless QuantumNaturalfPEPS conserves it in sampling (not checked), and hardcore needs site dim
  2; Γ = 0.1–1000 vs t = 1 is stiff. In its favour: diagonal long-range terms are cheap in VMC.
  First step for fPEPS would be an 8x4/8x5 ED benchmark (N, hardcore, D = 2–6, Γ_c).

Smaller options: Γ_c(Lx) figure at Ly = 4 (6x4, 8x4, 10x4); b(p) − b(0) / local-exponent panel
without Ly = 4; save energies (`all_nrgs`) in future sweeps to check E2 − E1 at Lx = 8.

**Lx sweep (done, 2026-10-14), for reuse:** runner `#= … if false` block "Γ sweeps on open
strips…" in `exact-diag/execute-ed.jl`, `julia execute-ed.jl "which_one" $SLURM_ARRAY_TASK_ID`,
`--array=0-34` (task = size × p, sizes outer: 0–6 = 8x4, …, 28–34 = 10x4), or
`"Lx" "Ly" "particles" ["exponent"]` by hand. Γ grid: 13 log points over 0.1–1000 + 12 in
[Γ_c/6, 10 Γ_c] around the 6×Ly GS Hill Γ_c.
- Returned data: one file per (size, p) holding the whole sweep in its metadata (see Layout).
  The repo runner's comments still describe per-point files.
- `if_reading = true`: the first job of a size builds + saves the undressed matrices to
  `cluster-data/exact-diag` while holding a lock file `building-undressed-N-…-Lx-…-Ly-…-obc.lock`
  (exclusive create), the others wait. A job killed while building leaves the lock: delete it
  before resubmitting. Laptop build times: 8x6 ~2 min, 10x4 ~30 min, 8x7 > 80 min. Per point
  (matrices read): 10x4 ~5 min, 8x7 ~3–4 min, 8x6 ~1.5 min; memory a few GB.
- Reading path (`getHamiltonian`) equals `buildHam` exactly (max |ΔH| = 0) at 8x4–8x6 for flat,
  p = 1, 3 and dd; interaction files have no BC in the name (torus ones are reused).

Notes for any new p:
- **Code (done):** `"powerlaw"` in `long_range_scaling` (`other-funcs/basic-2d-stuff.jl`,
  `strengths[x+1] = V / (s x)^p`), `"exponent"` in `get_normal_model_params_ed` and in filenames
  via `make_filename_dict` (`exact-diag/execute-ed.jl`). Also in time evolution (`timeham`
  rebuild, tevo filename, `get_critical_dt_at_spacing`), the TTN runner
  (`synth-dims/long-range-ttn.jl`: filename + `model_paras_dict`), `OPTIMIZATION_PROFILE_KEYS`,
  the intstren-ramp AD dressing and the QuOCS intstren-ramp configs (default 3.0). Missing
  `"exponent"` errors in `long_range_scaling` (no default). `"magnetic_gradient"` still hard-codes
  1/a³, and there's no default `dataloc` for `"powerlaw"` (ED and TTN fall back to `ulr-length`).
- **Recipe used for p = 0.5–2.5:** a no-save coarse scan (25–31 log points, Γ = 0.01–1000, ~6 min
  for 3 p × all Ly; reproduces the saved dd Γ_c to 1–4%) locates Γ_c, then the saved grid is Γ = 0.1, 1, 10,
  100, 1000 plus 10 log points in [Γ_c/6, 10 Γ_c] (15 files per Ly). Data goes to
  `cluster-data/exact-diag/obc/powerlaw-<p>` with `.` → `p` (folders made by the user), the analysis
  block finds it via `power_folder(p)`; add the new p to `ps`. Claude can't `mkdir` there
  (permission rules), the user makes the folders.

## Methods for the critical point

1. CDW contrast max − min of the manifold density, Hill fit `1/(1 + (s/x50)^α)` → l_c = x50.
2. Derivative of the contrast (ED; extra files at s + 0.0002).
3. Staggered order parameter |Σ occs e^{iπx}| / (Lx Ly), Hill fit, l_c / l_B (l_B = sqrt(Ly/2π))
   fit vs Ly to `A e^{-B Ly} + C` → superseded.
4. Same order parameter vs Γ, rising Hill fit, Γ_c vs Ly (torus, V = 300, 50).
5. **Current (open y):** from column correlations of the GS (`column_corr1`, `evencolumn_dist1`;
   until 2026-10-08 averaged over the two lowest states, `if_gs_only = false`):
   S(π) = ⟨|Σ_x e^{iπx} n_x|²⟩/N and the Binder cumulant of the staggered column magnetization.
   On the cylinder ⟨n_x⟩ is uniform, so 3./4. don't work there. Γ_c = x50 of a rising Hill fit
   to the min-max normalized curve from its minimum (the pre-transition dip) upward.
6. **Midpoint Γ_c (fit-free):** Γ where raw S(π) (or Binder) last crosses (S(Γ→0) + S(Γ→∞))/2,
   log-interpolated. Both endpoints are p-independent, so this doesn't depend on the dip.
7. **Path checks (2026-10-08, scratchpad only, not in the repo):** no-save ED (nev = 20) along
   Γ = 0.03–3000, 41 log points (81 for the reference paths p = 0, 3), Ly = 4–7, ~25 min in
   total. Per point: the two lowest states, 6 lowest energies, S(π), Binder, P(two particles in the
   same column) and their mean synthetic separation (from `full_basis` + `coordinate`).
   - GS best-match fidelity: for each state on path p, max over the reference path's grid of
     |⟨ψ1|φ1⟩|². The neighbour fidelity on the reference path is the grid floor.
   - GS fidelity susceptibility: (1 − |⟨ψ1(Γi)|ψ1(Γi+1)⟩|²) / Δ(ln Γ)².
   - Use GS-only quantities: the 2-state manifold version picks up excited-state crossings.
   Scripts and `paths.jld2` (290 MB) are in the session scratchpad
   (`/tmp/claude-1000/-home-patrick-fzj-main-git/884a0ff2-…/scratchpad/`: `paths.jl`,
   `paths-gs.jl`, `paths-analysis.jl`, `raw.jl`) and may be gone; rebuild from this description.

## Layout

- `exact-diag/daily-things.jl`: ED torus dd blocks; at the end the OBC/cylinder blocks ("look at
  CDW transition for Gamma with cylinder or open boundary conditions", "cdw transition for 6xi N=3
  scaling…", "cdw transition for 8xi N=4 scaling…" = current Γ_c fits and plots; second figure
  Γ_c(Ly) and b(p) with error bars).
- `exact-diag/execute-ed.jl`: `if false` runner block for the OBC/cylinder sweeps.
- `exact-diag/observables.jl`: `get_manifold_occupancy` (→ `occs_eig1/2`); `column_cdw_statistics`
  (⟨n_x n_x'⟩, P(k on even columns)) and `cdw_order_parameters` (S(π), correlation ratio, Binder).
- `synth-dims/ed-plus-ttns.jl`: ED+TTN torus analysis; scaling blocks at the end (s_c/l_B plateau,
  Γ with V = 300, 50, boundary-interaction check), all `#= =#`.
- `synth-dims/dd-paper-plots.jl`: paper figures (so far: shifted Zeeman minima + 1/(s r)³ profile).
- `fpeps/execute-fpeps.jl`: fPEPS (QuantumNaturalfPEPS) for dd on OBC, 8x4 V=50 s=0.5, bond dim 2;
  demo only (10 iterations), goal is larger open systems.
- `synth-dims/long-range-ttn.jl`: TTN runner. `synth-dims/daily-things.jl`: TTN scratch.
- Data:
  - ED torus `cluster-data/exact-diag/torus/new-gauge/dd-ints/`: V=300: 6x3 (42), 8x4 (65),
    10x5 (42); V=50: 6x3 (21), 8x4 (9), 10x5 (1 file holding the whole sweep).
  - TTN torus `cluster-data/synth-dims/torus/new-gauge/dd-ints/`: V=300, 12x6, 14x7, 16x8 (12 each).
  - ED `cluster-data/exact-diag/obc/dd-ints/`: V=50, 6×{4,5,6,7} N=3 (26 each: coarse s = 0.1–20
    plus a fine grid across the jump), 8x4 (11, coarse).
  - ED `cluster-data/exact-diag/obc/powerlaw-{0p5,1p0,1p5,2p0,2p5}/`: p = 0.5–2.5, V=50,
    6×{4,5,6,7} N=3 (15 each, see Next), `column_corr*`/`evencolumn_dist*` in metadata, 201 MB per
    folder.
  - ED `cluster-data/exact-diag/obc/flat-ints/`: p = 0 ("flat", lr = all), 6×{4,5,6,7} N=3,
    V = 0.1–1000 (15 each, same grid recipe in V), same metadata; filenames carry only
    `interaction_strength` (no `scaling`).
  - ED `cluster-data/exact-diag/cylinder/dd-ints/`: V=50, 6×{4,5,6} N=3 (15 each), 6x7 (11), 8x4 (22).
  - ED `cluster-data/exact-diag/cylinder/flat-ints/` (2026-10-15): flat Ui = 0, 6x4 N=3 and 8x4
    N=4, states + `ee_half` + `ee_cuts` in metadata; `full_basis` and `H` set to `nothing`.
  - Lx sweep (2026-10-14), all 7 p, in the same `obc/<folder>`: 8×{4,5,6,7} N=4 and 10x4 N=5, one
    file per (size, p), named after the smallest s (flat: V = 0.1), whose metadata holds the sweep:
    `all_gammas`, `all_columncorrs1/2`, `all_evencolumndist1/2`, `all_occs1/2` (25 points; no
    energies). 8×Ly 12–120 MB, 10x4 210 MB. Read the metadata only (`jldopen`, as the 8×Ly block).

## Open questions

- [ ] Does the exponent p change the physics? Hints: small p curves away from a power law in Ly,
  p ≥ 2 doesn't; b(p) − b(0) grows past p = 1. Needs more sizes → which family and method (Next).
- [ ] Why does the GS CDW signal at Γ → ∞ collapse with Ly (S(π)/N rise 0.40 → 0.035 for
  Ly = 4–7)? Same at N = 4 / Lx = 8, so it's Ly, not N. Is Γ_c(Ly) at fixed N meaningful there?
- [ ] GS-only: Hill and midpoint Γ_c disagree at Ly ≥ 6 (small signal). Which to use / show both?
- [ ] Γ_c converged in Lx at fixed Ly? Ly = 4: yes (Lx = 6–10). Ly ≥ 5: 8×Ly 10–25% above 6×Ly
  (only two Lx). Binder crossing in Lx with MPS/TTN?
- [ ] What is the pre-transition dip in S(π)/Binder? Known: it deepens with p and Ly, is absent
  for p ≲ 1, and comes with same-column pairs separating in y; in 6x4 there are also true GS
  level crossings near Γ → 0 and Γ_c.
- [ ] Is state 2 also excited (not a CDW partner) at 8×Ly, Ly ≥ 5? No energies in the sweep files.
- [ ] V = 50 TTN data at 12x6+ for the torus?
- [ ] Does fPEPS reach useful open sizes / accuracy (N conservation, D convergence of Γ_c)?
- [ ] C⁴ / PES vs s at TTN sizes (the paper's more robust diagnostics).
- [ ] Does the fitted c (vertical cuts, cylinder) converge with Lx, and to what, in the Laughlin
  and CDW regimes? 6x4 / 8x4 at Ui = 0 give 3.2 / 2.7, i.e. not converged.

## Gotchas

- `cylinder/flat-ints`: `full_basis` is `nothing` in the metadata, rebuild with
  `get_lattice_params_from_metadata(m)` (same basis order, checked against `ee_half`). 8x4
  `ee_cuts` has zeros for ℓ > Lx/2 (saved before mirroring; 6x4 is mirrored): mirror on reading.
- `compute_reduced_density_matrix` returns a dense ρ_A over all N_A ≤ N configs of A: keep
  ℓ ≤ Lx/2 (8x4 ℓ = 7 would be 24158², 9 GB).
- CurveFit is v1.4 (no `linear_fit`; the old TEE block at ~line 1853 fails): use LsqFit
  `curve_fit` or `design \ y`. CurveFit and LsqFit both export `stderror`.
- Manifold is 2-fold on the torus: densities from diagonalizing ⟨Ψi|n|Ψj⟩ (`occs_eig1`). On open
  y it isn't for Ly ≥ 5 (E2 − E1 ~ 0.03–0.09 at Γ → ∞): use the GS (`column_corr1`); averaging
  states 1, 2 (the old default) mixes in an excited state and inflates Γ_c at large Ly.
- Binder at fixed N has a disordered baseline 2/(3N) (0.222 for N=3, 0.167 for N=4); CDW plateau
  for perfect TT is 2/3 (S(π)/N → 1). With N ≤ 6 the offset distorts crossings. Growing Ly at fixed
  N doesn't lengthen the ordering direction, so there's no Binder crossing in Ly.
- The OBC CDW plateau is far from perfect TT and drops with Ly (S(π)/N 0.71, 0.53, 0.44, 0.40).
- OBC 6x6, 6x7 dip at s ≤ 0.15 and cylinder 8x4 is unconverged for s < 0.25 (|H| > 1e4,
  Lanczos); the analysis uses s ≥ 0.2 / 0.25. For p ≠ 3 the large-|H| limit moves (V/s^p).
- Open y: α = N / (ν Lx' (Ly − 1)) with Lx' = Lx − 1 on OBC, so 6×Ly N=3 files carry
  α = 1.2/(Ly−1) (0.4, 0.3, 0.24, 0.2); on OBC α drifts with Lx
  at fixed Ly (0.40, 0.38, 0.37 for Lx = 6, 8, 10 at Ly=4).
- `occs` layout: ED `[x,y]`, TTN `[y,x]`. TTN files without `"occs"` are skipped.
- ED torus 10x5 V=50 is one file with `all_magspacs` / `all_occs_eig1`, not one file per s.
- OBC 8x4 s = 0.1, 0.17 have no `column_corr*` in metadata.
- The "look at CDW transition for Gamma…" block computes `column_corr*` if missing and **writes
  them to the file metadata** with `modify_data` — ask first.
- Open y, Ly ≥ 5: states 1, 2 are not a degenerate CDW pair (E2 − E1 ~ 0.03–0.09 at Γ → ∞), so a
  2-state manifold fidelity picks up excited-state crossings; use GS-only fidelity there.
- `write_data_jld2` never overwrites: an existing name gets `-mk-N` appended (`prep_file`), and a
  missing target folder silently writes into `pwd()`. Write to a fresh name and `mv`; check the
  folder exists (don't create folders: the user manages them).
- `flat-ints` files store the defaults `magnetic_spacing = 1.0`, `exponent = 3.0` in metadata:
  for p = 0 use Γ = `interaction_strength`, and don't filter `find_data_file` on V.
- `"powerlaw"` filenames have a different key order than the dd ones (`ed-hopping_anisotropy-…`);
  `find_data_file` matches them fine.
- Hand-excluded TTN points: 14x7 s = 20.0, 12x6 s = 1.53.
- The 8x4 dd sweep sits in the old `magnetic_spacing-0.1` file and its `all_gammas` has 39 points:
  the 25 new ones plus the old coarse 8x4 points (incl. Γ = 10210, 50000, i.e. s = 0.17, 0.1) and
  3 duplicates. Γ_c is unaffected (17.2 vs 17.1 with Γ ≤ 1000).

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
- 2026-10-08 — `"powerlaw"` scaling (ED, tevo, TTN, control code); OBC 6×Ly N=3 data for
  p = 0.5–2.5 and p = 0 (flat-ints): Γ_c ∝ Ly^b, b ≈ 1.1 p + 1.5, all p diverge incl. p = 0.
  "Same physics?" check (user's doubt about the fits): same endpoints and crossover above Γ_c,
  p-dependent precursor below Γ_c for p ≳ 1.5, smooth GS path; fit-free Γ_c gives the same b(p).
  Found that open-y order parameters average the GS with an excited state for Ly ≥ 5.
- 2026-10-09 — GS-only order parameters (`if_gs_only` in the block) and p = 0 in the block. GS
  Γ_c much smaller at large Ly (p = 3, 6x7: 103 vs 200), b ≈ 0.6 p + 1.3. GS CDW signal collapses
  with Ly at fixed N (S(π)/N rise 0.40, 0.16, 0.064, 0.038) → Lx check at fixed Ly suggested next.
- 2026-10-09 — Reading path verified against buildHam (8x4–8x6 exact); GS signal holds with Lx at
  Ly = 4 (10x4: 0.365), collapses with Ly also at N = 4. Cluster runner block for the Lx sweep
  (8x4–8x7, 10x4 × 7 p) added to `execute-ed.jl`; one small file per (size, p) under the historic
  name of its first Γ point in the existing per-p folders. Not run yet.
- 2026-10-14 — Lx sweep back: 8×Ly N=4 and 10x4, all p, one file per (size, p) with `all_*`
  metadata. 8×Ly block: Γ_c(Ly) still grows, Ly = 4 Γ_c flat in Lx, GS signal collapses with Ly as
  at 6×Ly. Error bars on b; local exponents: single power law poor for small p, clean for p ≥ 2.
  User: larger sizes needed, idea fPEPS; size family / MPS-TTN vs fPEPS discussed (Next).
- 2026-10-15 — EE / central-charge side thread (`exact-diag/daily-things.jl`, end). Half-cut EE
  (x ≤ Lx/2), OBC 6x4 N=3 flat Ui: S = 1.749 (Ui = 0), 0.948 (Ui = 50). Periodic x / open y
  chosen for Calabrese–Cardy (cut scanned along the physical direction). Cut-size scan on the
  cylinder, Ui = 0: c = 3.2 ± 0.1 (6x4), 2.7 ± 0.07 (8x4). `compute_reduced_density_matrix`
  rewritten as ρ_A = Ψ Ψ† with sparse Ψ[a,b] (was an O(dim²) pair loop; same ρ_A to 1e-16).
  Next: finite-size scaling of c for Laughlin and CDW.
