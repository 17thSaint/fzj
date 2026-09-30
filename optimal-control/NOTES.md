# Optimal control & time evolution — technical notes

Long-form companion to [CONTEXT.md](CONTEXT.md): implementation details, measured numbers and
pitfalls by topic.

## 1. Infrastructure

**Python side — `quocs_common.py`** (2026-07-17 refactor, 2026-09-25 cleanup)
- Julia env bootstrap + `jl` handle, `JuliaFoM` base (`to_julia_dict`), `halfstep_bins`,
  `fourier_pulse` (`scaling_lambda`, `max_frequency`), `dcrab_algorithm_settings`,
  `run_optimization`, `load_best_controls`, `plot_best_pulses`, `finish_run`.
- Configs are import-safe (run code under `main()`), paths resolved from `__file__`.
- `finish_run` writes one figure per run, `local-figs/<name>_<run stamp>.png`, never overwritten
  (user wants a record of every run).
- `JULIAUP_CHANNEL=1.11` pin is conditional on 1.11 being installed (only 1.12.6 is now; the
  unconditional pin made juliaup hard-error and broke every config).
- magspacGrad subclasses magspacRamp; maggrad has its own numpy `retarded_response`.

**Julia side — setup/eval split**
- Every FoM is `setup_*(dict)` once (endpoint ED, constant matrices) + per-call
  `compute_fidelity_*(pulses, dict, setup)`; the FoM `__init__` holds the setup tuple
  (`self._setup`) via juliacall. At 4x4 N=2 the speed gain is minor for run_timeevo-based FoMs
  (per-half-step `buildHam` rebuild dominates, ~3 s/eval for intstren), but it makes the FoM
  deterministic (see txRamp).
- Shared helpers: `optimization_ed_params(pd, "key"=>v...)` (positional Pairs — Julia kwargs can't
  be strings), `groundstate_manifold_fidelity` = Σ|⟨c|t⟩|² / n, `magnetic_gradient_response(B,
  steps|h, w) -> (x, xdot)` (O(n)), `get_critical_dt_at_spacing`, `halfstep_times`,
  `zvd_impulses`. `make_tevo_params` precomputes `magnetic_gradient_integral_time`, so `timeham`
  no longer does O(n) work per half-step. `time_evolution` no longer `display`s t_evo_params at
  output_level 0 (it flooded every FoM eval).
- `ramptime`/`dt` in `compute_fidelity` are dict-overridable (were hardcoded).

**Verification method** (used for both cleanups; the frozen mirror has since been deleted)
- Frozen mirror of the repo (symlinks + copies of edited files) at a path whose *first* `fzj`
  substring is the real one (the scratchpad path contains "fzj" and breaks `include_other_files`),
  plus a harness that extracts a `tevo-daily-things.jl` block by its `###` header, forces
  `if true`/`if_all=true`, and runs it headless (`task_local_storage(:SOURCE_PATH)` so relative
  includes resolve). All replay/ZVD/posicast blocks reproduced baseline to ~1e-15.
- Reference results at the time: intstren FoM bitwise identical; pinned within 5e-12 (Lanczos
  noise in the degenerate target manifold; manifold fidelity is subspace-invariant).

## 2. QuOCS pitfalls

- **Never use `get_best_controls()`.** It rebuilds the pulse from the *last* super-iteration's
  random Fourier basis; if a super-iteration was cut off (e.g. `max_eval_total` reached),
  `best_xx` belongs to an earlier basis and the pulse is wrong. Read the `*_best_controls.npz`
  BestDump in `optimization_obj.results_path` (`load_best_controls`); it reproduces the logged
  best FoM exactly. Likely affects dCRAB runs terminated mid-SI too.
- **AD mode:** `get_FoM` receives the pulses as rows of a **complex64** 2D jnp array (tracers
  under `jax.grad`) — take `jnp.real(...)` and stay in jnp. The FoM must be pure JAX (cannot call
  Julia).
- **AD settings layout:** L-BFGS-B stopping settings (`ftol`, `gtol`, `maxls`, `max_eval_total`)
  go in `algorithm_settings["stopping_criteria"]`, not `dsm_settings` (dCRAB-only); the global stop
  reads `algorithm_settings["max_eval_total"]` — set both.
- `shrink_ampl_lim` must stay off in AD mode (value-dependent Python branches, not traceable);
  default `limit_pulse` clipping (jnp.maximum/minimum) is fine.
- **dCRAB endpoint pinning:** `fourier_pulse`'s default parabolic `scaling_function` pins the update
  at *both* ends, freezing e.g. B(T) at the guess. Pass `scaling_lambda` (e.g. `"lambda t: t /
  t[-1]"`) when an end must be free. Measured on maggrad (cost, lower better; ceiling 3.63):
  both-ends-pinned 2.40–3.00, start-only 2.26–2.32, unpinned 2.26–3.17; with the end free the guess
  barely matters (<0.1%).
- **Bang-bang outputs** are `limit_pulse` clipping Fourier coefficients that Nelder–Mead drove to
  ~1e6–1e7, not a bad guess. dCRAB NM coefficients are unbounded (`amplitude_variation` only sizes
  the simplex). A bang-bang *guess* would be a real problem (every update clipped back off).
- jax 0.11.0 installed into `quocs-env` via `uv pip install --python quocs-env/bin/python jax`.

## 3. Time-evolution discretization

- **Half-step grid.** Pulses are sampled on the RK4 dt/2 grid that `timeham` indexes. A switch in
  a piecewise pulse must land on it: snap dt via `nsw = ceil(tsw/(dt_crit/2)); dt = 2*tsw/nsw`.
- **Snap + average only work as a pair.** Snapping alone makes ringing *worse* (a sample lands
  exactly where the trapezoid error is largest and `t < tsw` gives it the full B2). Also set the
  single sample on the switch to `(B1+B2)/2` — exact for the trapezoid rule. Measured residual
  ringing (peak-to-peak of a after the switch): neither fix 6e-5 (luck), snap only 3e-4, both 2e-15.
- **Exception at t=0:** the first sample takes the full amplitude (no "before" side inside the
  domain); halving it silently loses half the first interval.
- **Power-of-two half-step counts:** T/(dt/2) = 256 gives exactly 256.0; 252 gives
  252.00000000000003 and the `ceil` adds a stray sample.
- **`isapprox` needs explicit `atol`** on trapezoid-derived spacings: default rtol ~1.5e-8 is
  tighter than the ~1.6e-7 O(h²) quadrature error (cost a 22-min rerun). Exact quantities
  (`spacings[1] == a0`, final U after `long_range_scaling`'s 5-digit rounding) compare fine.
- **dt_crit is from the t=0 Hamiltonian only**; `run_timeevo` checks nothing later. Any pulse with
  a(t) < a_start (U ~ 1/a³) needs dt sized from min(a) — otherwise silent divergence hidden by the
  per-step renormalization. At a=0.1, intstren 10: U1 = 10000, dt_crit ~ 2e-4.
- **Indexing:** `tevo_wavefunc` has `1+(nsteps+1)/2` columns (last never written → `[:,end-1]`);
  `instant_spec` has `1+(nsteps-1)/2`, all written, at t = k·dt for k=1..nsteps_rk4 — the
  instantaneous series starts one step in, not at t=0.
- **Runtime:** `if_instant_gs=true` forces a dense diagonalization every RK4 step (~30 min at 6284
  steps, 120×120). Setting it false is the big lever; shorter tmax the mild one.
- **Degenerate endpoints:** dd endpoints are exactly 2-fold degenerate and `run_normal_ed`'s Lanczos
  sometimes returns only one of the pair (random start) → target manifold = gs + excited,
  fidelity ~0.45. Diagonalize densely and assert an isolated speccount-fold level.
  `setup_intstren_ramp` still uses Lanczos.

## 4. Problems

### txRamp (`config_txRamp.py`, control-functions.jl) — broken
- FoM ill-defined: at tx=0.001 the gs is quasi-degenerate (decoupled columns), so each Lanczos draw
  gives a different state; the old per-eval-ED FoM gave 0.19–0.49 for the *same pulse*. Old records
  (0.7357, 2026-07-09) partly reflect lucky draws. With cached setup the FoM is reproducible but
  the ceiling depends on the draw (~0.09 in one smoke run). Fix needs a physics decision:
  degenerate-manifold fidelity or a well-defined start state.
- dt=0.05 > dt_crit 0.0067 → run_timeevo's guard rejects it; the FoM cannot run at all.
- Fixed: `time_txRamp.initial_value` was 1.0 while the ramp is 2.0 (plot axis only); magic 81 bins
  now derived.

### pinnedRamp (`config_pinnedRamp.py`)
- Two sequential hopping ramps (ty then tx), corner-pinned product state → FCI manifold. Run
  20260709_162520. Dormant.

### intstrenRamp, dCRAB and AD (`config_intstrenRamp{,_AD}.py`)
- 4x4 N=2, ramptime 1.0, dt 0.005, strong ULR → FCI. Linear ~0.889, dCRAB 0.9101 (2026-07-14),
  AD 0.9059 (200 evals, first run) → **0.9158** (154 evals, 88 s); eval-budget-limited.
- AD architecture (`intstren-ramp-ad-functions.jl`: `setup_intstren_ramp_ad`): Julia computes
  pulse-independent constants once — endpoint manifolds and H(u) = H_hop + u·H_int from the saved
  undressed matrices (`getHopping`/`getInteraction`, no `buildHam`). H_int is dressed at
  `intstren_max` (12.0) and divided back out, **not** at unit strength, else couplings below
  `interaction_cutoff` (1e-5) at u=1 but above it at large u would be missing. Checked via endpoint
  eigenstate residuals (~1e-14 flat; ~1e-6 non-flat from 5-digit rounding of U; tolerance
  1e-5·N(N−1)). Pinning/disorder/periodic potential rejected (don't scale linearly with u).
- Non-flat profiles (exp/gaussian/rydberg/dd) supported — every `long_range_scaling` branch is
  u·profile; both FoMs forward `scaling_type`/`corr_length`/`sigma`/`blockade_radius`/
  `magnetic_spacing`. Validated exp (1.5) and rydberg (1.2): FoM diff ~1e-8, grads ~1e-11.
- JAX propagation: RK4 half-step grid mirroring `runge_kutta_step` incl. per-step renormalization,
  `lax.scan`, fidelity 0.5·tr(F†F); interaction applied elementwise (diagonal); hopping dense below
  1024 dims, `BCOO` sparse above (both verified identical). FoM matches Julia to 1e-16 (linear) /
  1e-9 (wiggly); grad vs finite differences ~1e-11. ~0.03 s/eval vs 3–10 s in Julia.

### maggradPulse (`config_maggradPulse.py`, maggrad-pulse-control-functions.jl) — classical
- Shapes B(t) so a(t) = a0 + (1/w)∫B(t') sin(w(t−t')) dt' sits at a_target over the end window.
  No ED: Julia only supplies the RK4 dt (pulse must sit on the dt/2 grid). 4x4 N=2, a0=0.1,
  tmax 2.5 → dt 2e-4, 25001 samples; ~1 ms/eval.
- O(n) trick: sin(w(t−t')) = sin wt cos wt' − cos wt sin wt' turns the retarded integral into two
  cumulative trapezoids (vs O(n²)); agrees to ~4e-19; self-checked by
  `check_response_implementation`.
- **Two caps:** resonant pumping caps the window *mean* at ~b_max·T/(2w); holding steady needs a
  standing gradient B = w²(a−a0), so the largest settleable spacing is a0 + b_max/w². The second
  one binds. Below it, mean-on-target and no-ringing are irreconcilable (solution jumps between
  "on target, ringing" and "steady at cap"); above it (b_max ≥ w²(a_target−a0) = 190) an
  unweighted sum settles.

      b_max  cost               mean     std
         30  |Δmean| + std      2.0000   0.4391   (cap 0.40)
         30  |Δmean| + 5·std    0.4000   0.0000   abandons target
        200  |Δmean| only       2.0000   6.9648   std term is load-bearing
        200  |Δmean| + std      2.0000   0.0013   100k evals

- Settings (user's choices, 2026-09-22): b_max 200, w 10, a0 0.1, a_target 2.0, tmax = drive_time
  2.5, window 0.25, cost |Δmean| + std unweighted (variance_weight knob removed at user's request),
  `min_spacing_penalty` 50 on the time-averaged violation below a0.
- **Zero crossings are set by the guess** (U ~ 1/a³ diverges at a=0; window terms only constrain
  the end). Resonant guess: min(a) −3.99, 6 crossings; smooth ramp to the holding value: min(a) =
  a0, 0 crossings, std 1e-4 with no floor penalty. The floor penalty also removes crossings but
  costs ~10× in std. Final: mean 2.0000, std 0.0000, min(a) = 0.1.
- Analytic optimum (→ posicast): B = w²Δa/2 for π/w, then w²Δa; settles by t=0.314, std ~1e-5.
- Bug fixed: `build_optimization_dictionary` read module-level `drive_time`; now a parameter.

### Posicast / ZVD gradient blocks (`tevo-daily-things.jl`)
- Posicast (2-step): B1 = w²Δa/2 on [0, π/w), B2 = w²Δa after; w=10, a0=0.1→2.0: B1 95, B2 190,
  tsw 0.31416, run to 4·tsw. With snap + average: ringing 2e-15, deviation from the closed form
  a0 + (Δa/2)(1−cos wt) 1.6e-7 (pure O(h²)), min(a) = a0 so the t=0 dt_crit holds throughout.
- ZVD (3-step), written over the impulse train: amps [1/4, 1/2, 1/4] at [0,1,2]·π/w convolve to the
  staircase Bf·[1/4, 3/4, 1]; feeding [1/2, 1/2] reproduces the posicast [95, 142.5, 190]. Closed
  form a = a0 + Δa Σ A_i(1−cos w(t−t_i)) over fired impulses. One dt-snap covers all switches.
  tmax 3π/w, 4713 steps: ringing 2.4e-15, closed-form deviation 1.6e-7, energy drift after the
  last impulse exactly 0 (H genuinely static).
- Manifold population: **ZVD 0.876 vs posicast 0.857** (the extra half-period is a gentler quench).
- Caveat: ZVD's O(ε²) frequency robustness is asserted algebraically (Σ A_i e^{iwt_i} and its
  w-derivative), never demonstrated with a detuned oscillator.
- Verification idea worth reusing: check the O(n) integral against `get_magnetic_gradient_integral`
  itself at a few samples (dict with `magnetic_gradient_time`, `when_dt_ends` past the end, `dt`,
  `trap_frequency`), not against a second copy of the same arithmetic.

### magspacRamp — control a(t) (`config_magspacRamp.py`, magspac-ramp-control-functions.jl)
- Setup (user's choices, 2026-09-24): 4x4 N=2 pbc, dd intstren 10, a 0.5 → 2.0 (U1 80 → 1.25),
  T 1.0, dt 0.005, target = dd gs manifold at a_end, speccount 2. Baselines: quench 0.8963,
  linear 0.9253, ZVD (w=10, Bf=150, settles at 2π/w = 0.628 then holds) 0.9111 — ZVD front-loads
  the transfer and the static hold buys nothing; settling time is set by w, not T.
- FoM avoids run_timeevo: H = H0 + Σ_k U_k(a) diag(M_k) (M_k diagonal), adaptive Strang split,
  m = ceil(h·spread(D)/θ), θ=2 → ~10 ms/eval (vs ~2 s exact-exp); agrees with run_timeevo ~2.5e-6.
  Gotcha: applyHam's range is count(|U|>cutoff)−1, so extracting M_k with unit vectors e_k gives
  zero — perturb all-ones instead.
- Scaling 16 t²(1−t)² pins value and slope at both ends (a leaves/arrives at rest, so B is finite).
  Bounds a ≥ 0.1.
- Result (run 20260924_150154): **0.9392**, 4485 evals in 50 s, RK4 replay 0.93916 — but a(t) swings
  0.4–3.8 and read-back |B| peaks 3302 (ZVD 150). Low-passing (≤2 cycles/T, a ∈ [0.5, 2.4]) still
  gave |B| 8231 (F 0.9403): clipping a(t) makes slope kinks = delta kicks in B. → control B instead.

### magspacGrad — control B(t) (`config_magspacGrad.py`)
- Julia: `magspac_gradient_fom`, `magspac_spacing_from_gradient`, `zvd_gradient_pulse`. Hard QuOCS
  limit |B| ≤ 300 (2× ZVD). Update scaled by (1 − t/T), pinning B(T) at the holding value 150
  (user's choice). FoM = F − R, R = sqrt((a(T)−a_end)² + (adot(T)/w)²) (ringing amplitude after T),
  unweighted. T = 2π/w = 0.628 (user: ramp time must equal the ZVD pulse's own duration),
  dt = 2T/256.
- The pinned last sample makes the trapezoid ramp the final ZVD jump over one half-step: guess
  R 4.6e-3 (vs 7.5e-5); optimizer removes it.
- Run 20260925_102851: **F 0.9204, R 1.9e-6**, |B| touches 300, min a 0.5; run_timeevo via scaling
  "magnetic_gradient" agrees to 6e-7. T = 1.0 run (20260925_102017): 0.9221.

### magspacGradCRAB — plain CRAB for B(t) (`config_magspacGradCRAB.py`)
- Advised setup (2026-09-30): regular CRAB = dCRAB with `super_iteration_number=1`; 4 frequencies
  from `Uniform` [0.5, 4.5] (stratified: one in each [k−½, k+½] cycles/T); Nelder–Mead; zero
  guess dropped (user) for a smooth guess B = 150(1−cos πt/T)/2 so both ends sit at the holding
  gradients (B(0)=0 holds a_start, B(T)=w²Δa holds a_end), default parabolic envelope pins them.
  Everything else = magspacGrad (T = 2π/w, 256 half-steps, |B| ≤ 300, FoM = F − R).
  Optional seed argument (`python config_magspacGradCRAB.py <seed>`, folder `_seed<k>`).
- Guess: F 0.9048, R 0.50. Seeds 1–4 (runs 20260930_14344[6-7]): F 0.9048 / **0.9052** / 0.9049 /
  0.9002, R ~2.5e-6, each fatol-converged after ~265 evals (budget 5000); coefficients O(10–100),
  |B| ≤ 155, a(t) monotone. run_timeevo replay of seed 2: 0.905219 (fast 0.905222).
- Reading: NM removes the ringing and stops in the guess's basin; fidelity is set by the guess's
  shape, not found. Below ZVD 0.9111 and ZVD-guess dCRAB 0.9204.
- dCRAB option: `python config_magspacGradCRAB.py <seed> <super_iterations>`; SI > 1 sets
  `max_eval` 400 per SI (Nelder–Mead `stopping_criteria`) and names runs
  `magspacGradSmooth_dCRAB_seed<k>` (not matched by the CRAB block's loader). Runs 20260930_150652,
  10 SI: best FoM per SI (seed 1–4) 0.9194 / **0.9208** / 0.9186 / 0.9197, R 2.6–5.7e-6, 2500–3100
  evals, ~27 s; seeds 3, 4 touch |B| = 300. Pairwise RMS pulse difference 82–180 at equal F: many
  different pulses reach ~0.92, so no unique optimal shape at T = 2π/w.
- Best (seed 2) replayed by the last block of `tevo-daily-things.jl` (loads only the dCRAB runs,
  picks max FoM): run_timeevo 0.920779 vs fast 0.920780. a(t) is non-monotonic: 0.5 → 1.15 at
  t ≈ 0.18, back to 0.55 at t ≈ 0.32, then → 2.0 (overshoot 2.03); |B| ≤ 211, dips to −70.
  Energy plot with `nev` = 3 (only sets the instantaneous states; 2 manifold states are evolved):
  E1/E2 degenerate throughout; transported energy leaves them by t ≈ 0.05 and ends ~0.11 above
  (E3 ~0.28 above).
- Code B is δB from the physical holding gradient w²a0 (trap rest point a = 0; user, 2026-09-30),
  so the |B| ≤ 300 cap is −250…350 physical. Same for all magspac/maggrad configs.

### Larger systems (planned, 2026-09-30)
- Motivation: at 4x4 N=2 (dim C(16,2) = 120) a quench gives ~0.90 and control adds ~0.02; quench
  fidelity is expected to fall roughly exponentially with N (unverified), leaving room for control.
- 8x4 N=4: C(32,4) = 35,960; 10x5 N=5: C(50,5) ≈ 2.1M.
- `setup_magspac_ramp` / `evolve_magspac_ramp` are dense (H, exp(−iH0 d), `eigen`): ~20 GB per
  matrix at 8x4, so that route is out. Options: `run_timeevo` (sparse RK4) as the FoM, or a sparse
  Krylov version of the fast FoM. Endpoint manifolds need Lanczos with a check on the degenerate
  pair (§3 pitfall).
