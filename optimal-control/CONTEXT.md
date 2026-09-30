# Optimal control & time evolution — project context

Living status file. Update the **Current state** and **Log** sections whenever something changes.
Last updated: 2026-09-30

## Goal

State preparation in the synthetic-dimension / FQH lattice (small ED systems, mostly 4x4 N=2 pbc):
find time-dependent control pulses that carry an easily prepared groundstate manifold into a
target manifold (FCI / weak-dd state) with high fidelity. Pulse shapes are optimized
with QuOCS (dCRAB, and AD for one problem); the time evolution itself is exact (ED + RK4 or
split-operator).

Motivation (paper draft Sec. V + App. F, see `CLAUDE.md`): the Ui path from the trivial CDW to the
Laughlin FCI stays gapped, so it is a candidate finite-time preparation route; App. F finds that
higher excited states, not the lowest gap, dominate the adiabatic time. Experimentally Ui (here the
dd strength via the spacing a) is tuned by a Stern–Gerlach gradient B separating the synthetic levels.

Current physical focus: **controlling the dipolar spacing `a(t)` through a magnetic field
gradient `B(t)`** acting on a trapped mode (trap frequency `w`):

    a(t) = a0 + (1/w) int_0^t B(t') sin(w(t - t')) dt'        U_r = intstren / (r a)^3

## Layout

| Where | What |
|---|---|
| `optimal-control/quocs_common.py` | Shared: Julia bootstrap (`jl`), `JuliaFoM` base, `fourier_pulse` (`scaling_lambda`, `max_frequency`), `halfstep_bins`, `dcrab_algorithm_settings`, `run_optimization`, `load_best_controls` (reads BestDump npz), `finish_run` (per-run figure) |
| `optimal-control/config_*.py` | One optimization problem each (see table below); run code under `main()` |
| `optimal-control/QuOCS_Results/<stamp>_<name>_<alg>/` | Run outputs (gitignored); `*_best_controls.npz` is the source of truth for the best pulse |
| `optimal-control/local-figs/<name>_<stamp>.png` | One figure per run, never overwritten |
| `optimal-control/quocs-env/` | Python venv (quocslib, juliacall, jax 0.11) |
| `exact-diag/time-evolution.jl` | RK4 time evolution (`run_timeevo`, `timeham`, `make_tevo_params`), `long_range_scaling` incl. `"magnetic_gradient"`, `magnetic_gradient_response` (O(n)), `get_critical_dt(_at_spacing)`, `halfstep_times`, `zvd_impulses` |
| `exact-diag/control-functions.jl` | `optimization_ed_params`, `groundstate_manifold_fidelity`, txRamp setup/FoM |
| `exact-diag/*-control-functions.jl` | Per-problem Julia setup + FoM (intstren, intstren-AD, pinned, maggrad, magspac) |
| `exact-diag/tevo-daily-things.jl` | Scratch blocks (`### header` + `if true/false`): baselines, analytic pulses, replays of optimized pulses |
| `exact-diag/tests-tevo.jl` | Time-evolution tests |

## Optimization problems

| Config | Control | Transport | Best result | Status |
|---|---|---|---|---|
| `config_txRamp.py` | tx(t) | tx~0 gs → isotropic gs | (0.7357, not meaningful) | **Broken**: dt=0.05 > dt_crit 0.0067 (guard rejects); FoM ill-defined (quasi-degenerate start, Lanczos-draw dependent) |
| `config_pinnedRamp.py` | ty(t), then tx(t) | corner-pinned product state → FCI manifold | see run 20260709_162520 | Works, dormant |
| `config_intstrenRamp.py` | U(t) (dCRAB) | strong ULR → FCI, T=1, dt=0.005 | 0.9101 (linear 0.889) | Works, dormant |
| `config_intstrenRamp_AD.py` | U(t) (JAX AD, L-BFGS-B) | same | **0.9158** (154 evals) | Works, dormant; eval-budget-limited |
| `config_maggradPulse.py` | B(t), classical FoM | spacing a0=0.1 → a_target=2.0 held steady | mean 2.0000, std ~0 | Done; analytic posicast/ZVD solution is better anyway |
| `config_magspacRamp.py` | a(t) directly | dd gs a=0.5 → a=2.0, T=1 | 0.9392 | Done, but **unphysical**: read-back \|B\| ~3300 |
| `config_magspacGrad.py` | B(t), \|B\|≤300 | dd gs a=0.5 → a=2.0, T=2π/w=0.628 | **0.9204**, R=1.9e-6 | Current latest (run 20260925_102851) |

### Magnetic-spacing transport reference numbers (4x4 N=2 pbc, dd intstren 10, a 0.5 → 2.0)

| Pulse | T | Fidelity |
|---|---|---|
| Sudden quench | — | 0.8963 |
| ZVD gradient staircase (w=10, Bf=150) | 1.0 | 0.9111 |
| Linear a(t) (needs delta kicks in B → unphysical) | 1.0 | 0.9253 |
| dCRAB on B(t), \|B\|≤300 | 1.0 | 0.9221 |
| dCRAB on B(t), \|B\|≤300 | 2π/w | 0.9204 |
| dCRAB on a(t), unconstrained | 1.0 | 0.9392 (unphysical shaking) |

Conclusion so far (commit 71432d0, "optimizing posicast is bad"): with a physical gradient cap,
optimizing from the posicast/ZVD guess gains only ~0.01 over ZVD and does not beat the (unphysical)
linear ramp. Most of the fidelity budget is set by T and w, not pulse shape.

## Current state (2026-09-30)

- **Posicast (2-step) and ZVD (3-step) gradient pulses are the preferred pulses.** They move the
  spacing a0 → a_target cleanly (no residual ringing, a(t) never below a0, only one dt-snap needed).
- **Problem: in the quantum time evolution they don't follow the groundstate well.** They are
  designed to settle the classical trapped mode, not for adiabaticity: the transfer is done in
  ~π/w–2π/w and front-loaded, so the state is left behind (ZVD 0.9111 < linear a(t) 0.9253 at T=1;
  manifold population 0.876 ZVD / 0.857 posicast in the 0.1 → 2.0 maggrad blocks).
- **Next idea: use posicast/ZVD as the initial guess for QuOCS**, keeping their clean
  (non-ringing, physical-B) endpoints while shaping for groundstate following.
- Live blocks in `tevo-daily-things.jl`: magspac baseline (quench/linear/ZVD), a(t)-dCRAB replay,
  B(t)-dCRAB replay (last block).
- Last activity on this project: 2026-09-25 (library cleanup + magspacGrad run). Since then work
  moved to the dd critical spacing analysis in `synth-dims/`.
- No optimization currently running.

## Open questions / next steps

- [ ] **Posicast/ZVD as QuOCS starting point.** Already tried once each: `config_magspacRamp.py`
      (ZVD a(t) guess → 0.9392, but read-back |B| ~3300) and `config_magspacGrad.py` (ZVD B(t) guess,
      |B|≤300, T=2π/w → 0.9204). Gain over ZVD was only ~0.01. Knobs not yet explored: posicast
      (2-step) guess instead of ZVD; T longer than the pulse itself (ZVD stretched, or ZVD + shaped
      tail); lower w (slower mode, w = 2π/T); larger |B| cap; AD instead of dCRAB; FoM that
      rewards following the instantaneous groundstate along the path rather than only the end.
- [ ] Is the fidelity ceiling at fixed T physics (adiabaticity at this gap) or control? Sweep T
      (and w = 2π/T for a fair ZVD comparison) for the B-controlled problem.
- [ ] Larger gradient cap / longer T trade-off for magspacGrad.
- [ ] ZVD frequency robustness is only asserted algebraically — a detuned-w sweep (~12 classical
      lines, no ED) would demonstrate it.
- [ ] txRamp: fix dt, and decide physics of the FoM (degenerate-manifold fidelity or a
      well-defined start state).
- [ ] AD variant for the magspac problem? (FoM would need to be pure JAX.)

## Conventions & gotchas

- Pulses live on the RK4 **half-step grid** (dt/2) that `timeham` indexes. Switch times must land
  on it: snap dt so switch times are multiples of dt/2, **and** give the sample sitting on an
  interior switch the average `(B1+B2)/2` (not at t=0). Either fix alone is worse than neither.
- Choose half-step counts as powers of two (e.g. 256) to avoid `T/(dt/2)` = 252.00000000000003.
- dt_crit is computed from the t=0 Hamiltonian; `run_timeevo` only checks that. Any pulse with
  a(t) < a_start (U ~ 1/a^3) needs dt sized from min(a).
- dd endpoints are exactly 2-fold degenerate: Lanczos sometimes drops the partner → use dense
  `eigen` and assert an isolated speccount-fold level.
- Fidelity is the manifold fidelity `sum |<c|t>|^2 / n`.
- QuOCS: never use `get_best_controls()` (wrong basis if a super-iteration was cut) — read the
  BestDump npz. In AD mode the FoM gets complex64 jnp rows, must stay pure JAX; stopping settings
  go in `algorithm_settings["stopping_criteria"]` and `max_eval_total`.
- `fourier_pulse`'s default parabolic scaling pins both ends; pass `scaling_lambda` when an end
  must be free (e.g. B(T) of a holding gradient).
- `isapprox` on trapezoid-derived spacings needs an explicit `atol` (O(h^2) error ~1.6e-7 >
  default rtol).
- `if_instant_gs=true` in time evolution forces a dense diag every RK4 step — the dominant runtime.
- Verification harness: frozen mirror at `/home/patrick/fzj/.cleanup-mirror/main-git` + block
  extractor that runs a `tevo-daily-things.jl` block headless by its `###` header.

## Log

- 2026-04-15/20 — QuOCS connected to Julia ED; first txRamp dCRAB runs.
- 2026-07-09/10 — pinnedRamp (Goldman-style preparation) optimized.
- 2026-07-14 — intstren ULR → FCI linear vs dCRAB (0.889 → 0.9101).
- 2026-07-17 — AD (JAX) variant, 0.9158; refactor onto `quocs_common.py`; setup/eval split.
- 2026-09-04/08 — time-dependent magnetic gradient added to time evolution (trap frequency 09-21).
- 2026-09-22 — maggradPulse: classical dCRAB on B(t); steady-state cap `a0 + b_max/w^2`.
- 2026-09-23/24 — analytic posicast (2-step) and ZVD (3-step) gradient pulses, ringing 2e-15.
- 2026-09-24 — magspac a(t) ramp baseline + dCRAB (0.9392, unphysical |B|).
- 2026-09-25 — magspacGrad: B(t) control with |B|≤300, 0.9204; second library cleanup.
- 2026-09-30 — Status: posicast/ZVD preferred but poor groundstate following; next step is using
  them as QuOCS initial guesses. Context files + root `CLAUDE.md` created.
