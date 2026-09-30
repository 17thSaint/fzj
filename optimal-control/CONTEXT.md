# Optimal control & time evolution — status

Update **Current state**, **Open questions** and **Log** as things change. Technical detail:
[NOTES.md](NOTES.md). Last updated: 2026-09-30

## Goal

Prepare target groundstate manifolds (FCI / weak-dd) from easy ones with QuOCS-optimized pulses
(dCRAB; AD for one problem), exact ED time evolution, mostly 4x4 N=2 pbc. Motivation: the paper's
gapped Ui path (CDW → FCI) as a finite-time preparation route; App. F: higher excited states, not
the lowest gap, set the adiabatic time.

Current focus: dipolar spacing a(t) driven by a gradient B(t) on a trapped mode (frequency w),

    a(t) = a0 + (1/w) ∫_0^t B(t') sin(w(t - t')) dt',        U_r = intstren / (r a)^3

## Layout

- `optimal-control/`: `quocs_common.py` (shared helpers), `config_*.py` (one problem each),
  `QuOCS_Results/` (runs; best pulse = `*_best_controls.npz`), `local-figs/` (one figure per run).
- `exact-diag/`: `time-evolution.jl` (RK4, `magnetic_gradient` scaling), `control-functions.jl` +
  `*-control-functions.jl` (per-problem setup/FoM), `tevo-daily-things.jl` (baselines, analytic
  pulses, replays), `tests-tevo.jl`.

## Problems

| Config | Control → transport | Best | Status |
|---|---|---|---|
| `txRamp` | tx(t): tx≈0 gs → isotropic | — | **Broken**: dt > dt_crit; FoM ill-defined |
| `pinnedRamp` | ty then tx: pinned product → FCI | run 20260709_162520 | dormant |
| `intstrenRamp` / `_AD` | U(t): strong ULR → FCI, T=1 | 0.9101 / **0.9158** (linear 0.889) | dormant |
| `maggradPulse` | B(t), classical: a 0.1 → 2.0 held | mean 2.0, std ~0 | done; posicast/ZVD better |
| `magspacRamp` | a(t): dd a 0.5 → 2.0 | 0.9392 | unphysical (\|B\| ~3300) |
| `magspacGrad` | B(t), \|B\| ≤ 300: same | **0.9204** (T=2π/w), 0.9221 (T=1) | latest |

Magspac baselines (4x4 N=2, dd intstren 10, T=1): quench 0.8963, ZVD 0.9111, linear a(t) 0.9253
(unphysical: needs delta kicks in B).

## Current state (2026-09-30)

- **Posicast (2-step) and ZVD (3-step) gradient pulses are the preferred pulses**: clean transfer,
  no ringing, a(t) ≥ a0, one dt-snap.
- **But they don't follow the groundstate**: built to settle the classical mode in π/w–2π/w, not for
  adiabaticity (ZVD 0.9111 < linear 0.9253; manifold population ZVD 0.876, posicast 0.857).
- **Next: use them as QuOCS initial guesses.** Tried once (commit 71432d0, "optimizing posicast is
  bad"): ZVD guess gained only ~0.01 under the |B| cap.
- Last active 2026-09-25; work since moved to `synth-dims/`. Nothing running.

## Open questions

- [ ] Posicast/ZVD as guess — untried knobs: posicast instead of ZVD; T longer than the pulse (or
      pulse + shaped tail); lower w (w = 2π/T); larger |B| cap; AD instead of dCRAB; FoM rewarding
      groundstate following along the path, not just at the end.
- [ ] Is the fidelity ceiling at fixed T physics (adiabaticity) or control? Sweep T.
- [ ] Demonstrate ZVD frequency robustness (detuned-w sweep, classical, ~12 lines).
- [ ] txRamp: fix dt; decide FoM physics (degenerate-manifold fidelity or defined start state).

## Log

- 2026-04 — QuOCS connected to Julia ED; first txRamp runs.
- 2026-07-09 — pinnedRamp optimized.
- 2026-07-14/17 — intstren dCRAB 0.9101, AD 0.9158; refactor onto `quocs_common.py`.
- 2026-09-04/21 — time-dependent magnetic gradient + trap frequency in time evolution.
- 2026-09-22 — maggradPulse (classical); steady-state cap a0 + b_max/w².
- 2026-09-23/24 — analytic posicast and ZVD pulses (ringing 2e-15).
- 2026-09-24/25 — magspac a(t) dCRAB 0.9392 (unphysical); magspacGrad 0.9204; library cleanup.
- 2026-09-30 — status recorded; next step posicast/ZVD as QuOCS guesses.
