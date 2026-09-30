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

**Convention (user, 2026-09-30): B(0) = B_hold unless the user says otherwise.** No step in the
gradient at t = 0: before t = 0 the lattice is held at rest at a(0). B in code and results is δB, the
deviation from the physical holding gradient w²a0 (= 50), so B_hold = 0 at a0 and physical
B = 50 + δB. Posicast/ZVD and magspacGrad (B(0) free) break the rule; old benchmarks only.

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
| `magspacGrad` | B(t), \|B\| ≤ 300: same | **0.9204** (T=2π/w), 0.9221 (T=1) | old benchmark (breaks B(0) rule) |
| `magspacGradCRAB` | same, plain CRAB, 4 freqs, NM, smooth guess | 0.9052 (4 seeds 0.900–0.905) | stuck near guess |
| `magspacGradCRAB` + 10 SI | same, dCRAB 10 SI × 4 freqs, B(0)=0 | **0.9208** (4 seeds 0.9186–0.9208) | latest |

Magspac baselines (4x4 N=2, dd intstren 10, T=1): quench 0.8963, ZVD 0.9111, linear a(t) 0.9253
(unphysical: needs delta kicks in B).

## Current state (2026-09-30)

- Smooth-guess pulses only (B(0) rule). Plain CRAB (4 freqs, Nelder–Mead) stays at the guess,
  F ≈ 0.905; dCRAB 10 SI reaches 0.9186–0.9208 over 4 seeds, same as the old ZVD-guess dCRAB
  (0.9204). Replay: last block of `tevo-daily-things.jl`.
- **4x4 N=2 is probably too small to show what control can do**: a quench already gives ~0.90 and
  optimization adds only +0.02. The landscape is flat (pulses 80–180 apart in RMS B reach the
  same F), so shapes are neither unique nor interpretable, and hard to implement.
- **Next: larger system.** Cheap check first at 8x4 N=4: quench and smooth guess with the sparse
  `run_timeevo`, timing one evolution. The fast FoM is dense and cannot reach 36k states.

## Open questions

- [ ] 8x4 N=4: does the quench fidelity drop enough to give control room? If so, which FoM:
      `run_timeevo`, or a sparse Krylov version of the fast FoM?
- [ ] Implementable pulses: band limit (≤ 2 cycles/T) or a penalty on dB/dt.
- [ ] Untried knobs for smooth-guess dCRAB: longer T, lower w, larger |B| cap, AD, FoM rewarding
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
- 2026-09-30 — plain CRAB for B(t) 0.9052; dCRAB 10 SI 0.9208 (4 seeds); B(0) = B_hold rule; 4x4
  judged too small → plan 8x4 N=4.
