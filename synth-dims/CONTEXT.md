# Dipolar phase transition — status

Update **Current state**, **Open questions** and **Log** as things change. Last updated: 2026-10-27

## Goal

The paper leaves open *"the degree of nonlocality required"* to drive the FCI → CDW transition.
Replace infinite-range Ui by a finite-range interaction along the synthetic direction,

    U(r) = V / (r s)^p      (p = 3: scaling "dd", s = magnetic_spacing, V = intstren = 300 or 50)

and find the critical point (small s → CDW, large s → Laughlin) and whether it survives the
thermodynamic limit. Scaling variable Γ = V / s^p (nearest-neighbour synthetic coupling). Fixed: ν = 1/2.
Torus: ρ1D = 1/2, Lx × Lx/2, α = 1/Ly; sizes 6x3, 8x4, 10x5 (ED); 12x6, 14x7, 16x8 (TTN, H100).
Open y (OBC / cylinder, ED), V = 50: 6×Ly N=3 and 8×Ly N=4 with Ly = 4–7, 10x4 N=5.
Other exponents: scaling `"powerlaw"` with key `"exponent"` = p (p = 3 equals `"dd"`).

## Current state (2026-10-20)

**Edge / central-charge side thread (2026-10-15/16, ED, flat Ui = 0 vs 50, cylinder x periodic,
y open; blocks "try Calabrese-Cardy scaling cut size…" and "GS one-body correlations along x per
row" at the end of `exact-diag/daily-things.jl`).** Aim: Laughlin = gapped bulk + gapless edges,
CDW = gapped everywhere (c → 0). Established: the CDW edge is gapped under strong all-to-all Ui.
- Expected c: S = ((c_L + c_R)/6) log chord, i.e. 1/6 per chiral edge crossed (a chiral boson is
  half a non-chiral one: Dirac = R ⊕ L Weyl, S_R = S_L). The vertical cut crosses both edges
  (R + L) → c_fit = 1, not 2; a cut through one edge only → c_fit = 1/2. Gapped bulk adds only
  area law and −γ. Refs: Holzhey–Larsen–Wilczek NPB 424 (1994); Calabrese–Cardy JSTAT P06002
  (2004); Castro et al. JHEP 07, 114 (2014) (c_L ≠ c_R); Wen PRB 41, 12838 (1990); synthetic
  ladders: Calvanese Strinati et al. PRX 7, 021033 (2017).
- Diagnostics (GS only; excited states too expensive in TTN): S(ℓ) for vertical cuts (x ≤ ℓ) vs
  the chord variable log[(Lx/π) sin(πℓ/Lx)], slope c/3: linear = gapless edge, saturating =
  gapped edge (hybridized or CDW). Row correlations |⟨b†_{r,y} b_{0,y}⟩| vs chord distance: edge
  rows power law vs bulk rows exponential. Bulk gap from torus data, not the cylinder.
- Fitted c (ℓ ≥ 1 / ℓ ≥ 2; Lx = 6, 8, 10; 2–5 points):

    | c fit | Ly = 4, Ui = 0 | Ly = 4, Ui = 50 | Ly = 6, Ui = 0 | Ly = 6, Ui = 50 |
    |---|---|---|---|---|
    | ℓ ≥ 1 | 3.2, 2.7, 2.9 | 1.46, 1.23, 1.06 | 2.76, 2.63 | 1.76, 1.59 |
    | ℓ ≥ 2 | 3.9, 3.0, 3.1 | 1.37, 0.95, 0.74 | 3.06, 2.75 | 1.59, 1.36 |

- Spectra E − E1: Ly = 6 Ui = 0: E2 0.181/0.129, E3 0.199/0.143 (Lx = 6/8), gap·Lx ≈ const (v ≈ 0.17).
  Ly = 4 Ui = 0: E2 0.105/0.0158/0.0155, E3 0.112/0.059/0.056 (Lx = 6/8/10), flattening.
  Ui = 50: TT doublet (Ly = 4: 0.023/0.011/0.0054; Ly = 6: 0.080/0.056) under a gap ≈ 0.2 / 0.1.
- Row correlations, r = 1–4 at Lx = 8: Ui = 0 edge row stays large (8x6: 0.91, 0.72, 0.54, 0.46;
  8x4: 0.69, 0.47, 0.51, 0.49), middle rows fall (8x6 y=3: 0.71, 0.31, 0.12, 0.03). Ui = 50 all
  rows fall alike (≈ 0.6, 0.15, 0.03–0.06, 0.01–0.03).
- Bond currents (2026-10-16, blocks "currents on the cylinder…" and "currents and correlations on
  OBC…" before the CC blocks; earlier versions lacked the Peierls phase). Equilibrium
  (magnetization) currents exist in both phases (max 0.146 vs 0.095 at 8x4), so magnitude isn't
  diagnostic; the shape is.
  - Cylinder 4x8 N=2, J_x(y), y = 1–4 (mirror rows opposite): Ui = 0: −0.015, −0.094, +0.057,
    +0.049 (edge flow + bulk backflow, net ≈ 0 per half); Ui = 50: −0.006, −0.025, −0.026, −0.014
    (whole half flows one way, net −0.071, no separate edge channel). Inner rows reverse at
    Ui ≈ 5, near flat Γ_c (2–4 for Ly = 4–7): possible order-parameter-free marker (unchecked).
    Outer rows nearly empty (n = 0.007).
  - OBC 8x4 N=4: the takeaway is the correlations. Laughlin edge rows coherent over the whole
    length (0.5–0.8 to r = 7), inner rows decay to 0.24; CDW every row exponential, ξ ≈ 1 site
    (×10 per 2 sites), edge as fast as bulk → no edge channel. Currents: Laughlin lumpy density
    with local loops (likely N = 4 OBC finite size); CDW uniform density (GS = superposition of
    the TT patterns), one loop around the perimeter, bulk J_y ≈ 0 (same whole-width flow as the
    cylinder). Ly = 4, two bulk rows: confirm with the saved 8x6 N=4 cylinder densmats.
- **Reading:** CDW gapped everywhere (gap, S saturates, c falls, all rows decay fast). Laughlin at
  Ly = 6 qualitatively right (gap·Lx const, coherent edge vs short-ranged bulk) but not
  quantitative: c ≈ 3 not 1, edge decay much slower than the Laughlin r^-2. Ly = 4 ambiguous
  (gaps flatten, but S(ℓ) at 10x4 still rises and the edge doesn't decay). Too small: ℓ ≤ 4 holds
  0.5–2 particles (spacing 1/ρ1D = 2), so fits are UV-dominated; need Lx ≳ 16 (TTN).
- Wider Ly helps little: at fixed ρ1D, N_φ = Lx for every Ly (no extra orbitals); edge separation
  (Ly−1)/l_B = √(2π(Ly−1)) grows only 4.3 → 5.6 for Ly 4 → 6. Family commensurate for every even
  Lx (α = 1/(Ly−1), phase on x bonds depends on y only, N_φ = 2N; bosons, no N-parity effect).

**Bottom line:** OBC ground-state data for p = 0 (flat), 0.5, …, 3 at 6×Ly N=3 and 8×Ly N=4
(Ly = 4–7), plus 10x4 N=5; p = 4, 5, 6 at 8×Ly only. All p go from the same Laughlin state to the
same large-Γ state; below Γ_c, p ≳ 1.5 takes a detour (precursor), no new phase seen.

**Key result (2026-10-20, user: "this is big"): the larger the exponent, the stronger Γ_c diverges
with system size.** At 8×Ly, Γ_c ∝ Ly^b over Ly = 5–7 with b = 0.84, 1.49, 2.47, 3.66, 5.07, 6.96,
8.72 for p = 0, 1, 2, 3, 4, 5, 6 (S(π); Binder +0.2). Γ_c grows with Ly for every p, including
p = 0 (partly geometry). What's new is that b(p) gets steeper with p instead of saturating:
- b − p ≈ 0.4–0.8 for p ≤ 3. A plausible reading (unchecked): the coupling at range ~Ly,
  Γ_c Ly^−p, grows only slowly.
- b − p jumps to 1.1, 2.0, 2.7 for p = 4, 5, 6.
- For p ≥ 4, Γ_c(Ly) bends **up** (local slope rises with Ly), so the divergence may be faster than
  a power law (3 points can't tell a power law from an exponential).

So for short-range interactions the CDW gets ever harder to reach as the system grows. Caveats:
- Only Ly = 5–7, and Lx = 8 only for p ≥ 4.
- Fixed N: the GS CDW signal collapses with Ly (S(π)/N rise ≈ 0.4, 0.15, 0.06, 0.035 for
  Ly = 4–7, at both Lx).
- Hill vs midpoint Γ_c disagree at Ly ≥ 6.
- Ly = 4 is off the trend and is left out of b.

At Ly = 4, Γ_c does not grow with Lx (6, 8, 10). User (2026-10-14): these are good limited-size
studies, and larger sizes are needed to settle the p dependence (see Next).

- **8×Ly N=4 (cluster sweep, 2026-10-14):** block "cdw transition for 8xi N=4 scaling…" in
  `exact-diag/daily-things.jl` (copy of the 6×Ly block, reads the `all_*` metadata, GS only).
  GS Hill Γ_c from S(π), 6×Ly in brackets. b ± its standard error from the residuals of the
  log-log line, over Ly = 4–7 (S(π), 6×Ly in brackets) and Ly = 5–7 (S(π) / Binder; the
  block's current fit, `fit_ly_min = 5`, Ly = 4 still plotted):

    | Γ_c GS, Hill | 8x4 | 8x5 | 8x6 | 8x7 | b, Ly 4–7 | b, Ly 5–7 |
    |---|---|---|---|---|---|---|
    | p = 0 | 1.95 (1.94) | 3.33 (2.78) | 3.83 (3.41) | 4.42 (3.98) | 1.43 ± 0.30 (1.28) | 0.84 ± 0.05 / 1.08 |
    | p = 0.5 | 2.63 (2.75) | 5.04 (4.23) | 6.13 (5.35) | 7.36 (6.43) | 1.80 ± 0.33 (1.51) | 1.12 ± 0.03 / 1.36 |
    | p = 1 | 3.62 (4.0) | 7.63 (6.31) | 9.91 (8.4) | 12.6 (10.7) | 2.19 ± 0.34 (1.74) | 1.49 ± 0.04 / 1.73 |
    | p = 1.5 | 5.37 (5.86) | 11.5 (9.29) | 16.0 (13.4) | 22.1 (19.0) | 2.48 ± 0.27 (2.09) | 1.94 ± 0.07 / 2.17 |
    | p = 2 | 8.01 (8.6) | 16.9 (13.8) | 26.0 (22.0) | 39.0 (32.6) | 2.80 ± 0.16 (2.39) | 2.47 ± 0.08 / 2.70 |
    | p = 2.5 | 12.0 (12.8) | 25.1 (20.4) | 42.7 (35.2) | 70.2 (56.9) | 3.14 ± 0.06 (2.68) | 3.05 ± 0.09 / 3.27 |
    | p = 3 | 17.2 (18.9) | 37.0 (29.7) | 70.5 (56.3) | 127 (103) | 3.57 ± 0.06 (3.05) | 3.66 ± 0.08 / 3.88 |
    | p = 4 | 40.1 | 82.1 | 203 | 453 | 4.35 ± 0.35 | 5.07 ± 0.07 / 5.26 |
    | p = 5 | 86.7 | 179 | 596 | 1860 | 5.51 ± 0.72 | 6.96 ± 0.22 / 7.13 |
    | p = 6 | 186 | 385 | 1790 | 7270 | 6.61 ± 1.04 | 8.72 ± 0.19 / 8.89 |

  - p = 4, 5, 6 (2026-10-20): 8×Ly only, folders `powerlaw-{4,5,6}p0`. The Γ sweep reaches
    8–11 Γ_c at Ly = 6, 7. The Hill fit has `lower=[1e-3,0.1]` (p = 6, Ly = 4 otherwise walks x50
    negative; no change for p ≤ 3).
  - The Ly = 5–7 errors have 1 dof. They match the ≈ 0.1 that the Hill Γ_c errors alone give.

  - Ly ≥ 5: 8×Ly Γ_c is 10–25% above 6×Ly; still grows with Ly for every p, so not an N = 3
    artifact (Lx = 8 also changes α and Nφ).
  - Ly = 4, Lx = 6, 8, 10, Γ_c S(π) (Binder): p = 0: 1.94, 1.95, 1.73 (1.55, 1.70, 1.56);
    p = 1: 4.0, 3.62, 3.34; p = 2: 8.6, 8.01, 8.14; p = 3: 18.9, 17.2, 17.7 (17.4, 16.6, 16.3).
    Flat or slightly falling with Lx.
  - Ly = 4 transitions are steeper than the Hill shape (the fit misses the upper tail), as at 6x4.
- **Neighbour-difference order parameter (2026-10-20):** O = Σ_{x,y} |⟨n_{x,y}⟩ − ⟨n_{x+1,y}⟩| /
  ((Lx−1) Ly) from the GS occupancies `all_occs1`, block "cdw transition for 8xi N=4 scaling with
  the neighbour-difference order parameter…". The GS ⟨n⟩ is x-reflection symmetric (Lx even maps
  even ↔ odd columns), so the global staggered sum Σ_x (−1)^x ⟨N_x⟩ is exactly 0, but local
  alternation pinned at each open edge is not. O is that one-point, edge-pinned part; S(π) is the
  global staggered correlation (two-point, ⟨N⟩ part 0). Σ_x ⟨(N_x − N_{x+1})²⟩ stays ≈ 6–7 through
  the transition; its disconnected part Σ (⟨N_x⟩ − ⟨N_{x+1}⟩)² follows O. Column profiles ⟨N_x⟩,
  8x4: Laughlin 0.15 0.74 0.52 0.59 | sym. (edge depletion, overshoot, Friedel-like alternation
  from both edges, 50–75% of O on the two edge bonds); CDW 0.50 flat (Ly = 4: symmetric cat of
  the two TT patterns; Ly = 6: residual 0.48 0.52 0.49 0.51). **Precursor (p = 3, Ly = 4, Γ = 3–10):**
  0.14 0.82 0.35 0.70 | sym., i.e. each edge pins an opposite TT pattern (even Lx) with a domain
  wall at the centre: O peaks and S(π) dips because the wall cuts the staggered correlation. At
  Γ_c the GS switches to the flat cat: O drops, S(π) jumps. Edge-dominated, so O ~ 1/(Lx−1) at
  large Lx (boundary observable, unchecked). It still tracks the transition: O falls where S(π)
  rises and 1 − O (normalized) nearly collapses onto S(π) for Ly ≥ 5; at Ly = 4 O is a sharp
  step (GS level crossing) while S(π) keeps growing up to Γ ~ 10³. Falling Hill Γ_c(O) vs Γ_c(S(π)): same within 2–7% at Ly = 4 (all p), 5–20% high at
  Ly = 5, 6, matches again at Ly = 7 for p ≥ 1.5, but for p ≤ 1 at Ly = 7 it is non-monotonic in Ly
  (p = 0: 1.51, 4.04, 4.91, 3.49; no O peak, slow fall from Γ → 0). Ly = 4 large Γ (≳ 10³):
  O jumps back up (near-degenerate TT doublet mixed arbitrarily, breaks reflection), fit stops at
  O's minimum. O ≡ 0 on the cylinder (⟨n⟩ uniform in x).
- **Single power law Γ_c ∝ Ly^b over Ly = 4–7 is a poor fit (local exponents, 2026-10-14/20):**
  at 8×Ly the 4→5 step has slope 2.4–3.4 for every p (≈ 3.2 for p = 4–6), so Ly = 4 is off the
  trend for both small and large p. After that step the slope drops for p ≤ 2, stays flat for
  p = 2.5–3, and rises for p ≥ 4 (5→6: 5.0, 6.6, 8.4; 6→7: 5.2, 7.4, 9.1 for p = 4, 5, 6). The
  Ly = 4–7 b error bars show it: ±0.3 for p ≤ 1.5, ±0.06 for p = 2.5–3, ±0.35–1.0 for p ≥ 4.
  Hence b over Ly = 5–7 (table; b < 1 for small p). Midpoint Γ_c instead of Hill gives b(5–7) =
  1.26, 1.69, 2.09, 2.52, 3.00, 3.55 for p = 0–2.5 (never < 1). A local slope carries ±0.2–0.3, so "b crosses 1" isn't resolved. Both
  give b(p) − b(0) still growing past p = 1 (midpoint ≈ 0.85 p): no saturation or kink at p = 1,
  where Σ r^−p starts converging in 1D. The naive picture (pair spread over Ly feels
  Γ Ly^−min(p,1)) predicts saturation; the data don't show it, which fits the pairs separating in
  y (precursor). Hill vs midpoint differ mainly in the baseline (dip minimum vs Γ → 0 value).
  With p = 4–6 (2026-10-20), b(5–7) steepens from ≈ 0.6 to ≈ 1.8 per unit p (see Bottom line).
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

**Edge thread → TTN (2026-10-16).** Ly = 6, Lx = 12, 16 (+ Ly = 8 width check), Ui = 0 and 50,
GS only. First check which vertical cuts ℓ are tree bonds (compute all S(ℓ) on the same footing).
Then per size: S(ℓ) vs chord, c in a fixed ℓ/Lx window (e.g. 0.1–0.5) with fit rms, S_2 from the
Schmidt values (slope c/4, must agree with S_1), row correlations (edge vs bulk). Laughlin result
= c → 1 with power-law edge / exponential bulk; c → 0 means hybridized edges (go wider).
Single-edge cut A = {x ≤ ℓ, y ≤ Ly/2}: S = (1/6) log chord + a ℓ + const (fit all 0 < ℓ < Lx,
no ℓ ↔ Lx − ℓ symmetry); independent chiral edges ⇒ half the full-cut slope (c_fit 1/2 vs 1),
CDW both → 0. Not in ED (bulk 2–3 rows, 3 parameters on 3–5 points).
If gapped, topology is in the saturated S ≈ 2(a Ly − γ) → −2γ vs Ly (γ = ln√2 vs 0), needs
several Ly and a single-sector GS. Dropped: Rényi-2 in ED, Ly = 4 orbital-counting question.

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
  Γ_c); Γ = 0.1–1000 vs t = 1 is stiff. (Since addressed in the script: hardcore `"S=1/2"` sites,
  sampling restricted to the N sector via `max_counts`.) In its favour: diagonal long-range terms are cheap in VMC.
  First step for fPEPS would be an 8x4/8x5 ED benchmark (N, hardcore, D = 2–6, Γ_c).

Smaller options: Γ_c(Lx) figure at Ly = 4 (6x4, 8x4, 10x4); p = 4–6 at 6×Ly N=3 (cheap; second
Lx for the large-p divergence); Ly = 8 at 8×Ly N=4 for p ≥ 3 (power law vs exponential); save energies (`all_nrgs`) in future sweeps to check E2 − E1 at Lx = 8.

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

**p = 4, 5, 6 at 8×Ly N=4 (runner ready, not run, 2026-10-20):** active block "Γ sweeps on 8×Ly N=4
open strips for p = 4, 5, 6" in `exact-diag/execute-ed.jl`, guarded by `abspath(PROGRAM_FILE) == @__FILE__`
(runs only as the script, not when included). `--array=0-11` (sizes outer, 0–2 = 8x4), `julia
--threads=$SLURM_CPUS_PER_TASK execute-ed.jl "which_one" $SLURM_ARRAY_TASK_ID`, 14 cpus for 28 Γ points.
Γ points are threaded: getHamiltonian locked (JLD2 reads + the global `dressed_hopping_cache` aren't
thread-safe), eigsolves parallel, BLAS 1 thread. Writes the returned-sweep format (one file per (size, p),
named after the largest Γ, H/basis `nothing`) once at the end into `obc/powerlaw-{4p0,5p0,6p0}`
(must exist on the cluster), plus `all_nrgs`, `all_converged`. Grid: 16 coarse over 0.1–1e4 + 12 in
[Γ_g/10, 10 Γ_g]; Γ_g from a no-save 6×Ly N=3 31-point scan: p = 4: 40, 73, 152, 1400; p = 5: 86, 131,
421, 2020; p = 6: 182, 275, 1440, 7730 (Ly = 4–7; the scan's p = 3 gives 19, 30, 208, 392 vs fine-grid
18.9, 29.7, 56, 103, so Ly ≥ 6 is rough). Tested on 6x4 N=3 into the scratchpad (buildHam in place of
the reads): file read by the 8×Ly block, GS S(π) matches run_normal_ed to 1e-13. Add 4.0, 5.0, 6.0 to `ps`
when back.

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
- CDW order parameters, shared by ED and fPEPS in `other-funcs/basic-2d-observables.jl`:
  `column_cdw_statistics(ncols, probs)` (⟨n_x n_x'⟩, P(k on even columns) from column occupations
  `ncols[x,i]` of configurations with probabilities) and `cdw_order_parameters` (S(π), correlation
  ratio, Binder). Each method only supplies `column_occupations` + probabilities. Saving:
  `save_column_cdw_statistics(filepath, corr, dist, k)` / `read_column_cdw_statistics(metadata, k; if_sweep)`
  with the existing metadata keys (`column_corr<k>`, `evencolumn_dist<k>`; sweeps `all_columncorrs<k>`,
  `all_evencolumndist<k>`); old ED blocks still write their dicts inline.
- `exact-diag/observables.jl`: `get_manifold_occupancy` (→ `occs_eig1/2`); `column_occupations(full_basis, Lx)`
  and the ED wrapper `column_cdw_statistics(wavefunc, lattice_params)` (same call as before).
- `synth-dims/ed-plus-ttns.jl`: ED+TTN torus analysis; scaling blocks at the end (s_c/l_B plateau,
  Γ with V = 300, 50, boundary-interaction check), all `#= =#`.
- `synth-dims/dd-paper-plots.jl`: paper figures (so far: shifted Zeeman minima + 1/(s r)³ profile).
- `fpeps/execute-fpeps.jl`: fPEPS (QuantumNaturalfPEPS + QuantumNaturalGradient, env
  `--project=fpeps`), Hamiltonian from `fpeps/build-ham.jl` (`long_range_HH_ham`, hardcore bosons
  on S=1/2, takes ED key names `if_periodic_x/y`). Currently a 4x4 N=2 flat Ui = 50 test, D = 6,
  20 iterations; goal is larger open systems. Saves like the TTNs (2026-10-20, see Gotchas); saving,
  the per-iteration `saving_fpeps_callback` and `continue_fpeps(filepath, extra_iters)` (resumes a
  saved run in the same files) are in `fpeps/savings.jl` (2026-10-29). Loading: `fpeps/load-fpeps.jl`
  is the hub (packages + other-funcs + build-ham/savings/observables via `include_other_files`, so
  pwd inside fzj/ or local/<user>/); library files include nothing, drivers load
  `include-other-files.jl` then the hub.
  With open boundaries the 4x4 N=2 test trips `check_fluxes` (α = 2/(0.5·3·3) = 0.44, "Alpha is
  too large"); add `"if_check_fluxes" => false` to the params dict to skip the check.
- `fpeps/observables.jl` (no ED includes): `sample_fpeps_occupations` (`get_sample` + importance
  weights |ψ|²/p_c → `occs[x,y,i]`), `column_occupations(occs)`, `column_cdw_statistics(occs, weights)`;
  `fpeps_cdw_order_parameters(occs, weights)` adds jackknife errors (one pass over bins),
  `(filepath, nsamples; if_save_data)` samples a saved run and optionally saves into its fpeps-*.jld2
  metadata, `(metadata)` redoes the jackknife from saved bins. `output_level` 1/2 prints run info / sampling
  progress every 10% with time left (default 0, silent). `exact_fpeps_occupations(peps, N)` enumerates the sector
  exactly (small systems only).
- `fpeps/daily-things.jl`: fPEPS scratch blocks; first block tests the above against exact enumeration.
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
  - ED `cylinder/flat-ints/` (2026-10-15/16): flat Ui = 0 and 50, α = 1/3: 6x4 N=3, 8x4 N=4,
    10x4 N=5; α = 1/5: 6x6 N=3, 8x6 N=4 (nev = 10). States in data; `ee_cuts` (all) and GS
    `densmat` (not 10x4) in metadata; `full_basis` and `H` are `nothing`.
  - Lx sweep (2026-10-14), all 7 p, in the same `obc/<folder>`: 8×{4,5,6,7} N=4 and 10x4 N=5, one
    file per (size, p), named after the smallest s (flat: V = 0.1), whose metadata holds the sweep:
    `all_gammas`, `all_columncorrs1/2`, `all_evencolumndist1/2`, `all_occs1/2` (25 points; no
    energies). 8×Ly 12–120 MB, 10x4 210 MB. Read the metadata only (`jldopen`, as the 8×Ly block).
  - fPEPS `cluster-data/synth-dims/fpeps/` (default `dataloc` of `execute-fpeps.jl`): empty so far.

## Open questions

- [ ] Does the exponent p change the physics? Strongest hint so far (2026-10-20): b(5–7) − p
  ≈ 0.4–0.8 for p ≤ 3, then 1.1, 2.0, 2.7 for p = 4–6. Γ_c(Ly) curves down for small p and up for
  p ≥ 4. Is the large-p divergence faster than a power law (exponential in Ly)? Does it survive
  at fixed density (Lx growing with N)? Needs more sizes → which family and method (Next).
- [ ] Why does the GS CDW signal at Γ → ∞ collapse with Ly (S(π)/N rise 0.40 → 0.035 for
  Ly = 4–7)? Same at N = 4 / Lx = 8, so it's Ly, not N. Is Γ_c(Ly) at fixed N meaningful there?
- [ ] GS-only: Hill and midpoint Γ_c disagree at Ly ≥ 6 (small signal). Which to use / show both?
- [ ] Γ_c converged in Lx at fixed Ly? Ly = 4: yes (Lx = 6–10). Ly ≥ 5: 8×Ly 10–25% above 6×Ly
  (only two Lx). Binder crossing in Lx with MPS/TTN?
- [ ] What is the pre-transition dip in S(π)/Binder? Known: it deepens with p and Ly, is absent
  for p ≲ 1, and comes with same-column pairs separating in y; in 6x4 there are also true GS
  level crossings near Γ → 0 and Γ_c. Candidate (2026-10-20, 8x4 p = 3 column profile): opposite
  TT patterns pinned at the two open edges with a central domain wall. Check at other Ly/p and
  at odd Lx (no frustration, so the dip should vanish), and whether a dip exists on the cylinder.
- [ ] Is state 2 also excited (not a CDW partner) at 8×Ly, Ly ≥ 5? No energies in the sweep files.
- [ ] V = 50 TTN data at 12x6+ for the torus?
- [ ] Does fPEPS reach useful open sizes / accuracy (N conservation, D convergence of Γ_c)?
- [ ] C⁴ / PES vs s at TTN sizes (the paper's more robust diagnostics).
- [ ] Cylinder Laughlin: does c → 1 (and 1/2 for the single-edge cut) and the edge correlation
  become a power law at TTN Lx (Ly = 6)? ED (Lx ≤ 10) is qualitatively right but UV-dominated (c ≈ 3).
- [ ] Does the J_x sign reversal of the inner rows track Γ_c (cylinder, larger N, dd vs flat)?

## Gotchas

- `cylinder/flat-ints`: `full_basis` is `nothing` in the metadata, rebuild with
  `get_lattice_params_from_metadata(m)` (same basis order, checked against `ee_half`). 8x4
  `ee_cuts` has zeros for ℓ > Lx/2 (saved before mirroring; 6x4 is mirrored): mirror on reading.
  GS one-body density matrix is `m["densmat"]` (metadata); `d["densmat"]` is an empty vector.
- `physical_current` / `synthetic_current` (`exact-diag/observables.jl`) now take
  `(densmat, lattice_params, hamilt_params)` and return 2t Im(e^{iφ} ρ[s, s+e]) with the flux and
  twist phase (via `bond_current`); the old versions omitted the phase and divided by n_i + n_j.
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
- `find_eigenstates`' default Lanczos (krylovdim 30, maxiter 100) leaves the upper of 20 levels
  unconverged at large Γ: 6x4 N=3 p = 4 vs dense, Γ = 1e3: E15+ off by 1e-3–0.1; Γ = 1e4: E1–3 exact,
  E4–7 ~1e-6, E8+ 0.03–27. GS / state-2 observables unaffected. `Lanczos(krylovdim=60, maxiter=1000)`
  converges all 20 (used in the p = 4–6 runner). Saved `nrg` beyond the lowest few at large Γ in older
  files is suspect.
- `write_data_jld2` never overwrites: an existing name gets `-mk-N` appended (`prep_file`), and a
  missing target folder silently writes into `pwd()`. Write to a fresh name and `mv`; check the
  folder exists (don't create folders: the user manages them).
- `flat-ints` files store the defaults `magnetic_spacing = 1.0`, `exponent = 3.0` in metadata:
  for p = 0 use Γ = `interaction_strength`, and don't filter `find_data_file` on V.
- `"powerlaw"` filenames have a different key order than the dd ones (`ed-hopping_anisotropy-…`);
  `find_data_file` matches them fine.
- Hand-excluded TTN points: 14x7 s = 20.0, 12x6 s = 1.53.
- fPEPS files are JLD2 (the PEPS has no HDF5 serializer), split like the TTN ones:
  `wavefuncfpeps-<params>.jld2` holds `peps` (double-layer envs stripped, they're a rebuildable
  cache), `fpeps-<params>.jld2` holds `history` (DataFrame), `energy`, `rng` + metadata (params,
  `bond_dim`, `sample_nr`, `maxiter`, `learning_rate`, `niter`, `save_every`). Both are overwritten
  by the `evolve` callback every `save_every` iterations (params_dict key, default 10) and at maxiter, so
  a killed job loses < save_every iterations; that's what `evolve(...; misc_restart)` needs to resume.
  Find with `find_data_file(pdict, "fpeps", dataloc; file_type="jld2")` (default is `"h5"`).
  CDW statistics in fpeps metadata (GS, k = 1): `column_corr1`, `evencolumn_dist1` (all samples, ED
  names) plus `column_corr_bins1`, `evencolumn_dist_bins1`, `bin_weights1` (jackknife), `nsamples1` and
  `niter_sampled1`: if it differs from `niter`, the statistics predate a `continue_fpeps` and are stale.
  Learning rate is `learning_rate` in metadata: `"lr"` is the long-range setting.
  `evolve`'s `maxiter` counts restored iterations too (stops at `niter > maxiter`), so
  `continue_fpeps` passes `maxiter = niter + extra_iters`; it rebuilds the Hamiltonian and an Euler
  integrator from the metadata (only the integrator/solver type names are saved).
  Since 2026-10-29 metadata = `get_normal_model_paras` output (inputs `lr`, `ts`, `filling`,
  `if_check_fluxes`, scaling inputs + derived `t_x`, `long_range_strengths`, ...), so re-running it on
  the metadata gives the same model, for every scaling. Filenames now carry `lr` (+ `magnetic_spacing`
  for dd/powerlaw, `ts` if ≠ 1, `which_dir` if ≠ y). Earlier runs in `fpeps/flat-ints`: the lr-all
  46.4159 one has raw params only (`continue_fpeps` re-derives); the 18.90527 and 4x4 test files
  have `"N"` instead of `"particles"` and won't load (user: leave them).
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
- 2026-10-15 — Edge / central-charge side thread started (`exact-diag/daily-things.jl`, end):
  cylinder vertical-cut Calabrese–Cardy fits; `compute_reduced_density_matrix` rewritten as
  ρ_A = Ψ Ψ† with sparse Ψ (same ρ_A to 1e-16).
- 2026-10-16 — CC fits in one grid; Ly = 6 ED (6x6, 8x6); spectra and row correlations: CDW
  gapped everywhere, Laughlin qualitatively gapped bulk + gapless edge at Ly = 6 but UV-dominated
  (c ≈ 3); Ly = 4 ambiguous. Plan: TTN Ly = 6, Lx = 12, 16, GS-only diagnostics (Next).
- 2026-10-16 — Current functions fixed (Peierls phase); current blocks cleaned. CDW: one
  whole-width circulation, no edge/bulk counterflow, edge coherence as short as bulk (ξ ≈ 1).
  c_fit = 1 for Laughlin explained (two chiral edges, 1/6 each); single-edge cut test planned.
- 2026-10-20 — p = 4, 5, 6 (8×Ly N=4) added to the 8×Ly block (bounded Hill fit, skips missing
  sizes). Ly = 4 left out of the b fit (still plotted). **Key result:** Γ_c diverges with Ly ever
  more strongly as p grows: b(5–7) = 5.07, 6.96, 8.72 for p = 4, 5, 6 (b − p: 1.1, 2.0, 2.7 vs
  0.4–0.8 for p ≤ 3), and Γ_c(Ly) bends up for p ≥ 4.
- 2026-10-20 — fPEPS saving added (TTN-style wavefunc/metadata split, JLD2, every iteration,
  `cluster-data/synth-dims/fpeps/`); `fpeps` prefixes in `data-storage-funcs.jl`, JLD2 in the
  fpeps env. Tested end to end (2 iterations into the scratchpad). fPEPS dicts switched to
  `if_periodic_x/y` (build-ham ignored `_phys/_synth`, so runs were always open). `continue_fpeps`
  added, tested 5x5 N=2 D=2 (2 + 2 iterations into the scratchpad, history/niter carry over).
- 2026-10-20 — Neighbour-difference order parameter O from GS occupancies (8×Ly N=4, all p): tracks
  the transition via the open-x density profile (staggered ⟨n⟩ = 0 by reflection), Γ_c close to
  S(π) at Ly = 4 and for large p, unreliable for small p at Ly = 7.
- 2026-10-20 — Threaded cluster runner for p = 4, 5, 6 at 8×Ly N=4 in `execute-ed.jl` (not run yet);
  Γ_c guesses from a 6×Ly coarse scan; found the default Lanczos leaves upper levels unconverged at large Γ.
- 2026-10-27 — fPEPS CDW order parameters (`fpeps/observables.jl`, test in `fpeps/daily-things.jl`):
  random 4x4 N=2 D=3 PEPS, 16000 samples match exact enumeration within 1σ (S(π) 1.081 ± 0.014 vs
  1.094, Binder 0.384 ± 0.008 vs 0.391). The combinatorial sector proposal gives n_eff ≈ 0.4 × samples.
- 2026-10-27 — CDW statistics/order parameters moved to `other-funcs/basic-2d-observables.jl` (ED and
  fPEPS only supply column occupations). ED: identical to the old code and saved metadata to 1e-14
  (29 OBC files, Lx sweep), 60× faster (10x4 N=5: 0.07 s vs 4.4 s). fPEPS jackknife identical, 75× faster.

- 2026-10-27 — Common save/read for the column statistics (`basic-2d-observables.jl`, existing ED key
  names) and fPEPS saving with per-bin statistics; tested on scratchpad copies of ED point/sweep files
  and a scratchpad fPEPS file (round trips exact, other keys untouched, no -mk files).
