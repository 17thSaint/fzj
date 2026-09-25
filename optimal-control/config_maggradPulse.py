"""QuOCS dCRAB optimization of the time-dependent magnetic gradient pulse.

The gradient B(t) displaces the synthetic states in a trap of frequency w, so the
spacing entering the dipolar tail follows the driven response

    a(t) = a0 + (1/w) * int_0^t B(t') sin(w(t - t')) dt'

(exact-diag/time-evolution.jl: get_magnetic_gradient_integral, applied in the
"magnetic_gradient" branch of long_range_scaling). The figure of merit asks for the spacing
averaged over the last `window` of the run to sit at `a_target`, with a penalty on how much
it is still ringing there:

    cost = |mean(a over window) - a_target| + std(a over window)
    FoM  = -cost                                   (dcrab_algorithm_settings maximizes)

B is free over [0, drive_time]; anything between drive_time and tmax is a hold with the
gradient switched off. drive_time defaults to tmax, i.e. no hold -- the gradient stays
available right through the averaging window. Setting drive_time < tmax reinstates one.

Unlike the other configs here the FoM is classical -- no ED, no time evolution -- so it is
pure numpy and runs in ~1 ms. Julia is called once at setup, only to get the RK4 time step
(maggrad-pulse-control-functions.jl: setup_maggrad_pulse), because the pulse must be
sampled on the same half-step grid that a later run_timeevo will index it on.

Two things about this cost are worth knowing before reading a result (both are measured
and printed by main()):

* It is a *linear* functional of B, so the target-distance term alone is massively
  degenerate -- infinitely many pulses hit it exactly. The std term is what picks a
  particular one out, and the two enter unweighted, both being spacings.

* There are two separate caps, and which one bites decides whether the spacing SETTLES or
  merely averages right. Pumping the mode resonantly caps the reachable *mean* at roughly
  b_max * drive_time / (2 w) times a window-averaging factor. But holding the spacing
  steady, with no ringing left, needs a standing gradient B = w^2 * (a - a0), so the
  largest spacing that can be held at all is a0 + b_max/w^2 -- and that one is usually far
  more restrictive. At b_max = 30, w = 10 it is only 0.400.

Measured at a_target = 2.0, w = 10, over choices of b_max and of the cost itself:

    b_max   cost                          mean      std
       30   distance + std              2.0000   0.4391    ringing +-22%, cap is 0.40
       30   distance + 5*std            0.4000   0.0000    gives up on the target entirely
      200   distance only               2.0000   6.9648    348% -- the std term is load-bearing
      200   distance + std              2.0000   0.0198    <- defaults here, 20000 evals
      200   distance + std, 5x budget   2.0000   0.0021

At b_max = 30 the two terms are irreconcilable and no weighting helps: the solution jumps
discontinuously between "mean on target, ringing hard" and "steady at the 0.400 cap", with
nothing in between. Once b_max clears w^2 * (a_target - a0) = 190 they stop competing, an
unweighted sum settles the spacing to within ~0.1% and the remaining error is just
optimizer budget. The pulse that does it rings B hard to null the mode's velocity, then
parks near +190 for the last stretch, which is what flattens a(t) onto the target.

Zero crossings of a(t) are the other thing to watch, and they are controlled by the initial
guess far more than by the cost. Because U ~ 1/a^3, a trajectory that crosses zero makes the
interaction diverge and the eventual time evolution impossible; the window terms say nothing
about the route taken to get there, so the guess decides it. Measured at a_target = 2.0,
b_max = 200 (std / min(a) / number of zero crossings):

    drive T  floor  guess                 std      min(a)  crossings
        2.5    off  resonant pump      0.0120     -3.9886          6
        2.5    off  ramp to hold       0.0001      0.1000          0
        2.5     on  resonant pump      0.0293      0.1000          0
        2.5     on  ramp to hold       0.0005      0.1000          0
        1.0    off  resonant pump      0.0051     -2.5822          2
        0.5    off  resonant pump      0.0396      0.0844          0
        0.5    off  ramp to hold       0.0024      0.1000          0

Shortening the run does help -- the room for pumping excursions scales as b_max*T/(2w), so
by T = 0.5 even the resonant guess stops crossing zero -- but the guess is the stronger
lever and fixes it at every duration. min_spacing_penalty also removes the crossings, but as
a blunter instrument: it fights a bad guess rather than avoiding one, and costs an order of
magnitude in std. The defaults here use the ramp guess with the floor on as insurance, and
reach mean 2.0000 / std 0.0000 / min(a) = 0.100000 with a(t) never dipping below a0 -- so
max(U) over the whole run equals U(a0) and the t = 0 RK4 step stays valid throughout, which
is exactly what tevo-daily-things.jl's dt choice assumes.

For reference the analytic optimum is a two-step pulse: B = w^2*(a_target-a0)/2 for half a
trap period (pi/w), which brings the mode to rest exactly at the target, then B =
w^2*(a_target-a0) to hold it. That settles by t = 0.314 with std ~1e-5 and no overshoot at
all, so a run much longer than a trap period is not needed for the spacing's sake.

Note for reinstating a hold (drive_time < tmax): with B pinned to zero through the window
the mean and variance terms become strictly incompatible, since a free oscillator only
stops ringing at x = 0. A floor of a >= a0 then collapses the optimizer to B = 0 outright.
"""

import os

import numpy as np

from quocs_common import (FIGURES_DIR, JuliaFoM, best_pulse, dcrab_algorithm_settings,
                          fourier_pulse, halfstep_bins, include_julia, jl,
                          load_best_controls, run_optimization)

include_julia("maggrad-pulse-control-functions.jl")

PULSE_NAME = "maggradPulse"
TIME_NAME = "time_maggradPulse"


def retarded_response(b: np.ndarray, step: float, trap_frequency: float,
                      cos_wt: np.ndarray, sin_wt: np.ndarray) -> np.ndarray:
    """(1/w) int_0^t B(t') sin(w(t - t')) dt' at every sample of a uniform grid of spacing step.

    Uses the same trapezoid samples and weights as get_magnetic_gradient_integral
    (exact-diag/time-evolution.jl), but splitting the kernel as
    sin(w(t-t')) = sin(wt)cos(wt') - cos(wt)sin(wt') pulls the upper-limit dependence out of the
    integrand. What is left are two ordinary cumulative integrals, so the whole retarded response
    costs O(n) rather than redoing an O(n) sum at each of the n upper limits -- which matters,
    because n is ~25000 here. The numpy twin of magnetic_gradient_response (time-evolution.jl),
    kept in numpy so an evaluation stays ~1 ms; cos_wt / sin_wt are the kernel factors on the grid.
    """
    bc, bs = b * cos_wt, b * sin_wt
    cum_c = np.concatenate(([0.0], np.cumsum(0.5 * step * (bc[:-1] + bc[1:]))))
    cum_s = np.concatenate(([0.0], np.cumsum(0.5 * step * (bs[:-1] + bs[1:]))))
    return (sin_wt * cum_c - cos_wt * cum_s) / trap_frequency


class maggradPulse(JuliaFoM):

    def __init__(self, args_dict: dict = None):
        if args_dict is None:
            args_dict = {}

        # hamiltonian / lattice parameters. These only set the RK4 step, not the cost
        self.Lx = args_dict.setdefault("Lx", 4)
        self.Ly = args_dict.setdefault("Ly", 4)
        self.N = args_dict.setdefault("N", 2)
        self.lr = args_dict.setdefault("lr", "all")
        self.if_periodic_x = args_dict.setdefault("if_periodic_x", True)
        self.if_periodic_y = args_dict.setdefault("if_periodic_y", True)
        self.interaction_strength = args_dict.setdefault("interaction_strength", 10.0)
        self.speccount = args_dict.setdefault("speccount", 2)

        # the driven-response model
        self.a0 = args_dict.setdefault("a0", 0.1)                     # spacing at t = 0
        self.trap_frequency = args_dict.setdefault("trap_frequency", 10.0)

        # timing. The pulse QuOCS optimizes covers [0, drive_time]; the remaining
        # tmax - drive_time is the hold, where B is pinned to zero and only the ringing runs
        self.drive_time = args_dict.setdefault("drive_time", 2.5)
        self.tmax = args_dict.setdefault("tmax", 2.5)
        self.window = args_dict.setdefault("window", 0.25)            # averaged over the end

        # the cost
        self.a_target = args_dict.setdefault("a_target", 2.0)
        # Floor on a(t) over the WHOLE run. This is not cosmetic: the interaction scale goes
        # as U ~ 1/a^3, so a trajectory that dips toward zero spacing makes the eventual time
        # evolution arbitrarily stiff, and one that CROSSES zero makes U diverge outright.
        # The window terms only constrain the end of the run, so without this the optimizer
        # is free to pump through huge excursions on the way there -- which it does
        self.min_spacing = args_dict.setdefault("min_spacing", self.a0)
        self.min_spacing_penalty = args_dict.setdefault("min_spacing_penalty", 0.0)

        # one ED, only for the RK4 stability limit at the t = 0 (tightest, if B >= 0) profile
        self._setup = jl.setup_maggrad_pulse(self.to_julia_dict())
        self.dt = float(self._setup.dt)

        # the RK4 half-step grid (spacing dt/2) that pulse_ramp samples on and timeham
        # indexes directly; the optimized pulse occupies its leading drive_bins entries
        step = self.dt / 2
        self._t = np.arange(int(np.ceil(self.tmax / step)) + 1) * step
        self.drive_bins = halfstep_bins(self.drive_time, self.dt)
        assert self.drive_bins <= self._t.size, "drive window does not fit inside tmax"
        assert np.isclose(self._t[self.drive_bins - 1], self.drive_time), \
            "drive_time is not on the half-step grid, so the pulse cannot be zero-padded"

        # kernel factors and the averaging window, fixed for the whole optimization
        self._cos = np.cos(self.trap_frequency * self._t)
        self._sin = np.sin(self.trap_frequency * self._t)
        self._in_window = self._t >= self.tmax - self.window
        assert self._in_window.any(), "averaging window contains no grid points"

    def spacings(self, drive_pulse: np.ndarray) -> np.ndarray:
        """a(t) on the full tmax half-step grid, for a gradient given over the drive window
        (retarded_response; verified against the naive rule to ~1e-19 in main())."""
        b = np.zeros(self._t.size)
        b[:drive_pulse.size] = drive_pulse            # the gradient is off through the hold
        return self.a0 + retarded_response(b, self.dt / 2, self.trap_frequency, self._cos, self._sin)

    def cost(self, drive_pulse: np.ndarray) -> float:
        a = self.spacings(drive_pulse)
        in_window = a[self._in_window]
        # unweighted sum: the distance and the residual ringing enter on equal terms, both
        # in units of spacing. There is no coefficient to tune -- once b_max clears
        # w^2*(a_target - a0) the two are not in conflict and a weight would buy nothing
        cost = abs(in_window.mean() - self.a_target) + in_window.std()
        if self.min_spacing_penalty:
            # time-averaged violation rather than the single worst point: the depth of the
            # excursion and how long it lasts both matter for the eventual RK4 step, and a
            # mean gives Nelder-Mead a signal that keeps improving as the violating stretch
            # shrinks, where a bare min() only ever reports its deepest sample
            violation = np.maximum(0.0, self.min_spacing - a)
            cost += self.min_spacing_penalty * violation.mean()
        return float(cost)

    def get_FoM(self, pulses: list = [], parameters: list = [], timegrids: list = []) -> dict:
        drive_pulse = np.real(np.asarray(pulses[0])).ravel()
        return {"FoM": -self.cost(drive_pulse)}


# must match the maggradPulse defaults above; dt is not here because it comes from the ED
b_max = 200.0
trap_frequency = 10.0      # must match the maggradPulse defaults
a0 = 0.1                   # spacing at t = 0
tmax = 2.5
drive_time = tmax          # no hold: the gradient is free right through the averaging window
window = 0.25
a_target = 2.0

# The pulse QuOCS builds is (base + fourier_update) * scaling(t) + initial_guess(t), clipped
# to [-b_max, b_max]. So the guess is what the optimizer starts from AND what the frozen
# endpoints are held at, while the scaling picks which endpoints those are.
#
# The standing gradient that holds a_target with no ringing at all. Everything here is
# organised around it: it has to fit inside b_max, and it is where the pulse should end up
hold_gradient = trap_frequency ** 2 * (a_target - a0)
assert hold_gradient <= b_max, (
    f"a_target = {a_target} needs a standing gradient of {hold_gradient}, above b_max = "
    f"{b_max}; the spacing cannot be held steady there (see the module docstring)")

# A smooth ramp from zero up to that holding value. This matters far more than it looks:
# the window terms only constrain the END of the run, so the route taken to get there is
# whatever the guess biases the optimizer toward. A resonant guess (b_max/2 * sin(w*t)) pumps
# the mode, and the optimizer then polishes a trajectory that swings a(t) to -4 and through
# zero six times -- fatal for a later time evolution, where U ~ 1/a^3 diverges at a crossing.
# This guess instead drives the mode quasi-statically, so a(t) rises monotonically from a0
# and never dips below it. Measured at tmax 2.5 with no floor penalty: resonant guess gives
# min(a) = -3.99 with 6 zero crossings and std 0.0120, this one gives min(a) = 0.1000 with 0
# crossings and std 0.0001. It must also not be bang-bang -- a guess sitting on the amplitude
# limits would have every update clipped straight back off by limit_pulse
initial_guess_lambda = f"lambda t: {hold_gradient} * (t / t[-1])"
# Pin only t = 0, leaving B(tmax) free. With no hold the gradient should still be on at the
# end -- that is the only way to hold a displaced spacing steady -- so freezing the final
# value at the guess (which the default parabolic scaling in fourier_pulse would do) throws
# away the low-variance solutions entirely. Measured: with both ends pinned the optimizer
# stalls at a window mean of 2.6-3.4, with the end free it reaches 3.62 against a ceiling of
# 3.63, and which guess is used then barely matters (zero / half-sine / resonant sine all
# land within 0.1% of each other). The endpoint pinning, not the guess, was the constraint
scaling_lambda = "lambda t: t / t[-1]"

fom_params = {
    "a_target": a_target,
    # cheap insurance against the optimizer wandering back into large excursions: costs
    # ~0.0004 of std with the ramp guess, and forbids the zero crossings outright
    "min_spacing": a0,
    "min_spacing_penalty": 50.0,
    "drive_time": drive_time,
    "tmax": tmax,
    "window": window,
}


def build_optimization_dictionary(dt: float, drive_time: float = drive_time) -> dict:
    """The pulse's bins_number follows from the RK4 step, which only the ED knows, so
    unlike the other configs this cannot be a module-level constant.

    drive_time is a parameter rather than read from module scope so that a caller
    overriding it through fom_params gets a pulse of the matching length -- otherwise the
    bins_number and the FoM's grid silently disagree.
    """
    return {
        "optimization_client_name": "maggradPulse_dCRAB",
        # an evaluation is ~1 ms here (no ED, no propagation), so the budget is set by
        # patience rather than by cost -- this is ~100 s and leaves dCRAB room for the
        # super iterations it needs to redraw the Fourier frequencies toward the resonance.
        # Once b_max clears w^2*(a_target - a0) the residual ringing is budget-limited
        # rather than physics-limited: 20000 evals leave std ~0.02, 100000 leave ~0.002
        "algorithm_settings": dcrab_algorithm_settings(super_iteration_number=20,
                                                       max_eval_total=100000),
        "pulses": [
            # The Fourier frequencies are 2*pi*w_i*t/drive_time with w_i ~ U(0.01, 10), so
            # the trap resonance (w_i = trap_frequency*drive_time/(2*pi) ~ 4.0) is inside
            # the basis' reach -- which is what lets dCRAB pump the oscillator at all
            fourier_pulse(pulse_name=PULSE_NAME,
                          time_name=TIME_NAME,
                          ramptime=drive_time,
                          dt=dt,
                          lower_limit=-b_max,
                          upper_limit=b_max,
                          amplitude_variation=0.3 * b_max,
                          initial_guess_lambda=initial_guess_lambda,
                          scaling_lambda=scaling_lambda),
        ],
        "parameters": [],
        "times": [{"time_name": TIME_NAME, "initial_value": drive_time}],
    }


def check_response_implementation(fom: maggradPulse, n_coarse: int = 2000) -> float:
    """Reproduce get_magnetic_gradient_integral's O(n^2) rule on a coarse grid and compare
    against the cumulative form retarded_response uses. Cheap insurance that the two agree."""
    step = fom.dt / 2
    t = np.arange(n_coarse) * step
    b = np.random.default_rng(0).normal(size=t.size)
    naive = np.zeros(t.size)
    for k in range(1, t.size):
        integrand = b[:k + 1] * np.sin(fom.trap_frequency * (t[k] - t[:k + 1]))
        naive[k] = 0.5 * step * np.sum(integrand[:-1] + integrand[1:]) / fom.trap_frequency
    fast = retarded_response(b, step, fom.trap_frequency,
                             np.cos(fom.trap_frequency * t), np.sin(fom.trap_frequency * t))
    return float(np.abs(naive - fast).max())


def reachable_mean(fom: maggradPulse) -> float:
    """Upper bound on the achievable window mean, ignoring the variance term.

    The mean is linear in B, so with a box constraint the maximizer is bang-bang: B =
    b_max * sign(d mean / d B). That derivative is the window-averaged retarded kernel.
    """
    weights = np.zeros(fom._t.size)
    step = fom.dt / 2
    for k in np.where(fom._in_window)[0]:
        seg = np.sin(fom.trap_frequency * (fom._t[k] - fom._t[:k + 1])) / fom.trap_frequency
        trapezoid = np.full(k + 1, step)
        trapezoid[0] = trapezoid[k] = step / 2
        weights[:k + 1] += seg * trapezoid
    bang_bang = b_max * np.sign(weights[:fom.drive_bins])
    return float(fom.spacings(bang_bang)[fom._in_window].mean())


def report(fom: maggradPulse, drive_pulse: np.ndarray, label: str) -> np.ndarray:
    a = fom.spacings(drive_pulse)
    in_window = a[fom._in_window]
    print(f"{label}: mean(a) over last {fom.window} = {in_window.mean():.4f} "
          f"(target {fom.a_target}, distance {abs(in_window.mean() - fom.a_target):.4f}), "
          f"std = {in_window.std():.4f}, cost = {fom.cost(drive_pulse):.4f}")
    print(f"    a over the whole run: min {a.min():.4f}, max {a.max():.4f}; "
          f"|B| max {np.abs(drive_pulse).max():.2f}")
    if a.min() <= 0.0:
        print(f"    a(t) passes through zero (min {a.min():.4f}); U ~ 1/a^3 diverges at the "
              f"crossing, so a time evolution of this pulse needs the RK4 step checked there.")
    elif a.min() < fom.a0:
        tighter = float(jl.get_critical_dt_at_spacing(a.min(), fom.tmax,
                                                      fom._setup.lattice_params,
                                                      fom._setup.hamilt_params))
        print(f"    NOTE: a dips below a0, so U ~ 1/a^3 exceeds the t=0 scale the grid was "
              f"built for. A time evolution of this pulse needs dt <= {tighter:.2e} "
              f"(grid was built at dt = {fom.dt:.2e}); resample or set min_spacing_penalty.")
    return a


def plot_result(fom: maggradPulse, drive_pulse: np.ndarray, a: np.ndarray, cost: float,
                results_path: str) -> str:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    b_full = np.zeros(fom._t.size)
    b_full[:drive_pulse.size] = drive_pulse

    fig, (ax_b, ax_a) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    ax_b.plot(fom._t, b_full, c="b")
    ax_b.axhline(b_max, ls=":", c="k")
    ax_b.axhline(-b_max, ls=":", c="k")
    if fom.drive_time < fom.tmax:
        ax_b.axvline(fom.drive_time, ls="--", c="g", label="gradient switched off")
        ax_b.legend()
    ax_b.set_ylabel("Gradient B")
    ax_b.grid()

    ax_a.plot(fom._t, a, c="r")
    ax_a.axhline(fom.a_target, ls="--", c="k", label=f"target {fom.a_target}")
    ax_a.axvspan(fom.tmax - fom.window, fom.tmax, color="orange", alpha=0.25,
                 label=f"averaging window ({fom.window})")
    ax_a.set_xlabel("Time")
    ax_a.set_ylabel("Spacing a")
    ax_a.legend()
    ax_a.grid()

    fig.suptitle(f"Optimized magnetic gradient pulse, cost {cost:.4f} "
                 f"(mean {a[fom._in_window].mean():.3f}, std {a[fom._in_window].std():.3f})")
    os.makedirs(FIGURES_DIR, exist_ok=True)
    # named after the run's QuOCS_Results folder so every run keeps its own figure rather
    # than overwriting the last one
    figure_path = os.path.join(FIGURES_DIR,
                               f"maggradPulse_{os.path.basename(results_path)}.png")
    fig.savefig(figure_path, dpi=120)
    print(f"Saved {figure_path}")
    return figure_path


def save_for_julia(fom: maggradPulse, drive_pulse: np.ndarray, results_path: str) -> str:
    """Write the pulse zero-padded onto the full tmax half-step grid, which is the array
    tevo_params_maggrad["magnetic_gradient_time"] wants as the third tuple entry."""
    b_full = np.zeros(fom._t.size)
    b_full[:drive_pulse.size] = drive_pulse
    path = os.path.join(results_path, "maggradPulse_for_julia.npz")
    np.savez(path, bpulse=b_full, bgrid=fom._t, spacings=fom.spacings(drive_pulse),
             dt=fom.dt, tmax=fom.tmax, drive_time=fom.drive_time,
             trap_frequency=fom.trap_frequency, a0=fom.a0, a_target=fom.a_target)
    print(f"Saved {path} (bpulse has {b_full.size} entries on the dt/2 grid, dt = {fom.dt})")
    return path


def main():
    fom = maggradPulse(dict(fom_params))
    print(f"RK4 step from the t=0 ED: dt = {fom.dt} -> {fom.drive_bins} pulse bins over the "
          f"{fom.drive_time} drive, {fom._t.size} samples over tmax = {fom.tmax}")
    print(f"Retarded-response self-check vs the O(n^2) rule: "
          f"max deviation {check_response_implementation(fom):.2e}")

    steady_b = fom.trap_frequency ** 2 * (fom.a_target - fom.a0)
    steady_cap = fom.a0 + b_max / fom.trap_frequency ** 2
    print(f"Steady-spacing check: holding a = {fom.a_target} without ringing needs a standing "
          f"gradient B = w^2*(a - a0) = {steady_b:.1f}, vs b_max = {b_max}. "
          f"{'OK' if steady_b <= b_max else 'OUT OF REACH'} -- with this b_max and w the "
          f"largest spacing that can be held steady at all is {steady_cap:.3f} "
          f"(b_max {fom.trap_frequency ** 2 * (fom.a_target - fom.a0):.0f} or "
          f"w {(b_max / (fom.a_target - fom.a0)) ** 0.5:.2f} would make {fom.a_target} settle).")

    ceiling = reachable_mean(fom)
    print(f"Bang-bang ceiling on the window mean at |B| <= {b_max}: {ceiling:.4f} "
          f"(target {fom.a_target} is {'reachable' if fom.a_target <= ceiling else 'NOT reachable'}; "
          f"that ceiling ignores the variance term, which pulls the mean back toward a0)")

    optimization_obj = run_optimization(build_optimization_dictionary(fom.dt), fom)

    best_controls = load_best_controls(optimization_obj)
    _, drive_pulse = best_pulse(best_controls, PULSE_NAME, TIME_NAME)
    a = report(fom, drive_pulse, "Optimized pulse")
    plot_result(fom, drive_pulse, a, fom.cost(drive_pulse), optimization_obj.results_path)
    save_for_julia(fom, drive_pulse, optimization_obj.results_path)


if __name__ == "__main__":
    main()
