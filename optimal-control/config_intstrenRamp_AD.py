"""QuOCS AD (automatic differentiation) optimization of the interaction strength ramp.

Gradient-based counterpart of config_intstrenRamp.py: same figure of merit and pulse
parametrization, but optimized with QuOCS' AD algorithm (L-BFGS-B on the Fourier
coefficients with exact JAX gradients) instead of dCRAB's Nelder-Mead.

QuOCS' ADAlgorithm pushes jax.grad/jax.jit straight through get_FoM, so the FoM must
be pure JAX and cannot call Julia. Julia therefore computes only the pulse-independent
constants, once at construction (setup_intstren_ramp_ad in
exact-diag/intstren-ramp-ad-functions.jl), and the propagation is reimplemented in
jnp below.
"""

import math

import numpy as np

from quocs_common import (JuliaFoM, fourier_pulse, include_julia, jl,
                          linear_ramp_lambda, load_best_controls, plot_best_pulses,
                          run_optimization)

include_julia("intstren-ramp-ad-functions.jl")

# quocslib's ADAlgorithm also enables x64 at import time; do it before building the
# constant arrays so they are complex128 from the start
import jax
jax.config.update("jax_enable_x64", True)
import jax.numpy as jnp

# above this Hilbert space dimension the hopping matrix is kept sparse (BCOO);
# below it a dense matmul is faster and the memory (<= dim^2 * 16 B) is negligible
DENSE_HOPPING_DIM_LIMIT = 1024


class intstrenRampAD(JuliaFoM):
    """
    Same figure of merit as intstrenRamp in config_intstrenRamp.py, but written in pure
    JAX so QuOCS' AD algorithm can push jax.grad/jax.jit straight through get_FoM.
    The pulse-independent constants (endpoint ED manifolds and the H(u) = H_hop + u*H_int
    decomposition, taken from the saved undressed matrices via getHopping/getInteraction)
    are computed once in Julia at construction; the RK4 propagation and the ground-state
    manifold fidelity are reimplemented here in jnp, mirroring runge_kutta_step /
    run_ramp_stages (time-evolution.jl) and groundstate_manifold_fidelity
    (control-functions.jl). The interaction is diagonal, so it is applied elementwise
    from its diagonal; the hopping is applied dense below DENSE_HOPPING_DIM_LIMIT and
    as a sparse BCOO matvec above, keeping memory at O(nnz) for large Hilbert spaces.
    """

    def __init__(self, args_dict: dict = None):
        if args_dict is None:
            args_dict = {}

        # hamiltonian / lattice parameters
        self.Lx = args_dict.setdefault("Lx", 4)
        self.Ly = args_dict.setdefault("Ly", 4)
        self.N = args_dict.setdefault("N", 2)
        self.lr = args_dict.setdefault("lr", "all")
        self.if_periodic_x = args_dict.setdefault("if_periodic_x", True)
        self.if_periodic_y = args_dict.setdefault("if_periodic_y", True)

        # interaction profile over distance: "flat" (default), "exp", "gaussian",
        # "rydberg" or "dd", with their respective length/width parameters (defaults
        # mirror get_normal_model_params_ed); all profiles are linear in the overall
        # strength, so the H_hop + u*H_int split holds for each of them
        self.scaling_type = args_dict.setdefault("scaling_type", "flat")
        self.corr_length = args_dict.setdefault("corr_length", float(self.Ly))
        self.sigma = args_dict.setdefault("sigma", 1.0)
        self.blockade_radius = args_dict.setdefault("blockade_radius", 1.0)
        self.magnetic_spacing = args_dict.setdefault("magnetic_spacing", 1.0)

        # ramp targets. The ramp duration is fixed (only the pulse shape is optimized)
        # and must match the ramptime value the pulse dictionary's bins are derived from
        self.intstren_start = args_dict.setdefault("intstren_start", 10.0)
        self.intstren_end = args_dict.setdefault("intstren_end", 0.0)
        self.ramptime = args_dict.setdefault("ramptime", 1.0)
        # largest strength the pulse can reach -- must match the pulse dictionary's
        # upper_limit; H_int is dressed at this value in Julia so that no decaying-profile
        # coupling above interaction_cutoff anywhere along the ramp is missing from it
        self.intstren_max = args_dict.setdefault("intstren_max", 12.0)

        # figure of merit / time evolution settings
        self.speccount = args_dict.setdefault("speccount", 2)
        self.dt = args_dict.setdefault("dt", 0.005)

        # pulse_ramp's half-step grid (spacing dt/2): 2*n_steps + 1 samples, RK4 step m
        # uses samples 2m, 2m+1, 2m+2 as H(t), H(t + dt/2), H(t + dt)
        self.n_steps = math.ceil(self.ramptime / self.dt)
        self.n_bins = 2 * self.n_steps + 1

        (hop_rows, hop_cols, hop_vals, hint_diag, dim,
         starting_states, target_states, max_residual) = jl.setup_intstren_ramp_ad(self.to_julia_dict())

        # the split is checked in Julia against the ED eigenpairs of both endpoint
        # Hamiltonians, so this bounds any mismatch with what run_normal_ed diagonalized.
        # long_range_scaling rounds U to 5 digits at every strength, so for non-flat
        # profiles u*H_int legitimately deviates from the rounded ED couplings by up to
        # ~1e-5 per interacting pair; the tolerance sits above that but well below the
        # ~1e-4 signature of a coupling dropped by interaction_cutoff
        max_residual = float(max_residual)
        residual_tol = 1e-5 * max(1, self.N * (self.N - 1))
        assert max_residual < residual_tol, \
            f"H(u) = H_hop + u*H_int decomposition failed: max eigenstate residual {max_residual}"
        self.dim = int(dim)
        print(f"Julia setup done. dim = {self.dim}, max endpoint eigenstate residual {max_residual:.2e}")

        # Julia COO indices are 1-based
        hop_rows = np.asarray(hop_rows, dtype=np.int64) - 1
        hop_cols = np.asarray(hop_cols, dtype=np.int64) - 1
        hop_vals = np.asarray(hop_vals, dtype=np.complex128)
        if self.dim <= DENSE_HOPPING_DIM_LIMIT:
            ham_hop = np.zeros((self.dim, self.dim), dtype=np.complex128)
            ham_hop[hop_rows, hop_cols] = hop_vals
            self._ham_hop = jnp.asarray(ham_hop)
        else:
            from jax.experimental import sparse as jsparse
            self._ham_hop = jsparse.BCOO(
                (jnp.asarray(hop_vals), jnp.asarray(np.stack([hop_rows, hop_cols], axis=1))),
                shape=(self.dim, self.dim))
        # the interaction matrix is diagonal, so keep only its diagonal and apply it
        # elementwise (column shape broadcasts over the states)
        self._hint_diag = jnp.asarray(np.asarray(hint_diag), dtype=jnp.complex128)[:, None]
        self._starting_states = jnp.asarray(np.asarray(starting_states), dtype=jnp.complex128)
        self._target_states = jnp.asarray(np.asarray(target_states), dtype=jnp.complex128)
        self._fidelity = jax.jit(self._make_fidelity())

    def _make_fidelity(self):
        ham_hop, hint_diag = self._ham_hop, self._hint_diag
        starting_states, target_states = self._starting_states, self._target_states
        dt = self.dt

        # @ works for both the dense and the BCOO hopping representation
        def apply_ham(u, psi):
            return ham_hop @ psi + u * (hint_diag * psi)

        # mirrors runge_kutta_step (time-evolution.jl), including the per-step
        # renormalization time_evolution applies after each RK4 step
        def rk4_step(psi, u3):
            k1 = -1j * apply_ham(u3[0], psi)
            k2 = -1j * apply_ham(u3[1], psi + (0.5 * dt) * k1)
            k3 = -1j * apply_ham(u3[1], psi + (0.5 * dt) * k2)
            k4 = -1j * apply_ham(u3[2], psi + dt * k3)
            psi_new = psi + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
            psi_new = psi_new / jnp.linalg.norm(psi_new, axis=0, keepdims=True)
            return psi_new, None

        def fidelity(pulse):
            # group the half-step samples into (n_steps, 3) triples (t, t+dt/2, t+dt)
            u_triples = jnp.stack([pulse[0:-2:2], pulse[1:-1:2], pulse[2::2]], axis=1)
            final_states, _ = jax.lax.scan(rk4_step, starting_states, u_triples)
            # groundstate_manifold_fidelity: 0.5 * tr(F^dag F), F_ij = <psi_i|target_j>
            overlap_matrix = final_states.conj().T @ target_states
            return 0.5 * jnp.sum(jnp.abs(overlap_matrix) ** 2)

        return fidelity

    def get_FoM(self, pulses: list = [], parameters: list = [], timegrids: list = []) -> dict:
        # in AD mode Controls hands the pulses over as rows of a complex64 jnp array
        # (possibly as jax tracers during jax.grad) -- stay in jnp throughout
        pulse = jnp.real(jnp.asarray(pulses[0])[:self.n_bins])
        return {"FoM": self._fidelity(pulse)}


# must match intstrenRampAD's intstren_start / intstren_end / ramptime / dt defaults
# above; intstren_max is both the pulse's upper amplitude limit and the strength H_int
# is dressed at in the Julia setup
intstren_start = 10.0
intstren_end = 0.0
intstren_max = 12.0
ramptime = 1.0
dt = 0.005

optimization_dictionary = {
    "optimization_client_name": "intstrenRamp_AD",
    # each super iteration redraws the random Fourier frequencies (as in dCRAB) and then
    # runs L-BFGS-B on the basis coefficients with jax gradients instead of NelderMead.
    # ADAlgorithm reads its L-BFGS-B stopping settings from stopping_criteria (not from
    # dsm_settings, which is dCRAB-only); max_eval_total is read from both levels
    "algorithm_settings": {
        "algorithm_name": "AD",
        "optimization_direction": "maximization",
        "super_iteration_number": 5,
        "max_eval_total": 200,
        "stopping_criteria": {
            "max_eval_total": 200,
            "ftol": 1e-10,
            "gtol": 1e-8,
        },
    },
    "pulses": [
        fourier_pulse(pulse_name="intstrenRamp",
                      time_name="time_intstrenRamp",
                      ramptime=ramptime,
                      dt=dt,
                      lower_limit=0.0,
                      upper_limit=intstren_max,
                      amplitude_variation=3.0,
                      initial_guess_lambda=linear_ramp_lambda(intstren_start, intstren_end)),
    ],
    "parameters": [],
    "times": [{"time_name": "time_intstrenRamp", "initial_value": ramptime}],
}


def main():
    optimization_obj = run_optimization(optimization_dictionary, intstrenRampAD({
        "intstren_start": intstren_start,
        "intstren_end": intstren_end,
        "intstren_max": intstren_max,
        "ramptime": ramptime,
        "dt": dt,
    }))

    best_controls = load_best_controls(optimization_obj)
    fidelity = abs(optimization_obj.opt_alg_obj.best_FoM)
    plot_best_pulses(best_controls,
                     [("intstrenRamp", "time_intstrenRamp", "Interaction strength")],
                     title=f"AD-Optimized Intstren-Ramp Fidelity: {fidelity:.4f}",
                     filename="intstrenRamp_AD_optimized.png")


if __name__ == "__main__":
    main()
