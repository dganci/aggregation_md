#!/usr/bin/env python3
"""metad_quality: known populations, known surfaces, known counts."""
import numpy as np

from analysis_common import check, close, finish  # noqa: F401
import metad_quality as mdq  # noqa: E402
from constants import KB_KJ_PER_MOL_K  # noqa: E402

KT_310 = KB_KJ_PER_MOL_K * 310.0


heights = np.concatenate([np.full(50, 0.8), np.full(50, 0.08)])
decay = mdq.hill_decay(np.arange(100.0) * 500.0, heights)
check("the hill decay ratio is final over initial", close(decay["final_over_initial"], 0.1))
check("the half-height time is when it first halved", close(decay["half_height_time_ps"], 25000.0))
check("too few hills is an error, not a number", "error" in mdq.hill_decay([0.0], [1.0]))

check("equal weights give the full sample size",
      close(mdq.effective_sample_size(np.zeros(100))["ess_fraction"], 1.0))
dominant = np.concatenate([[50.0], np.zeros(999)])
ess = mdq.effective_sample_size(dominant)
check("one frame carrying the weight leaves an ESS of one", ess["ess"] < 1.01)
check("the dominant weight fraction is reported", ess["max_weight_fraction"] > 0.99)
check("log weights are never exponentiated before subtracting the maximum",
      np.isfinite(mdq.effective_sample_size(np.full(10, 5000.0))["ess"]))

check("blocks never straddle a segment boundary",
      mdq.block_ranges(20, 2, [(0, 10), (10, 20)]) == [(0, 5), (5, 10), (10, 15), (15, 20)])
check("without segments the whole array is one stretch",
      mdq.block_ranges(10, 2) == [(0, 5), (5, 10)])
check("a segment shorter than the block count is dropped",
      mdq.block_ranges(11, 5, [(0, 3), (3, 11)]) == [(3, 4), (4, 6), (6, 7), (7, 9), (9, 11)])

values = np.concatenate([np.full(400, 0.25), np.full(200, 0.75)])
values = values[np.argsort(np.tile(np.arange(200), 3), kind="stable")]
centres, fes, sigma = mdq.block_free_energy(values, np.zeros(len(values)), 310.0,
                                            bins=2, n_blocks=5, value_range=(0.0, 1.0))
check("the reweighted free energy recovers a known 2:1 population",
      close(fes[1] - fes[0], KT_310 * np.log(2.0), 1e-9))
check("the lower-population bin is the higher free energy", fes[0] == 0.0 and fes[1] > 0)
check("an exactly reproducible split has no block-to-block error",
      close(sigma[0], 0.0, 1e-9) and close(sigma[1], 0.0, 1e-9))

steps, delta = mdq.running_free_energy_difference(values, np.zeros(len(values)), 0.5, 310.0,
                                                  n_points=5)
check("the running difference converges to the same kT ln 2",
      close(delta[-1], -KT_310 * np.log(2.0), 1e-9))
check("the running difference is reported against frames used", steps[-1] == len(values))

flat = np.linspace(0.0, 100.0, 40)
check("identical surfaces have zero RMSD",
      close(mdq.fes_rmsd([flat, flat.copy()])["rmsd_kJ_per_mol"][0], 0.0))
check("a constant offset is not a change in the surface",
      close(mdq.fes_rmsd([flat + 17.0, flat])["rmsd_kJ_per_mol"][0], 0.0, 1e-9))
moved = flat.copy()
moved[-1] += 500.0
check("a change in unexplored grid does not count as drift",
      close(mdq.fes_rmsd([moved, flat], window_kJ=30.0)["rmsd_kJ_per_mol"][0], 0.0, 1e-9))
shifted = flat.copy()
shifted[0] += 4.0
check("a change where the surface is populated does count",
      mdq.fes_rmsd([shifted, flat], window_kJ=30.0)["rmsd_kJ_per_mol"][0] > 0.1)
check("a single surface cannot be compared with anything",
      "error" in mdq.fes_rmsd([flat]))

bias = np.array([1.0, 3.0, 6.0, 10.0])
rbias = np.array([1.0, 2.0, 3.0, 4.0])
ct = mdq.ct_curve(bias, rbias)
check("c(t) is bias minus rbias", close(ct["curve"], [0.0, 1.0, 3.0, 6.0]))
check("a well-ordered c(t) never goes backwards", close(ct["ct_decreasing_fraction"], 0.0))
check("a c(t) that goes backwards is caught",
      mdq.ct_curve(np.array([1.0, 3.0, 2.0]), np.zeros(3))["ct_decreasing_fraction"] > 0)

balance = mdq.walker_balance({"walker0": 100, "walker1": 90, "walker2": 0})
check("a walker that deposited nothing is named", balance["empty_walkers"] == ["walker2"])
check("walker balance is the smallest count over the largest", close(balance["min_over_max"], 0.0))
check("the hills are totalled across walkers", balance["total_hills"] == 190)


finish()
