#!/usr/bin/env python3
"""aggregation_metrics against arithmetic done by hand."""
import numpy as np

from analysis_common import check, close, finish  # noqa: F401
import aggregation_metrics as agg  # noqa: E402
from constants import KB_KJ_PER_MOL_K  # noqa: E402

KT_310 = KB_KJ_PER_MOL_K * 310.0


fields = ["time", "cn_1_2", "cn_1_3", "cn_2_3", "cn_total"]
pairs, cols = agg.pair_columns(fields, 3)
check("pair_columns finds every pair present", pairs == [(0, 1), (0, 2), (1, 2)])
check("pair_columns returns the matching columns", cols == [1, 2, 3])
check("a missing pair is left out rather than guessed",
      agg.pair_columns(["cn_1_2", "cn_2_3"], 3)[0] == [(0, 1), (1, 2)])

contact = np.array([
    [False, False, False],
    [True, False, False],
    [True, False, True],
    [False, True, False],
])
labels = agg.cluster_labels(contact, pairs, 3)
sizes = agg.cluster_sizes(labels, 3)
check("no contacts means three monomers", sorted(sizes[0][sizes[0] > 0]) == [1, 1, 1])
check("one contact gives a dimer and a monomer", sorted(sizes[1][sizes[1] > 0]) == [1, 2])
check("a chain of contacts merges all three", sorted(sizes[2][sizes[2] > 0]) == [3])
check("contacts need not be adjacent protomers", sorted(sizes[3][sizes[3] > 0]) == [1, 2])
check("every frame accounts for every chain", bool(np.all(sizes.sum(axis=1) == 3)))
check("the largest oligomer is the row maximum", list(sizes.max(axis=1)) == [1, 2, 3, 2])

counts = agg.size_histogram(sizes, 3)
check("size histogram counts oligomers, not chains", list(counts) == [0.0, 5.0, 2.0, 1.0])

dist = agg.size_distribution(counts)
check("number-average size is chains over oligomers", close(dist["mean_size_number"], 12 / 8))
check("weight-average size is the chain-weighted one",
      close(dist["mean_size_weight"], (5 * 1 + 2 * 4 + 1 * 9) / 12))
check("the number fractions add to one", close(sum(dist["number_fraction"]), 1.0))
check("the mass fractions add to one", close(sum(dist["mass_fraction"]), 1.0))

g = agg.size_free_energy(counts, 310.0)
check("the monomer is the zero of the size free energy", close(g[0], 0.0))
check("a rarer size costs free energy", close(g[1], -KT_310 * np.log(2.0 / 5.0)))
check("a size never seen is NaN, not infinity",
      bool(np.isnan(agg.size_free_energy(np.array([0.0, 4.0, 0.0, 1.0]), 310.0)[1])))

lengths, censored = agg.contact_runs([True, True, False, True, False, False, True])
check("contact_runs measures every uninterrupted run", list(lengths) == [2, 1, 1])
check("a run starting at frame 0 is censored", bool(censored[0]))
check("a run ending at the last frame is censored", bool(censored[-1]))
check("a run inside the trajectory is not censored", not bool(censored[1]))
check("one run spanning everything is censored once, from both ends",
      list(agg.contact_runs([True, True, True])[1]) == [True])
check("a series with no contact has no runs", len(agg.contact_runs([False, False])[0]) == 0)

series = np.array([[True], [True], [False], [True], [True], [True], [False], [True]])
life = agg.contact_lifetimes(series, dt_ps=10.0)
check("the mean lifetime uses only the complete runs", close(life["mean_lifetime_ps"], 30.0))
check("the censored runs are counted and reported", life["censored_runs"] == 2)
check("the complete runs are counted", life["complete_runs"] == 1)

always = np.ones((64, 1), dtype=bool)
acf = agg.contact_autocorrelation(always, max_lag=8)
check("C(0) is 1 by construction", close(acf[0], 1.0, 1e-12))
check("a contact that never breaks has C(t) = 1", close(acf, 1.0, 1e-9))

alternating = (np.arange(64) % 2 == 0).reshape(-1, 1)
acf_alt = agg.contact_autocorrelation(alternating, max_lag=4)
check("an alternating contact is anticorrelated at lag 1", close(acf_alt[1], -1.0, 1e-9))
check("an alternating contact is back in phase at lag 2", close(acf_alt[2], 1.0, 1e-9))

noise = np.random.default_rng(7).random((4000, 6)) < 0.35
tail = agg.contact_autocorrelation(noise, max_lag=50)
check("the correlation starts at one", close(tail[0], 1.0, 1e-9))
check("independent frames decorrelate to zero, not to the occupancy",
      float(np.abs(tail[1:]).max()) < 0.1)

check("the correlation time stops at the first zero crossing",
      close(agg.correlation_time([1.0, 0.5, 0.0, 0.4], dt_ps=2.0), 2 * (0.5 + 0.5 * 0.5)))
check("a correlation that never crosses zero is integrated whole",
      close(agg.correlation_time([1.0, 0.5], dt_ps=2.0), 2 * 0.75))

times = np.array([0.0, 10.0, 20.0, 30.0])
check("first passage reports when the size was first reached",
      close(agg.first_passage([1, 2, 3, 2], 3, times), 20.0))
check("a size never reached has no first-passage time",
      agg.first_passage([1, 2, 2, 1], 4, times) is None)

sphere = agg.shape_descriptors([1.0, 1.0, 1.0])
check("a sphere has no asphericity", close(sphere["asphericity"], 0.0))
check("a sphere has zero anisotropy", close(sphere["anisotropy"], 0.0))
check("a sphere has NPMI 1", close(sphere["npmi"], 1.0))
rod = agg.shape_descriptors([1.0, 0.0, 0.0])
check("a rod has anisotropy 1", close(rod["anisotropy"], 1.0))
check("a rod has NPMI 0", close(rod["npmi"], 0.0))
oblate = agg.shape_descriptors([2.0, 2.0, 1.0])
check("NPMI is not the ratio of gyration eigenvalues",
      not close(oblate["npmi"], 1.0 / 2.0, 1e-3))
check("NPMI is I_min over I_max", close(oblate["npmi"], (5 - 2) / (5 - 1)))

square = agg.shape_descriptors([0.25, 0.25, 0.0])
check("a flat square is not reported as a rod", close(square["npmi"], 0.5))
check("a flat square has the anisotropy a flat square has",
      close(square["anisotropy"], 0.25))

cens = agg.censoring_report([5, 5, 3, 5], 5)
check("censoring reports how often the box was full", close(cens["fraction_at_n_prot"], 0.75))


finish()
