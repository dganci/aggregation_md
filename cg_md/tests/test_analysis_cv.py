#!/usr/bin/env python3
"""cv_quality against hand-made rows and series."""
import numpy as np

from analysis_common import check, close, finish  # noqa: F401
import cv_quality as cvq  # noqa: E402


check("tied ranks are averaged", close(cvq.rank_data([3.0, 1.0, 1.0, 2.0]), [3.0, 0.5, 0.5, 2.0]))
check("distinct values rank in order", close(cvq.rank_data([5.0, 7.0, 6.0]), [0.0, 2.0, 1.0]))

x = np.arange(1.0, 21.0)
corr = cvq.correlation(x, x ** 3)
check("Spearman is 1 for any increasing relation", close(corr["spearman"], 1.0, 1e-9))
check("Pearson is below 1 for a non-linear one", corr["pearson"] < 0.98)
check("a constant series gives no correlation",
      bool(np.isnan(cvq.correlation(x, np.ones_like(x))["pearson"])))

produced = np.array([[0.0], [0.5], [1.5], [-1.0]])
cover = cvq.training_coverage(produced, train_min=[0.0], train_max=[1.0], train_std=[0.5],
                              names=["mycv.node-0"])
check("coverage counts the frames outside the training range",
      close(cover[0]["fraction_outside"], 0.5))
check("coverage measures the excursion in training standard deviations",
      close(cover[0]["max_excess_in_train_sd"], 1.0 / 0.5))
check("a frame exactly on the training edge is inside",
      close(cvq.training_coverage(np.array([[0.0], [1.0]]), [0.0], [1.0], [0.5])[0]
            ["fraction_outside"], 0.0))

grid = cvq.grid_coverage(produced, [-2.0], [2.0], names=["mycv.node-0"])
check("nothing outside the grid when the grid is wide enough",
      close(grid[0]["fraction_outside_grid"], 0.0))
check("grid occupancy is the sampled span over the grid span",
      close(grid[0]["fraction_of_grid_used"], 2.5 / 4.0))

square = np.array([0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0])
trans = cvq.basin_transitions(square, lo=0.2, hi=0.8)
check("crossings are counted between the two basins", trans["crossings"] == 3)
check("two crossings make one round trip", trans["round_trips"] == 1)
check("rattling in the middle is not a crossing",
      cvq.basin_transitions(np.array([0.0, 0.5, 0.45, 0.55, 0.5, 0.0]),
                            lo=0.2, hi=0.8)["crossings"] == 0)

independent = np.column_stack([np.sin(np.arange(500)), np.cos(np.arange(500))])
check("orthogonal components report a small off-diagonal",
      cvq.component_independence(independent)["max_abs_offdiagonal"] < 0.2)
duplicated = np.column_stack([x, 2 * x])
check("a duplicated component is caught",
      close(cvq.component_independence(duplicated)["max_abs_offdiagonal"], 1.0, 1e-9))

its_rows = [
    {"lag_frames": 5.0, "lag_ps": 50.0, "its1_ps": 900.0, "in_plateau": 0.0, "selected": 0.0},
    {"lag_frames": 50.0, "lag_ps": 500.0, "its1_ps": 500.0, "in_plateau": 1.0, "selected": 1.0},
    {"lag_frames": 100.0, "lag_ps": 1000.0, "its1_ps": 505.0, "in_plateau": 1.0, "selected": 0.0},
]
plateau = cvq.plateau_summary(its_rows, selected_lag=50)
check("the plateau is read back from the scan", plateau["plateau_lags"] == [50, 100])
check("a lag inside the plateau is recognised", plateau["selected_lag_in_plateau"])
check("a lag outside the plateau is flagged",
      not cvq.plateau_summary(its_rows, selected_lag=5)["selected_lag_in_plateau"])
check("a scan with no plateau warns",
      "warning" in cvq.plateau_summary([dict(r, in_plateau=0.0) for r in its_rows]))

lag_rows = [{"lag": 50.0, "score": 0.8, "valid": 1.0, "vamp2_train": 1.9, "vamp2_test": 1.8,
             "vamp_gap": 0.1, "ck_spectral_error": 0.05, "spectral_gap": 0.3}]
check("the lag table keeps the numbers that decide trust",
      cvq.lag_score_summary(lag_rows, 50)["selected"]["overfit_gap"] == 0.1)
check("a table where nothing passed is flagged",
      "warning" in cvq.lag_score_summary([dict(lag_rows[0], valid=0.0)]))

check("a floored hill width is reported",
      cvq.sigma_vs_spread({"sigma": [0.001, 0.5], "std": [0.0001, 1.0]})
      ["floored_components"] == [0])


fv_fields = ["time", "d_1_2", "cn_1_2", "cn_total", "rg1", "rg2", "rg_com"]
fv_data = np.array([[0.0, 3.0, 5.0, 5.0, 2.5, 2.0, 1.5],
                    [1.0, 4.0, 1.0, 1.0, 1.0, 2.0, 2.0]])
fv_names, fv = cvq.permutation_invariant_features(fv_fields, fv_data, n_prot=2)
check("feature blocks come in cg_cvgen's order",
      fv_names == ["sorted_d.1", "sorted_cn.1", "sorted_rg.1", "sorted_rg.2", "rg_com", "cn_total"])
check("a permutable block is sorted row-wise", close(fv[:, 2:4], [[2.0, 2.5], [1.0, 2.0]]))
check("the invariant scalars pass through", close(fv[:, 4:6], [[1.5, 5.0], [2.0, 1.0]]))
try:
    cvq.permutation_invariant_features(fv_fields, fv_data, n_prot=3)
    check("three protomers need three pair columns", False)
except ValueError:
    check("three protomers need three pair columns", True)


finish()
