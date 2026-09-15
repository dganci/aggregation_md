#!/usr/bin/env python3
"""analyze_offline.py end to end: the three subcommands on files written for the purpose,
including a COLVAR whose clock restarts.
"""
import json
import tempfile
from pathlib import Path

import numpy as np

from analysis_common import check, close, finish, run_cli, write_colvar  # noqa: F401


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    walker = root / "walker0"
    walker.mkdir()
    out = root / "report"

    n_frames = 400
    times = np.arange(n_frames, dtype=float) * 10.0
    cn12 = np.where(times >= 2000.0, 1.4, 0.2)
    cn13 = np.full(n_frames, 0.1)
    cn23 = np.full(n_frames, 0.1)
    rows = np.column_stack([times, cn12, cn13, cn23, cn12 + cn13 + cn23,
                            3.0 + 0.1 * np.sin(times / 200.0)])
    write_colvar(root / "COLVAR", ["time", "cn_1_2", "cn_1_3", "cn_2_3", "cn_total", "rg_com"],
                 rows)

    result = run_cli("aggregation", "--colvar", str(root / "COLVAR"), "--n-prot", "3",
                     "--out-dir", str(out))
    check("the aggregation subcommand exits cleanly", result.returncode == 0)
    summary = json.loads((out / "aggregation_summary.json").read_text()) \
        if (out / "aggregation_summary.json").is_file() else {}
    check("aggregation reports every frame", summary.get("frames") == n_frames)
    check("aggregation finds the dimer", summary.get("clustering", {}).get("largest_max") == 2)
    check("aggregation times the association",
          summary.get("first_passage_ps", {}).get("2") == 2000.0)
    check("aggregation does not invent a trimer",
          summary.get("first_passage_ps", {}).get("3") is None)

    cv0 = np.sin(times / 300.0)
    cv1 = np.cos(times / 700.0)
    rbias = 2.0 + 0.5 * np.sin(times / 150.0)
    write_colvar(walker / "COLVAR",
                 ["time", "mycv.node-0", "mycv.node-1", "metad.bias", "metad.rbias"],
                 np.column_stack([times, cv0, cv1, rbias + np.linspace(0, 5, n_frames), rbias]))
    write_colvar(walker / "COLVAR_monitor",
                 ["time", "cn_total", "rg_com", "metad.rbias"],
                 np.column_stack([times, cn12 + cn13 + cn23, rows[:, 5], rbias]))
    hill_times = np.arange(40, dtype=float) * 100.0
    write_colvar(walker / "HILLS",
                 ["time", "mycv.node-0", "mycv.node-1", "sigma_mycv.node-0",
                  "sigma_mycv.node-1", "height", "biasf"],
                 np.column_stack([hill_times, cv0[:40], cv1[:40], np.full(40, 0.1),
                                  np.full(40, 0.1), 0.8 * np.exp(-hill_times / 1500.0),
                                  np.full(40, 10.0)]))

    cv_dir = root / "cvs"
    cv_dir.mkdir()
    (cv_dir / "cv_params.json").write_text(json.dumps({
        "n_cvs": 2, "min": [-1.0, -1.0], "max": [1.0, 1.0], "std": [0.7, 0.7],
        "grid_min": [-1.5, -1.5], "grid_max": [1.5, 1.5], "sigma": [0.35, 0.35]}))
    (cv_dir / "cv_manifest.json").write_text(json.dumps({"selected_lag": 50, "n_cvs": 2}))
    (cv_dir / "its_scan.csv").write_text(
        "lag_frames,lag_ps,its1_ps,its2_ps,in_plateau,selected\n"
        "5,50,900,300,0,0\n50,500,500,200,1,1\n100,1000,505,205,1,0\n")
    (cv_dir / "lag_scores.csv").write_text(
        "lag,score,valid,vamp2_train,vamp2_test,vamp_gap,ck_spectral_error,"
        "spectral_gap,timescale_plateau,final_rank\n"
        "50,0.8,True,1.9,1.8,0.1,0.05,0.3,0.02,2\n")

    result = run_cli("cv", "--cv-dir", str(cv_dir), "--colvar", str(walker / "COLVAR"),
                     "--monitor", str(walker / "COLVAR_monitor"), "--out-dir", str(out))
    check("the cv subcommand exits cleanly", result.returncode == 0)
    cv_summary = json.loads((out / "cv_quality.json").read_text()) \
        if (out / "cv_quality.json").is_file() else {}
    check("cv reads both CV components", len(cv_summary.get("training_coverage", [])) == 2)
    check("cv confirms the run stayed inside the grid",
          all(r["fraction_outside_grid"] == 0 for r in cv_summary.get("grid_coverage", [{}])
              if "fraction_outside_grid" in r))
    check("cv says the selected lag was in the plateau",
          cv_summary.get("lag_plateau", {}).get("selected_lag_in_plateau") is True)
    check("cv correlates the components against the observables",
          len(cv_summary.get("cv_interpretation", [])) == 4)

    result = run_cli("metad", "--hills", str(walker / "HILLS"),
                     "--colvar", str(walker / "COLVAR"),
                     "--monitor", str(walker / "COLVAR_monitor"), "--out-dir", str(out))
    check("the metad subcommand exits cleanly", result.returncode == 0)
    md_summary = json.loads((out / "metad_quality.json").read_text()) \
        if (out / "metad_quality.json").is_file() else {}
    check("metad sees the hills decay",
          md_summary.get("hill_decay", {}).get("walker0", {}).get("final_over_initial", 1.0) < 0.2)
    check("metad reports the effective sample size",
          0 < md_summary.get("reweighting", {}).get("ess_fraction", 0) <= 1.0)
    check("metad builds a reweighted surface per observable",
          set(md_summary.get("reweighted_fes", {})) == {"cn_total", "rg_com"})
    check("metad writes JSON a strict parser accepts",
          "NaN" not in (out / "metad_quality.json").read_text())

    restarted = np.vstack([rows, np.column_stack([times, cn12, cn13, cn23,
                                                  cn12 + cn13 + cn23, rows[:, 5]])])
    write_colvar(root / "COLVAR_batched",
                 ["time", "cn_1_2", "cn_1_3", "cn_2_3", "cn_total", "rg_com"], restarted)
    result = run_cli("aggregation", "--colvar", str(root / "COLVAR_batched"), "--n-prot", "3",
                     "--out-dir", str(root / "report2"))
    batched = json.loads((root / "report2" / "aggregation_summary.json").read_text()) \
        if (root / "report2" / "aggregation_summary.json").is_file() else {}
    check("a restarted clock keeps both batches", batched.get("frames") == 2 * n_frames)
    check("a restarted clock doubles the sampled time",
          close(batched.get("sampled_time_us", 0.0), 2 * summary.get("sampled_time_us", 0.0), 1e-9))


finish()
