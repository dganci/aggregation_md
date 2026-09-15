#!/usr/bin/env python3
"""cg_cvgen's Python backend: the entry point, and nothing else."""
import argparse
import json
from pathlib import Path

import numpy as np

import cvgen_deps
from cvgen_clock import rebase_clocks  # noqa: F401
from cvgen_data import Data
from cvgen_deps import deps
from cvgen_export import Exporter
from cvgen_metrics import Metrics  # noqa: F401
from cvgen_plateau import PlateauScan
from cvgen_score import Scorer
from cvgen_train import Trainer


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--input-colvar", nargs="+", default=None,
                    help="Override the config's input_colvar; accepts several paths or globs, "
                         "e.g. an adaptive run's COLVAR_chunk_*.dat series.")
    args = ap.parse_args()
    with open(args.config) as f:
        cfg = json.load(f)
    if args.input_colvar:
        cfg["input_colvar"] = args.input_colvar
    Path(cfg["output_dir"]).mkdir(parents=True, exist_ok=True)
    cvgen_deps.load()
    np.random.seed(cfg["seed"])
    deps.torch.manual_seed(cfg["seed"])

    data = Data(cfg).load()
    print(f"features={len(data.feature_cols)} {data.feature_cols}")

    if cfg.get("auto_lags", True):
        auto = PlateauScan.run(cfg, data, cfg["output_dir"])
        if auto:
            cfg["lags"] = auto

    results, outputs = {}, {}
    for lag in cfg["lags"]:
        print(f"lag={lag}")
        out = Trainer(cfg, data, lag, cfg["split_ratio"]).run()
        results[int(lag)] = out.metrics
        outputs[int(lag)] = out
    Scorer.apply(results, cfg["n_cvs"])
    selected = Scorer.choose(results, cfg.get("selected_lag", 0))
    print(f"selected_lag={selected} score={results[selected]['score']:.4f}")

    final = Trainer(cfg, data, selected, 1.0).run()
    exporter = Exporter(cfg, data)
    exporter.write_scores(results)
    params = exporter.metad_params(final.s_train)
    exporter.final_model(final.model, selected, final, params, results)
    print("written", cfg["output_dir"])


if __name__ == "__main__":
    main()
