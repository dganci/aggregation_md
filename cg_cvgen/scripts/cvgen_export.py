"""Writing the model and its parameters in the form cg_md reads back."""
import csv
import json
import math
import pickle
from pathlib import Path

import numpy as np

from cvgen_deps import deps


class Exporter:
    def __init__(self, cfg, data):
        self.cfg = cfg
        self.data = data
        self.out = Path(cfg["output_dir"])
        self.out.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def clean(obj):
        if isinstance(obj, dict):
            return {str(k): Exporter.clean(v) for k, v in obj.items() if not str(k).startswith("_")}
        if isinstance(obj, (list, tuple)):
            return [Exporter.clean(v) for v in obj]
        if isinstance(obj, np.ndarray):
            return Exporter.clean(obj.tolist())
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            obj = float(obj)
        if isinstance(obj, float):
            return obj if math.isfinite(obj) else None
        return obj

    def write_scores(self, results):
        fields = ["lag", "score", "valid", "vamp2_train", "vamp2_test", "vamp_gap", "ck_spectral_error", "spectral_gap", "timescale_plateau", "final_rank"]
        with open(self.out / "lag_scores.csv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields)
            w.writeheader()
            for lag in sorted(results):
                row = {"lag": lag, **{k: results[lag].get(k) for k in fields if k != "lag"}}
                w.writerow(row)
        with open(self.out / "lag_scores.json", "w") as f:
            json.dump(Exporter.clean(results), f, indent=2)

    def metad_params(self, s):
        mins = s.min(axis=0)
        maxs = s.max(axis=0)
        stds = s.std(axis=0)
        span = maxs - mins
        margin = self.cfg["margin"]

        sigma = 0.5 * stds
        floor = np.maximum(1e-3, np.where(span > 0, span / 1000.0, 1e-3))
        sigma = np.round(np.maximum(sigma, floor), 3)
        degenerate = [i for i, v in enumerate(0.5 * stds) if v < floor[i]]
        if degenerate:
            print(f"WARNING: CV component(s) {degenerate} have a near-zero spread; "
                  f"sigma floored to {sigma[degenerate]}. The CV is likely degenerate - "
                  f"check lag_scores.csv before biasing with it.")

        grid_min = np.round(mins - margin * np.where(span > 0, span, 1.0), 3)
        grid_max = np.round(maxs + margin * np.where(span > 0, span, 1.0), 3)

        params = {
            "n_cvs": int(s.shape[1]),
            "min": np.round(mins, 3),
            "max": np.round(maxs, 3),
            "std": np.round(stds, 3),
            "grid_min": grid_min,
            "grid_max": grid_max,
            "sigma": sigma,
            "feature_cols": self.data.feature_cols,
            "permutation_invariant": bool(self.cfg.get("permutation_invariant", True)),
        }
        with open(self.out / "cv_params.pkl", "wb") as f:
            pickle.dump(params, f)
        with open(self.out / "cv_params.json", "w") as f:
            json.dump(Exporter.clean(params), f, indent=2)
        return params

    def torchscript(self, model, x_ref):
        import torch
        from torch import nn

        children = list(model.children())
        if len(children) < 4:
            raise RuntimeError("Unsupported DeepTICA layout; cannot export plain TorchScript model")
        norm, ff, tica = children[1], children[2], children[3]

        class PlainNorm(nn.Module):
            def __init__(self, mean, range_):
                super().__init__()
                self.register_buffer("mean", mean.clone().detach())
                self.register_buffer("range_", range_.clone().detach())
            def forward(self, x):
                return (x - self.mean) / self.range_

        class PlainTICA(nn.Module):
            def __init__(self, mean, evecs):
                super().__init__()
                self.register_buffer("mean", mean.clone().detach())
                self.register_buffer("evecs", evecs.clone().detach())
            def forward(self, x):
                return torch.matmul(x - self.mean, self.evecs)

        class CVModel(nn.Module):
            def __init__(self, net):
                super().__init__()
                self.net = net
            def forward(self, x):
                return self.net(x)

        net = nn.Sequential(PlainNorm(norm.mean, norm.range), ff.nn, PlainTICA(tica.mean, tica.evecs)).eval()
        with torch.no_grad():
            ref = model(x_ref[: min(16, len(x_ref))]).detach().cpu()
            got = net(x_ref[: min(16, len(x_ref))]).detach().cpu()
            if not torch.allclose(ref, got, atol=1e-5):
                raise RuntimeError("Plain TorchScript export does not match DeepTICA output")
        scripted = torch.jit.script(CVModel(net).eval())
        scripted.save(str(self.out / "CVs_torchscript.pt"))

    def final_model(self, model, selected_lag, train_output, params, results):
        checkpoint = {
            "state_dict": model.state_dict(),
            "n_cvs": self.cfg["n_cvs"],
            "input_dim": len(self.data.feature_cols),
            "feature_cols": self.data.feature_cols,
            "selected_lag": int(selected_lag),
        }
        deps.torch.save(checkpoint, self.out / f"CVs_model_lag{selected_lag}.pt")
        with open(self.out / f"model_lag{selected_lag}.pkl", "wb") as f:
            pickle.dump(model, f)
        self.torchscript(model, train_output.x_train)
        if self.cfg.get("save_embeddings"):
            with open(self.out / f"CVs_emb_lag{selected_lag}.pkl", "wb") as f:
                pickle.dump(train_output.s_train, f)
        manifest = {
            "model": "CVs_torchscript.pt",
            "cv_params_pkl": "cv_params.pkl",
            "cv_params_json": "cv_params.json",
            "feature_schema": "feature_schema.json",
            "selected_lag": int(selected_lag),
            "n_cvs": self.cfg["n_cvs"],
            "feature_cols": self.data.feature_cols,
            "sigma": Exporter.clean(params["sigma"]),
            "grid_min": Exporter.clean(params["grid_min"]),
            "grid_max": Exporter.clean(params["grid_max"]),
            "lag_scores": "lag_scores.csv",
        }
        with open(self.out / "cv_manifest.json", "w") as f:
            json.dump(manifest, f, indent=2)
        with open(self.out / "feature_schema.json", "w") as f:
            json.dump({"feature_cols": self.data.feature_cols, "n_features": len(self.data.feature_cols)}, f, indent=2)
