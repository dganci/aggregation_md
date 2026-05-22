#!/usr/bin/env python3
import argparse
import csv
import json
import math
import os
import pickle
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np


def load_deps():
    global pd, torch, nn, mlcolvar, L, EarlyStopping, Callback, create_timelagged_dataset
    import pandas as pd
    import torch
    from torch import nn
    import mlcolvar
    import lightning.pytorch as L
    from lightning.pytorch.callbacks import EarlyStopping, Callback
    from mlcolvar.utils.timelagged import create_timelagged_dataset


class EpochLogger:
    def __init__(self, interval=25, enabled=False):
        self.interval = interval
        self.enabled = enabled

    def callback(self):
        parent = self

        class _Logger(Callback):
            def on_train_epoch_end(self, trainer, pl_module):
                if not parent.enabled or trainer.current_epoch % parent.interval:
                    return
                tr = trainer.callback_metrics.get("train_loss")
                va = trainer.callback_metrics.get("valid_loss")
                print(f"epoch={trainer.current_epoch} train={tr} valid={va}")

        return _Logger()


class MetricsLogger:
    def __init__(self, probe_x, lag):
        self.probe_x = probe_x
        self.lag = lag
        self.history = []
        self.vamp2_history = []
        self.its_history = []

    def callback(self):
        parent = self

        class _Metrics(Callback):
            def on_train_epoch_end(self, trainer, pl_module):
                with torch.no_grad():
                    out = pl_module.forward_nn(parent.probe_x).detach().cpu().numpy()
                parent.history.append(out)
                if len(out) <= parent.lag:
                    parent.vamp2_history.append(np.nan)
                    parent.its_history.append(np.nan)
                    return
                z0, z1 = out[:-parent.lag], out[parent.lag:]
                score, svals = Metrics.vamp2_score(z0, z1, return_singular=True)
                parent.vamp2_history.append(score)
                eig = np.sort(np.real(svals))[::-1][1:]
                valid = (eig > 1e-6) & (eig < 1 - 1e-6)
                parent.its_history.append(float(np.mean(-parent.lag / np.log(eig[valid]))) if np.any(valid) else np.nan)

        return _Metrics()


class Metrics:
    @staticmethod
    def make_pairs(x, lag):
        return x[:-lag], x[lag:]

    @staticmethod
    def inv_sqrt(c, eps=1e-6):
        w, v = np.linalg.eigh(c + eps * np.eye(c.shape[0]))
        return v @ np.diag(1.0 / np.sqrt(np.maximum(w, eps))) @ v.T

    @staticmethod
    def vamp2_score(z0, z1, eps=1e-6, return_singular=False):
        z0 = np.asarray(z0) - np.asarray(z0).mean(axis=0, keepdims=True)
        z1 = np.asarray(z1) - np.asarray(z1).mean(axis=0, keepdims=True)
        n = len(z0)
        c00 = (z0.T @ z0) / n
        c11 = (z1.T @ z1) / n
        c01 = (z0.T @ z1) / n
        k = Metrics.inv_sqrt(c00, eps) @ c01 @ Metrics.inv_sqrt(c11, eps)
        s = np.linalg.svd(k, compute_uv=False)
        score = float(np.sum(s ** 2))
        return (score, s) if return_singular else score

    @staticmethod
    def koopman(z0, z1, eps=1e-6):
        z0 = np.asarray(z0)
        z1 = np.asarray(z1)
        n = len(z0)
        c00 = (z0.T @ z0) / n
        c11 = (z1.T @ z1) / n
        c01 = (z0.T @ z1) / n
        return Metrics.inv_sqrt(c00, eps) @ c01 @ Metrics.inv_sqrt(c11, eps)

    @staticmethod
    def diagnostics(s, lag, kstep=2):
        if len(s) <= lag * kstep or len(s) < 5 * lag:
            return {}
        z0, z1 = Metrics.make_pairs(s, lag)
        z0k, z1k = Metrics.make_pairs(s, lag * kstep)
        k = Metrics.koopman(z0, z1)
        kk = Metrics.koopman(z0k, z1k)
        sv = np.sort(np.linalg.svd(k, compute_uv=False))[::-1]
        svk = np.sort(np.linalg.svd(kk, compute_uv=False))[::-1]
        m = min(len(sv), len(svk))
        ck = float(np.linalg.norm((sv[:m] ** kstep) - svk[:m])) if m else np.nan
        eig = np.sort(np.real(np.linalg.eigvals(k)))[::-1][1:]
        valid = (eig > 1e-12) & (eig < 1 - 1e-12)
        times = (-lag / np.log(eig[valid])).astype(float) if np.any(valid) else np.array([])
        return {
            "singular_vals": sv.astype(float).tolist(),
            "ck_spectral_error": ck,
            "spectral_gap_0-1": float(max(1.0 - sv[1], 0.0)) if len(sv) > 1 else np.nan,
            "spectral_gap_1-2": float(max(sv[1] - sv[2], 0.0)) if len(sv) > 2 else np.nan,
            "timescales": times[np.isfinite(times)].tolist(),
        }

    @staticmethod
    def bootstrap_vamp2(s, lag, n_boot=20):
        block = max(2 * lag, 1)
        n_blocks = max(len(s) // block, 1)
        scores = []
        for _ in range(n_boot):
            starts = np.random.randint(0, max(len(s) - block + 1, 1), size=n_blocks)
            sample = np.concatenate([s[i:i + block] for i in starts], axis=0)
            if len(sample) > lag:
                z0, z1 = Metrics.make_pairs(sample, lag)
                scores.append(Metrics.vamp2_score(z0, z1))
        return (float(np.mean(scores)), float(np.std(scores))) if scores else (np.nan, np.nan)

    @staticmethod
    def stability(history, last_k=5):
        if len(history) < 2:
            return np.nan
        arr = np.stack(history[-min(last_k, len(history)):])
        return float(np.mean(np.std(arr, axis=1)))


class Data:
    def __init__(self, cfg):
        self.cfg = cfg
        self.df = None
        self.feature_cols = []

    def load(self):
        with open(self.cfg["input_colvar"]) as f:
            fields = next(line.split()[2:] for line in f if line.startswith("#! FIELDS"))
        df = pd.read_csv(self.cfg["input_colvar"], sep=r"\s+", comment="#", names=fields)
        df = df.loc[:, ~df.columns.str.contains(r"\.")]
        if "walker" not in df.columns:
            df["walker"] = 0
        if "time" in df.columns:
            df = df[df["time"] > self.cfg["equilibration_time_ps"]].reset_index(drop=True)
        rx = re.compile(self.cfg["feature_regex"])
        self.feature_cols = [c for c in df.columns if rx.search(c)]
        if not self.feature_cols:
            raise ValueError(f"No features matched regex: {self.cfg['feature_regex']}")
        x = df[self.feature_cols].to_numpy(dtype=np.float32)
        if not np.all(np.isfinite(x)):
            raise ValueError("NaN/inf found in selected features")
        self.df = df
        return self

    def tensor(self):
        return torch.tensor(self.df[self.feature_cols].to_numpy(dtype=np.float32), dtype=torch.float32)

    def time_tensor(self):
        values = self.df["time"].to_numpy(dtype=np.float64) if "time" in self.df.columns else np.arange(len(self.df), dtype=np.float64)
        return torch.tensor(values, dtype=torch.float64)


@dataclass
class TrainOutput:
    model: object
    x_train: object
    x_test: object
    s_train: object
    s_test: object
    metrics: dict


class Trainer:
    def __init__(self, cfg, data, lag, split_ratio):
        self.cfg = cfg
        self.data = data
        self.lag_ps = int(lag)
        self.time_step = self._time_step()
        self.lag_frames = max(1, int(round(self.lag_ps / self.time_step)))
        self.split_ratio = split_ratio

    def _time_step(self):
        if self.data.df is None or "time" not in self.data.df.columns or len(self.data.df) < 2:
            return 1.0
        dt = np.diff(self.data.df["time"].to_numpy(dtype=float))
        dt = dt[np.isfinite(dt) & (dt > 0)]
        return float(np.median(dt)) if len(dt) else 1.0

    def split(self):
        x = self.data.tensor()
        t = self.data.time_tensor()
        cut = len(x) if self.split_ratio >= 1.0 else int(self.split_ratio * len(x))
        return x[:cut], x[cut:] if cut < len(x) else None, t[:cut], t[cut:] if cut < len(t) else None

    def dataset(self, x, t):
        try:
            return create_timelagged_dataset(x, t=t, lag_time=self.lag_ps)
        except TypeError:
            return create_timelagged_dataset(x, lag_time=self.lag_frames)

    def run(self):
        x_train, x_test, t_train, _ = self.split()
        if x_train.shape[1] == 0:
            raise ValueError("No input features")
        ds = self.dataset(x_train, t_train)
        dm = mlcolvar.data.DictModule(ds, lengths=[0.8, 0.2], batch_size=0, random_split=False, shuffle=False)
        layers = [x_train.shape[1], *self.cfg["hidden_layers"]]
        model = mlcolvar.cvs.DeepTICA(layers=layers, n_cvs=self.cfg["n_cvs"])
        if hasattr(model, "set_regularization"):
            model.set_regularization(c0_reg=1e-4)
        probe_n = min(self.cfg["probe_size"], len(x_train))
        probe_ids = np.random.default_rng(self.cfg["seed"]).choice(len(x_train), size=probe_n, replace=False)
        mlog = MetricsLogger(x_train[probe_ids], self.lag_frames)
        callbacks = [
            EarlyStopping(monitor="valid_loss", patience=self.cfg["patience"], min_delta=1e-4, mode="min"),
            EpochLogger(enabled=self.cfg.get("verbose", False)).callback(),
            mlog.callback(),
        ]
        trainer = L.Trainer(
            max_epochs=self.cfg["max_epochs"],
            callbacks=callbacks,
            logger=False,
            enable_checkpointing=False,
            enable_progress_bar=self.cfg.get("verbose", False),
            gradient_clip_val=1.0,
        )
        trainer.fit(model, dm)
        with torch.no_grad():
            s_train = model(x_train).detach().cpu().numpy()
            s_test = model(x_test).detach().cpu().numpy() if x_test is not None and len(x_test) > self.lag_frames else None
        metrics = self.evaluate(s_train, s_test, mlog)
        return TrainOutput(model, x_train, x_test, s_train, s_test, metrics)

    def evaluate(self, s_train, s_test, mlog):
        out = {}
        cov = np.cov(s_train.T)
        eig = np.linalg.eigvalsh(cov)
        out["final_rank"] = int(np.linalg.matrix_rank(cov))
        out["min_eig"] = float(np.min(eig))
        out["stability"] = Metrics.stability(mlog.history)
        out["vamp2_logger"] = float(np.nanmean(mlog.vamp2_history)) if mlog.vamp2_history else np.nan
        out["timescale_logger"] = float(np.nanmean(mlog.its_history)) if mlog.its_history else np.nan
        z0, z1 = Metrics.make_pairs(s_train, self.lag_frames)
        out["vamp2_train"] = Metrics.vamp2_score(z0, z1)
        if s_test is not None:
            z0t, z1t = Metrics.make_pairs(s_test, self.lag_frames)
            out["vamp2_test"] = Metrics.vamp2_score(z0t, z1t)
        else:
            out["vamp2_test"] = np.nan
        out["vamp_gap"] = float(abs(out["vamp2_train"] - out["vamp2_test"])) if np.isfinite(out["vamp2_test"]) else np.nan
        out["vamp2_boot_mean"], out["vamp2_boot_std"] = Metrics.bootstrap_vamp2(s_train, self.lag_frames)
        out.update(Metrics.diagnostics(s_train, self.lag_frames))
        out["lag_ps"] = self.lag_ps
        out["lag_frames"] = self.lag_frames
        out["time_step_ps"] = self.time_step
        return out


class Scorer:
    better = {
        "vamp2_test": True,
        "vamp_gap": False,
        "ck_spectral_error": False,
        "timescale_plateau": False,
        "spectral_gap": True,
        "vamp2_boot_rel_std": False,
        "stability": False,
    }

    @staticmethod
    def rank_norm(vals, higher=True):
        vals = np.asarray(vals, dtype=float)
        out = np.full(len(vals), 0.5, dtype=float)
        mask = np.isfinite(vals)
        if mask.sum() < 2:
            return out
        order = np.argsort(vals[mask], kind="mergesort")
        ranks = np.empty(mask.sum(), dtype=float)
        ranks[order] = np.arange(mask.sum(), dtype=float)
        ranks = ranks / max(mask.sum() - 1, 1)
        if not higher:
            ranks = 1.0 - ranks
        out[mask] = ranks
        return out

    @staticmethod
    def valid(m, n_cvs):
        ts = np.asarray(m.get("timescales", []), dtype=float)
        ts_ok = len(ts) >= 1 and np.all(np.isfinite(ts)) and np.all(ts > 0)
        if len(ts) > 1 and np.max(ts) / max(np.min(ts), 1e-12) > 1e4:
            ts_ok = False
        return (
            m.get("final_rank", 0) >= min(2, n_cvs)
            and m.get("spectral_gap_1-2", 0.0) >= 1e-4
            and m.get("ck_spectral_error", np.inf) <= 0.25
            and m.get("vamp2_test", 0.0) >= 1e-8
            and ts_ok
        )

    @staticmethod
    def score(m, n_cvs):
        if not Scorer.valid(m, n_cvs):
            return 0.0
        kinetic = m["vamp2_test_norm"]
        markov = 1.0 - m["ck_spectral_error_norm"]
        stability = 1.0 - m["vamp2_boot_rel_std_norm"]
        metastability = 0.5 * m["spectral_gap_norm"] + 0.5 * (1.0 - m["vamp_gap_norm"])
        consistency = 1.0 - m["timescale_plateau_norm"] if np.isfinite(m["timescale_plateau"]) else 0.5
        return float(0.40 * kinetic + 0.20 * markov + 0.15 * metastability + 0.15 * stability + 0.10 * consistency)

    @staticmethod
    def apply(results, n_cvs):
        lags = sorted(results)
        for lag, m in results.items():
            m["vamp2_boot_rel_std"] = m["vamp2_boot_std"] / m["vamp2_boot_mean"] if m.get("vamp2_boot_mean", 0) > 0 else np.nan
            m["spectral_gap"] = m.get("spectral_gap_1-2", np.nan)
        for i, lag in enumerate(lags):
            if i == 0:
                results[lag]["timescale_plateau"] = np.nan
                continue
            a = np.asarray(results[lags[i - 1]].get("timescales", []), dtype=float)
            b = np.asarray(results[lag].get("timescales", []), dtype=float)
            k = min(len(a), len(b), 2)
            results[lag]["timescale_plateau"] = float(np.mean(np.abs(np.log(b[:k]) - np.log(a[:k])))) if k else np.nan
        for key, higher in Scorer.better.items():
            vals = [results[l].get(key, np.nan) for l in lags]
            norm = Scorer.rank_norm(vals, higher)
            for lag, val in zip(lags, norm):
                results[lag][key + "_norm"] = float(val)
        for lag in lags:
            results[lag]["score"] = Scorer.score(results[lag], n_cvs)
            results[lag]["valid"] = Scorer.valid(results[lag], n_cvs)
        return results

    @staticmethod
    def choose(results, forced=0):
        if forced:
            if forced not in results:
                raise ValueError(f"Forced lag {forced} was not evaluated")
            return forced
        ranked = sorted(results, key=lambda lag: (results[lag]["score"], results[lag].get("vamp2_test", -np.inf)), reverse=True)
        return ranked[0]


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
        if isinstance(obj, (np.integer,)):
            return int(obj)
        if isinstance(obj, (np.floating,)):
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
        params = {
            "n_cvs": int(s.shape[1]),
            "min": np.round(mins, 3),
            "max": np.round(maxs, 3),
            "std": np.round(stds, 3),
            "grid_min": np.round(mins - margin * span, 3),
            "grid_max": np.round(maxs + margin * span, 3),
            "sigma": np.round(0.5 * stds, 3),
            "feature_cols": self.data.feature_cols,
        }
        with open(self.out / "cv_params.pkl", "wb") as f:
            pickle.dump(params, f)
        with open(self.out / "cv_params.json", "w") as f:
            json.dump(Exporter.clean(params), f, indent=2)
        return params

    def torchscript(self, model, x_ref):
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
        torch.save(checkpoint, self.out / f"CVs_model_lag{selected_lag}.pt")
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    args = ap.parse_args()
    with open(args.config) as f:
        cfg = json.load(f)
    Path(cfg["output_dir"]).mkdir(parents=True, exist_ok=True)
    load_deps()
    np.random.seed(cfg["seed"])
    torch.manual_seed(cfg["seed"])

    data = Data(cfg).load()
    print(f"features={len(data.feature_cols)} {data.feature_cols}")
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
