"""Training one DeepTICA model at one lag, and scoring what came out."""
from dataclasses import dataclass

import numpy as np

from cvgen_deps import deps
from cvgen_metrics import Metrics


class EpochLogger:
    def __init__(self, interval=25, enabled=False):
        self.interval = interval
        self.enabled = enabled

    def callback(self):
        parent = self

        class _Logger(deps.Callback):
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

        class _Metrics(deps.Callback):
            def on_train_epoch_end(self, trainer, pl_module):
                with deps.torch.no_grad():
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
        self.lag_frames = max(1, int(lag))
        self.time_step = self._time_step()
        self.lag_ps = self.lag_frames * self.time_step
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
        w = self.data.walker_array()
        cut = len(x) if self.split_ratio >= 1.0 else int(self.split_ratio * len(x))
        self._walker_train = None if w is None else w[:cut]
        return x[:cut], x[cut:] if cut < len(x) else None, t[:cut], t[cut:] if cut < len(t) else None

    def dataset(self, x, t, walker=None):
        """Time-lagged pairs, never spanning two independent trajectories."""
        try:
            return deps.create_timelagged_dataset(x, t=t, lag_time=self.lag_ps, walker=walker)
        except TypeError:
            try:
                return deps.create_timelagged_dataset(x, lag_time=self.lag_frames, walker=walker)
            except TypeError:
                if walker is not None and len(set(np.asarray(walker).tolist())) > 1:
                    raise RuntimeError(
                        "This mlcolvar build does not accept the `walker` argument, so time-lagged "
                        "pairs cannot be prevented from spanning independent replicas. Train on one "
                        "run directory at a time, or upgrade deps.mlcolvar.")
                return deps.create_timelagged_dataset(x, lag_time=self.lag_frames)

    def run(self):
        x_train, x_test, t_train, _ = self.split()
        if x_train.shape[1] == 0:
            raise ValueError("No input features")
        ds = self.dataset(x_train, t_train, getattr(self, '_walker_train', None))
        dm = deps.mlcolvar.data.DictModule(ds, lengths=[0.8, 0.2], batch_size=0, random_split=False, shuffle=False)
        layers = [x_train.shape[1], *self.cfg["hidden_layers"]]
        model = deps.mlcolvar.cvs.DeepTICA(layers=layers, n_cvs=self.cfg["n_cvs"])
        if hasattr(model, "set_regularization"):
            model.set_regularization(c0_reg=1e-4)
        probe_n = min(self.cfg["probe_size"], len(x_train))
        probe_ids = np.random.default_rng(self.cfg["seed"]).choice(len(x_train), size=probe_n, replace=False)
        mlog = MetricsLogger(x_train[probe_ids], self.lag_frames)
        callbacks = [
            deps.EarlyStopping(monitor="valid_loss", patience=self.cfg["patience"], min_delta=1e-4, mode="min"),
            EpochLogger(enabled=self.cfg.get("verbose", False)).callback(),
            mlog.callback(),
        ]
        trainer = deps.L.Trainer(
            max_epochs=self.cfg["max_epochs"],
            callbacks=callbacks,
            logger=False,
            enable_checkpointing=False,
            enable_progress_bar=self.cfg.get("verbose", False),
            gradient_clip_val=1.0,
        )
        trainer.fit(model, dm)
        with deps.torch.no_grad():
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
