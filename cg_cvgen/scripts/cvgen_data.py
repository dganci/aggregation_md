"""Loading the COLVAR(s) a CV is trained on."""
import glob
import re
from pathlib import Path

import numpy as np

from cvgen_clock import rebase_clocks
from cvgen_deps import deps
from cvgen_features import check_against_cg_md_schema, sort_permutable_blocks


class Data:
    def __init__(self, cfg):
        self.cfg = cfg
        self.df = None
        self.feature_cols = []

    def _paths(self):
        """Every COLVAR file to train on, in simulated-time order."""
        raw = self.cfg["input_colvar"]
        items = raw if isinstance(raw, (list, tuple)) else [raw]
        paths = []
        for item in items:
            item = str(item)
            matches = sorted(Path(p) for p in glob.glob(item)) if any(c in item for c in "*?[") \
                else [Path(item)]
            paths.extend(matches)
        missing = [p for p in paths if not p.is_file()]
        if missing:
            raise FileNotFoundError(f"COLVAR file(s) not found: {', '.join(map(str, missing))}")
        if not paths:
            raise FileNotFoundError(f"No COLVAR file matched: {raw}")
        return paths

    @staticmethod
    def _read_one(path):
        with open(path) as f:
            fields = next((line.split()[2:] for line in f if line.startswith("#! FIELDS")), None)
        if fields is None:
            raise ValueError(f"No '#! FIELDS' header in {path}")
        df = deps.pd.read_csv(path, sep=r"\s+", comment="#", names=fields)

        n_before = len(df)
        while len(df) and df.iloc[-1].isna().any():
            df = df.iloc[:-1]
        if len(df) < n_before:
            print(f"note: dropped {n_before - len(df)} incomplete trailing row(s) from "
                  f"{Path(path).name} (file still being written)")
        return fields, df.reset_index(drop=True)

    @staticmethod
    def _make_clock_continuous(frames):
        """Shifts each COLVAR onto a global clock before they are concatenated."""
        timed = [i for i, df in enumerate(frames) if "time" in df.columns and not df.empty]
        if not timed:
            return frames

        rebased = rebase_clocks([frames[i]["time"].to_numpy(dtype=float) for i in timed])

        out = list(frames)
        for i, t in zip(timed, rebased):
            d = frames[i].copy()
            d["time"] = t
            out[i] = d
        return out

    def load(self):
        paths = self._paths()
        walker_of, frames, reference = {}, [], None
        for path in paths:
            walker_of.setdefault(path.parent, len(walker_of))
        if len(walker_of) > 1:
            print(f"{len(walker_of)} independent trajectories (one per run directory); "
                  f"time-lagged pairs will not cross between them")
        for path in paths:
            fields, df = self._read_one(path)
            if reference is None:
                reference = fields
            elif fields != reference:
                raise ValueError(f"COLVAR header of {path} differs from {paths[0]}; "
                                 "all training files must describe the same descriptors")
            if "walker" not in df.columns:
                df = df.copy()
                df["walker"] = walker_of[path.parent]
            frames.append(df)

        frames = self._make_clock_continuous(frames)
        df = deps.pd.concat(frames, ignore_index=True) if len(frames) > 1 else frames[0]

        if "time" in df.columns:
            df = df.sort_values("time", kind="mergesort")
            df = df[~df["time"].duplicated(keep="first")]
            eq = float(self.cfg["equilibration_time_ps"])
            if eq > 0 and "walker" in df.columns:
                starts = df.groupby("walker")["time"].transform("min")
                df = df[df["time"] > starts + eq]
            elif eq > 0:
                df = df[df["time"] > df["time"].min() + eq]
            df = df.reset_index(drop=True)

        df = df.loc[:, ~df.columns.str.contains(r"\.")]

        rx = re.compile(self.cfg["feature_regex"])
        self.feature_cols = [c for c in df.columns if rx.search(c)]
        if not self.feature_cols:
            raise ValueError(f"No features matched regex: {self.cfg['feature_regex']}")
        check_against_cg_md_schema(paths, self.feature_cols)

        if self.cfg.get("permutation_invariant", True):
            df, self.feature_cols = sort_permutable_blocks(df, self.feature_cols)

        x = df[self.feature_cols].to_numpy(dtype=np.float32)
        if not np.all(np.isfinite(x)):
            raise ValueError("NaN/inf found in selected features")
        if len(df) < 2:
            raise ValueError("Fewer than two usable frames after equilibration filtering; "
                             "lower --equilibration-time-ps or sample longer")
        self.df = df
        print(f"loaded {len(df)} frames from {len(paths)} COLVAR file(s)")
        return self

    def tensor(self):
        return deps.torch.tensor(self.df[self.feature_cols].to_numpy(dtype=np.float32), dtype=deps.torch.float32)

    def walker_array(self):
        """Trajectory id per frame, aligned with tensor()/time_tensor()."""
        if "walker" not in self.df.columns:
            return None
        return self.df["walker"].to_numpy()

    def time_tensor(self):
        values = self.df["time"].to_numpy(dtype=np.float64) if "time" in self.df.columns else np.arange(len(self.df), dtype=np.float64)
        return deps.torch.tensor(values, dtype=deps.torch.float64)
