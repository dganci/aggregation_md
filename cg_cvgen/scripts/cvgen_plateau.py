"""The implied-timescale scan that chooses the lags to train at."""
import csv
from pathlib import Path

import numpy as np


class PlateauScan:
    """Picks the lag times from the implied-timescale plateau of the data itself."""

    @staticmethod
    def implied_timescales(x, lag, dt_ps, walker=None, cutoff=1e-8):
        """Two slowest implied timescales, in ps, from linear TICA at `lag`."""
        idx = np.arange(len(x) - lag)
        if walker is not None:
            idx = idx[walker[idx] == walker[idx + lag]]
        if len(idx) < 10 * x.shape[1]:
            return []

        a, b = x[idx], x[idx + lag]
        n = len(a)
        c0 = (a.T @ a + b.T @ b) / (2.0 * n)
        ct = (a.T @ b + b.T @ a) / (2.0 * n)

        w, v = np.linalg.eigh(c0)
        keep = w > cutoff * w.max()
        if not np.any(keep):
            return []
        isq = v[:, keep] / np.sqrt(w[keep])
        m = isq.T @ ct @ isq
        ev = np.sort(np.linalg.eigvalsh(0.5 * (m + m.T)))[::-1]
        ev = ev[(ev > 0.0) & (ev < 1.0)]
        return [float(-lag * dt_ps / np.log(e)) for e in ev[:2]]

    @staticmethod
    def sweep(x, dt_ps, n_frames, walker=None):
        """ITS over a geometric ladder of lags, capped so pairs stay plentiful."""
        ladder = [5, 10, 25, 50, 100, 250, 500, 1000, 2000, 4000, 8000]
        rows = []
        for lag in ladder:
            if lag >= max(n_frames // 10, 2):
                break
            try:
                its = PlateauScan.implied_timescales(x, lag, dt_ps, walker=walker)
            except Exception:                       # noqa: BLE001 - a singular block is not fatal
                continue
            if len(its) >= 1 and np.isfinite(its[0]):
                rows.append({"lag_frames": lag,
                             "lag_ps": lag * dt_ps,
                             "its1_ps": its[0],
                             "its2_ps": its[1] if len(its) > 1 else float("nan")})
        return rows

    @staticmethod
    def plateau(rows, tol=0.03):
        """The longest run of consecutive lags whose ITS-1 agree to within `tol`."""
        if not rows:
            return []
        values = [r["its1_ps"] for r in rows]

        def agrees(indices):
            centre = float(np.median([values[k] for k in indices]))
            return all(abs(values[k] - centre) <= tol * centre for k in indices)

        best = []
        for start in range(len(rows)):
            run = [start]
            for j in range(start + 1, len(rows)):
                if not agrees(run + [j]):
                    break
                run.append(j)
            if len(run) > len(best):
                best = run
        return [rows[k] for k in best]

    @staticmethod
    def choose(rows, plateau, n_lags=4):
        """Up to `n_lags` lags spread geometrically across the plateau."""
        if len(plateau) < 3:
            return []
        lags = [r["lag_frames"] for r in plateau]
        if len(lags) <= n_lags:
            return lags
        idx = np.unique(np.round(np.linspace(0, len(lags) - 1, n_lags)).astype(int))
        return [lags[i] for i in idx]

    @staticmethod
    def run(cfg, data, out_dir):
        """Returns the lags to score, and writes the sweep to its_scan.csv."""
        x = data.df[data.feature_cols].to_numpy(dtype=np.float64)
        x = (x - x.mean(axis=0)) / np.maximum(x.std(axis=0), 1e-12)
        dt_ps = float(np.median(np.diff(data.df["time"].to_numpy()))) if "time" in data.df else 1.0

        walker = data.walker_array()
        if walker is None:
            shortest = len(x)
        else:
            counts = np.bincount(walker)
            counts = counts[counts > 0]
            shortest = int(counts.min()) if counts.size else len(x)
        rows = PlateauScan.sweep(x, dt_ps, shortest, walker=walker)
        band = PlateauScan.plateau(rows)
        chosen = PlateauScan.choose(rows, band)

        out = Path(out_dir) / "its_scan.csv"
        with open(out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["lag_frames", "lag_ps", "its1_ps", "its2_ps",
                                              "in_plateau", "selected"])
            w.writeheader()
            band_lags, chosen_lags = {r["lag_frames"] for r in band}, set(chosen)
            for r in rows:
                w.writerow({**r,
                            "in_plateau": int(r["lag_frames"] in band_lags),
                            "selected": int(r["lag_frames"] in chosen_lags)})

        print(f"implied-timescale scan ({dt_ps:.1f} ps per frame) -> {out}")
        for r in rows:
            mark = "*" if r["lag_frames"] in chosen_lags else (
                "+" if r["lag_frames"] in band_lags else " ")
            print(f"  {mark} lag {r['lag_frames']:>5} frames ({r['lag_ps'] / 1000:6.2f} ns)"
                  f"  ITS-1 {r['its1_ps'] / 1000:8.2f} ns")
        if chosen:
            print(f"plateau: {band[0]['lag_frames']}-{band[-1]['lag_frames']} frames; "
                  f"scoring lags {chosen}")
        else:
            print("WARNING: no implied-timescale plateau found (the leading timescale never "
                  "stops depending on the lag). Falling back to the configured --lags. This "
                  "usually means the trajectory is too short for its own slowest mode - check "
                  "its_scan.csv before trusting any CV trained on it.")
        return chosen
