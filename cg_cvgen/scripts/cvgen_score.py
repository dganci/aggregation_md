"""Ranking the candidate lags and choosing one."""
import numpy as np


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
        gap = m.get("spectral_gap", np.nan)
        return bool(
            m.get("final_rank", 0) >= min(2, n_cvs)
            and np.isfinite(gap) and gap >= 1e-4
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
            m.setdefault("spectral_gap", np.nan)
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
        if not any(results[lag]["valid"] for lag in results):
            print("WARNING: no lag passed the validity checks (rank, spectral gap, "
                  "Chapman-Kolmogorov error, timescales). Falling back to highest VAMP-2, which "
                  "systematically favours the shortest lag. Inspect lag_scores.csv before using "
                  "this CV - the trajectory is probably too short, or the CVs are degenerate.")
        ranked = sorted(results, key=lambda lag: (results[lag]["score"], results[lag].get("vamp2_test", -np.inf)), reverse=True)
        return ranked[0]
