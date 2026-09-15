"""VAMP-2 scoring and implied timescales for one pair of lagged blocks."""
import numpy as np


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

        gap_01 = float(max(1.0 - sv[1], 0.0)) if len(sv) > 1 else np.nan
        gap_12 = float(max(sv[1] - sv[2], 0.0)) if len(sv) > 2 else np.nan
        return {
            "singular_vals": sv.astype(float).tolist(),
            "ck_spectral_error": ck,
            "spectral_gap_0-1": gap_01,
            "spectral_gap_1-2": gap_12,
            "spectral_gap": gap_12 if np.isfinite(gap_12) else gap_01,
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
