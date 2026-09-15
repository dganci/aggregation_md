"""The `analyze_offline.py compare` subcommand: do the replicas agree?"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

import figures_aggregation as figs
from offline_io import dump, expand
from report_io import HAVE_MPL


def jensen_shannon_bits(p, q):
    """JSD of two probability vectors on the same support, in bits."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / p.sum() if p.sum() > 0 else p
    q = q / q.sum() if q.sum() > 0 else q
    m = 0.5 * (p + q)

    def kl(a):
        keep = a > 0
        return float(np.sum(a[keep] * np.log2(a[keep] / m[keep])))

    return 0.5 * kl(p) + 0.5 * kl(q)


def load_summaries(paths):
    """(label, summary) for every readable aggregation_summary.json."""
    out = []
    for path in paths:
        try:
            summary = json.loads(Path(path).read_text())
        except (OSError, json.JSONDecodeError) as error:
            print(f"analyze_offline: cannot read {path}: {error}", file=sys.stderr)
            continue
        if "size_distribution" not in summary or "error" in summary["size_distribution"]:
            print(f"analyze_offline: {path} has no size distribution - skipped", file=sys.stderr)
            continue
        out.append((label_for(path), summary))
    return out


def label_for(path):
    """The run's name: the file is always aggregation_summary.json and the pipeline puts it in
    <run>/report/, so neither of those tells runs apart.
    """
    parent = Path(path).resolve().parent
    if parent.name == "report" and parent.parent.name:
        parent = parent.parent
    return parent.name


def run(args) -> int:
    replicas = load_summaries(expand(args.summaries))
    if len(replicas) < 2:
        print("analyze_offline: compare needs at least two summaries", file=sys.stderr)
        return 1

    n_prot = {s["n_prot"] for _, s in replicas}
    if len(n_prot) != 1:
        print(f"analyze_offline: the summaries have different --n-prot ({sorted(n_prot)}); "
              "their size distributions are not comparable", file=sys.stderr)
        return 1
    n_prot = n_prot.pop()

    labels = [label for label, _ in replicas]
    number = np.array([s["size_distribution"]["number_fraction"] for _, s in replicas])
    mass = np.array([s["size_distribution"]["mass_fraction"] for _, s in replicas])
    weight_avg = np.array([s["size_distribution"]["mean_size_weight"] for _, s in replicas])
    censoring = np.array([s.get("censoring", {}).get("fraction_at_n_prot", np.nan)
                          for _, s in replicas])
    sampled = np.array([s.get("sampled_time_us", np.nan) for _, s in replicas])

    pairs = []
    for i in range(len(replicas)):
        for j in range(i + 1, len(replicas)):
            pairs.append({"a": labels[i], "b": labels[j],
                          "jsd_bits": jensen_shannon_bits(number[i], number[j])})
    worst = max(pairs, key=lambda p: p["jsd_bits"])

    n = len(replicas)
    sem = lambda a: np.nanstd(a, axis=0, ddof=1) / np.sqrt(n)   # noqa: E731

    summary = {
        "replicas": labels,
        "n_prot": n_prot,
        "sampled_time_us": sampled.tolist(),
        "pairwise_jsd_bits": pairs,
        "max_pairwise_jsd_bits": worst["jsd_bits"],
        "worst_pair": [worst["a"], worst["b"]],
        "jsd_threshold_bits": args.max_jsd,
        "replicas_agree": bool(worst["jsd_bits"] <= args.max_jsd),
        "pooled": {
            "sizes": list(range(1, n_prot + 1)),
            "number_fraction_mean": number.mean(axis=0).tolist(),
            "number_fraction_sem": sem(number).tolist(),
            "mass_fraction_mean": mass.mean(axis=0).tolist(),
            "mass_fraction_sem": sem(mass).tolist(),
            "mean_size_weight": float(weight_avg.mean()),
            "mean_size_weight_sem": float(sem(weight_avg)),
            "fraction_at_n_prot": float(np.nanmean(censoring)),
            "fraction_at_n_prot_sem": float(sem(censoring)),
        },
        "per_replica": {label: {"mean_size_weight": float(w), "fraction_at_n_prot": float(c),
                                "sampled_time_us": float(t)}
                        for label, w, c, t in zip(labels, weight_avg, censoring, sampled)},
    }

    out_dir = Path(args.out_dir)
    produced: list[str] = []
    if HAVE_MPL:
        figs.plot_replica_agreement(labels, number, summary["pooled"], pairs, out_dir, produced)
    summary["figures"] = produced
    dump(summary, out_dir, "replica_comparison.json")

    verdict = "agree" if summary["replicas_agree"] else "DISAGREE"
    print(f"analyze_offline: {n} replicas {verdict}: worst pair {worst['a']} / {worst['b']} "
          f"at {worst['jsd_bits']:.3f} bits (threshold {args.max_jsd}); "
          f"<s>_w = {summary['pooled']['mean_size_weight']:.2f} "
          f"+/- {summary['pooled']['mean_size_weight_sem']:.2f}")
    return 0
