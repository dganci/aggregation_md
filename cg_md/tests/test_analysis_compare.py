#!/usr/bin/env python3
"""analyze_offline.py compare: replicas that agree, and ones that do not."""
import json
import tempfile
from pathlib import Path

import numpy as np

from analysis_common import check, close, finish, run_cli  # noqa: F401


with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)

    def summary_for(name, number_fraction, weight):
        d = root / name / "report"
        d.mkdir(parents=True)
        (d / "aggregation_summary.json").write_text(json.dumps({
            "n_prot": 3, "sampled_time_us": 1.0,
            "size_distribution": {"sizes": [1, 2, 3], "number_fraction": number_fraction,
                                  "mass_fraction": number_fraction, "mean_size_weight": weight},
            "censoring": {"fraction_at_n_prot": number_fraction[-1]}}))
        return d / "aggregation_summary.json"

    same = [summary_for("rep1", [0.6, 0.3, 0.1], 1.5),
            summary_for("rep2", [0.58, 0.32, 0.10], 1.55),
            summary_for("rep3", [0.62, 0.28, 0.10], 1.45)]
    result = run_cli("compare", "--summaries", *map(str, same), "--out-dir", str(root / "cmp"))
    check("compare runs on agreeing replicas", result.returncode == 0)
    agree = json.loads((root / "cmp" / "replica_comparison.json").read_text())
    check("agreeing replicas are called agreeing", agree.get("replicas_agree") is True)
    check("the pooled weight-average is the mean of the replicas",
          close(agree["pooled"]["mean_size_weight"], 1.5))
    check("the between-replica error is the SEM",
          close(agree["pooled"]["mean_size_weight_sem"], np.std([1.5, 1.55, 1.45], ddof=1) / np.sqrt(3)))
    check("three replicas give three pairs", len(agree["pairwise_jsd_bits"]) == 3)

    odd = summary_for("rep4", [0.1, 0.1, 0.8], 2.7)
    result = run_cli("compare", "--summaries", *map(str, same + [odd]), "--out-dir", str(root / "cmp2"))
    disagree = json.loads((root / "cmp2" / "replica_comparison.json").read_text())
    check("a replica that assembled where the rest did not is caught",
          disagree.get("replicas_agree") is False and "rep4" in disagree["worst_pair"])

    result = run_cli("compare", "--summaries", str(same[0]), "--out-dir", str(root / "cmp3"))
    check("one summary is not a comparison", result.returncode == 1)


finish()
