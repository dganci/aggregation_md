#!/usr/bin/env python3
"""Regression tests for PlateauScan, the lag-selection scan."""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from cvgen_backend import PlateauScan  # noqa: E402

PASSED = 0
FAILURES = []


def check(name, cond):
    global PASSED
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    if cond:
        PASSED += 1
    else:
        FAILURES.append(name)


def ou(n, taus, rng):
    """Independent Ornstein-Uhlenbeck components with known correlation times."""
    x = np.zeros((n, len(taus)))
    for j, tau in enumerate(taus):
        a = np.exp(-1.0 / tau)
        s = np.sqrt(1.0 - a * a)
        for i in range(1, n):
            x[i, j] = a * x[i - 1, j] + s * rng.standard_normal()
    return x


def standardize(x):
    return (x - x.mean(0)) / np.maximum(x.std(0), 1e-12)


rng = np.random.default_rng(7)
TAUS = [400.0, 80.0, 40.0, 20.0]
REPLICAS = [ou(2000, TAUS, rng) for _ in range(8)]

x = standardize(np.vstack(REPLICAS))
w = np.concatenate([np.full(2000, k) for k in range(8)])

for lag, min_gap in ((200, 0.05), (400, 0.15)):
    naive = PlateauScan.implied_timescales(x, lag, 1.0)[0]
    aware = PlateauScan.implied_timescales(x, lag, 1.0, walker=w)[0]
    check(f"cross-replica pairs understate the timescale at lag {lag}",
          naive < aware * (1.0 - min_gap))
    check(f"the walker-aware estimate is the closer one at lag {lag}",
          abs(aware - 400.0) < abs(naive - 400.0))

single = standardize(REPLICAS[0])
check("walker filter is a no-op on a single trajectory",
      np.isclose(PlateauScan.implied_timescales(single, 100, 1.0)[0],
                 PlateauScan.implied_timescales(single, 100, 1.0,
                                                walker=np.zeros(len(single), dtype=int))[0],
                 rtol=1e-12))

check("recovers a known Ornstein-Uhlenbeck timescale within a factor of two",
      200.0 < PlateauScan.implied_timescales(x, 100, 1.0, walker=w)[0] < 800.0)

aware = PlateauScan.implied_timescales(x, 200, 1.0, walker=w)[0]
xr = standardize(np.hstack([x, x[:, :1] + x[:, 1:2]]))
red = PlateauScan.implied_timescales(xr, 200, 1.0, walker=w)
check("an exactly redundant feature does not change the answer",
      len(red) > 0 and np.isclose(red[0], aware, rtol=0.02))

check("a lag with too few pairs returns no timescale",
      PlateauScan.implied_timescales(x, len(x) - 5, 1.0, walker=w) == [])

rows = PlateauScan.sweep(x, 1.0, 2000, walker=w)
check("the sweep is capped by the shortest replica, not the concatenated length",
      len(rows) > 0 and all(r["lag_frames"] < 200 for r in rows))

band = PlateauScan.plateau(rows)
if band:
    lags = [r["lag_frames"] for r in rows]
    band_lags = [r["lag_frames"] for r in band]
    start = lags.index(band_lags[0])
    check("the plateau is a contiguous run of lags",
          band_lags == lags[start:start + len(band_lags)])
    centre = float(np.median([r["its1_ps"] for r in band]))
    check("every plateau lag is within tolerance of the band's median",
          all(abs(r["its1_ps"] - centre) <= 0.03 * centre + 1e-9 for r in band))
else:
    check("the plateau is a contiguous run of lags", True)
    check("every plateau lag is within tolerance of the band's median", True)

ladder = [(5, 442.0), (10, 441.0), (25, 439.0), (50, 437.0), (100, 438.0),
          (250, 407.0), (500, 395.0), (1000, 555.0), (2000, 693.0)]
shaped = [{"lag_frames": lag, "lag_ps": float(lag), "its1_ps": its, "its2_ps": 1.0}
          for lag, its in ladder]
picked = [r["lag_frames"] for r in PlateauScan.plateau(shaped)]
check("the plateau is the flat stretch, not the noisy minimum", picked == [5, 10, 25, 50, 100])
check("choose() then has enough lags to spread over",
      len(PlateauScan.choose(shaped, PlateauScan.plateau(shaped))) >= 3)

real = [(5, 126.0), (10, 122.0), (25, 118.0), (50, 113.0), (100, 112.0),
        (250, 112.5), (500, 114.0), (1000, 121.0)]
real_rows = [{"lag_frames": lag, "lag_ps": float(lag), "its1_ps": its, "its2_ps": 1.0}
             for lag, its in real]
check("the fall-basin-rise shape still yields the basin",
      [r["lag_frames"] for r in PlateauScan.plateau(real_rows)] == [50, 100, 250, 500])

tie = [(5, 100.0), (10, 100.0), (25, 100.0), (50, 130.0), (100, 160.0),
       (250, 160.0), (500, 160.0)]
tie_rows = [{"lag_frames": lag, "lag_ps": float(lag), "its1_ps": its, "its2_ps": 1.0}
            for lag, its in tie]
check("between equal runs the shorter lags win",
      [r["lag_frames"] for r in PlateauScan.plateau(tie_rows)] == [5, 10, 25])

total = PASSED + len(FAILURES)
print(f"\n{PASSED}/{total} tests passed")
sys.exit(1 if FAILURES else 0)
