#!/usr/bin/env bash
# Reports the measured throughput of every run under runs/, and what a given
# budget buys at that rate.
#
#   tools/throughput.sh [node-hours]        (default: 62.5, the try26 allocation)
#
# Two sources, in order of trust:
#   - a COMPLETED chunk writes "Performance: <ns/day>" into its md_chunk_NNN.log.
#     That is a real measurement over the whole chunk and is preferred.
#   - a chunk still in flight has no such line, so the rate is derived from
#     mdrun -v's own "step N, will finish <date>" against the log's mtime. That
#     is an estimate: mdrun's projection, read at the moment it was written.
#
# Note that a chunk killed by SIGKILL leaves neither. SLURM sends SIGTERM first
# only if the job asks; see --signal in the sbatch script.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUDGET_NODE_H="${1:-62.5}"
GPUS_PER_NODE="${GPUS_PER_NODE:-4}"
CHUNK_US="${CHUNK_US:-0.1}"

printf '  %-22s %-10s %12s  %s\n' entry chunk "ns/day" source
found=0; total=0; n=0
for d in "$ROOT"/runs/*/; do
    name="$(basename "$d")"
    [[ "$name" == _* || "$name" == relax_* ]] && continue
    last_done=""; rate=""; src=""
    for log in "$d"md_chunk_*.log; do
        [[ -e "$log" ]] || continue
        if p="$(grep -a 'Performance:' "$log" 2>/dev/null | tail -1 | awk '{print $2}')" && [[ -n "$p" ]]; then
            last_done="$(basename "$log" .log)"; rate="$p"; src="completed chunk"
        fi
    done
    if [[ -z "$rate" ]]; then
        # Newest chunk log, still running: use mdrun's own projection.
        log="$(ls -t "$d"md_chunk_*.log 2>/dev/null | head -1 || true)"
        [[ -n "$log" ]] || continue
        read -r step when < <(tail -c 4000 "$log" | tr '\r' '\n' |
            sed -n 's/^step \([0-9]*\), will finish \(.*\)$/\1 \2/p' | tail -1) || true
        [[ -n "${step:-}" ]] || continue
        rate="$(python3 - "$step" "$when" "$log" "$CHUNK_US" <<'PY'
import sys, datetime, os
step, when, log, chunk_us = int(sys.argv[1]), sys.argv[2], sys.argv[3], float(sys.argv[4])
total = int(chunk_us * 1e6 / 0.01)          # dt = 10 fs
end = datetime.datetime.strptime(when.strip(), "%a %b %d %H:%M:%S %Y")
now = datetime.datetime.fromtimestamp(os.path.getmtime(log))
left = (end - now).total_seconds()
print(f"{(total-step)/left*0.01*86400/1000:.1f}" if left > 0 and total > step else "")
PY
)"
        last_done="$(basename "$log" .log) (in flight)"; src="mdrun projection"
    fi
    [[ -n "$rate" ]] || continue
    printf '  %-22s %-10s %12s  %s\n' "$name" "$last_done" "$rate" "$src"
    total="$(python3 -c "print($total + $rate)")"; n=$((n+1)); found=1
done

[[ "$found" -eq 1 ]] || { echo "  no chunk has reported a rate yet"; exit 0; }

python3 - "$total" "$n" "$BUDGET_NODE_H" "$GPUS_PER_NODE" "$CHUNK_US" <<'PY'
import sys
total, n, budget, gpus, chunk = float(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5])
mean = total / n
print(f"\n  mean over {n} run(s): {mean:.0f} ns/day per replica")
print(f"  one {chunk} us chunk:  {chunk*1000/mean*24:.1f} h")
print(f"\n  at {gpus} replicas per node, {budget} node-hours buy:")
agg = mean * gpus * budget / 24 / 1000
print(f"    {agg:.1f} us of aggregate sampling")
for label, need in (("10 entries x 1 us  (minimum scenario)", 10.0),
                    ("15 entries x 1 us", 15.0),
                    ("15 entries x 10 us (ceiling)", 150.0)):
    print(f"    {label:<38} {'covered' if agg >= need else f'{agg/need*100:.0f}% of it'}")
print(f"\n  size the next job to a multiple of the chunk time: "
      f"{chunk*1000/mean*24:.1f} h, {2*chunk*1000/mean*24:.1f} h, {3*chunk*1000/mean*24:.1f} h ...")
PY
