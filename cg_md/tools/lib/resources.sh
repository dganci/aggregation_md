# shellcheck shell=bash
resources_detect() {
  cores="$(nproc 2>/dev/null || echo 4)"
  [[ "$JOBS" -ge 1 ]] || { echo "JOBS must be >= 1" >&2; exit 2; }
  if [[ "$NTOMP" -le 0 ]]; then
    usable=$(( cores > 4 ? cores - 2 : cores ))
    NTOMP=$(( usable / JOBS ))
    [[ "$NTOMP" -ge 1 ]] || NTOMP=1
  fi
  mem_gb="$(awk '/MemTotal/ {printf "%d", $2/1048576}' /proc/meminfo 2>/dev/null || echo 0)"
  if [[ "$JOBS" -gt 1 ]]; then
    echo "resources: $cores cores, ${mem_gb} GB RAM -> $JOBS concurrent runs x --ntomp $NTOMP"
    echo "           (aggregate throughput; each entry takes longer than it would alone)"
  else
    echo "resources: $cores cores, ${mem_gb} GB RAM -> --ntomp $NTOMP, one run at a time"
  fi
  [[ "$mem_gb" -gt 0 && "$mem_gb" -lt 8 ]] && echo "WARNING only ${mem_gb} GB visible to the container; raise Docker's memory limit"

  simd="$("${CG_MD_GMX:-gmx_mpi}" --version 2>/dev/null | sed -n 's/^SIMD instructions: *//p' | tr -d '\r')"
  if [[ -z "$simd" ]]; then
    echo "WARNING could not determine GROMACS's SIMD level"
  elif [[ "$simd" != AVX* ]]; then
    echo "WARNING GROMACS SIMD is '$simd', not an AVX level: a compute-bound run will be"
    echo "        several times slower than it should be. See tools/batch_host.sh for"
    echo "        the measured emulated-vs-native numbers."
  else
    echo "simd: GROMACS $simd"
  fi
}
