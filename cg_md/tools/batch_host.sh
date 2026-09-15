#!/usr/bin/env bash
#
# Mac-side launcher: starts the container and runs the whole batch inside it.
# This is the ONE command to leave running.
#

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CG_MD_DIR="$(cd "$HERE/.." && pwd)"

IMAGE="${IMAGE:-dganci/aggregation-md:alma9-metatomic}"
CONTAINER="${CONTAINER:-cg_md_batch}"
MANIFEST="${MANIFEST:-tools/sims.tsv}"
CPUS="${CPUS:-8}"
MEMORY="${MEMORY:-24g}"
JOBS="${JOBS:-2}"
NTOMP="${NTOMP:-}"
PER_SIM_TIMEOUT="${PER_SIM_TIMEOUT:-0}"
PLATFORM="${PLATFORM:-linux/amd64}"

MODE=""; ONLY=""; FOREGROUND=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --check)      MODE="check"; FOREGROUND=1; shift ;;
    --foreground) FOREGROUND=1; shift ;;
    --only)       ONLY="${2:-}"; shift 2
                  [[ -n "$ONLY" ]] || { echo "--only needs an entry name" >&2; exit 2; } ;;
    *) echo "usage: $0 [--check|--foreground] [--only NAME]" >&2; exit 2 ;;
  esac
done

command -v docker >/dev/null || { echo "docker not found; start Docker Desktop" >&2; exit 2; }
docker info >/dev/null 2>&1 || { echo "Docker is not running; start Docker Desktop" >&2; exit 2; }
[[ -f "$CG_MD_DIR/$MANIFEST" ]] || { echo "manifest not found: $CG_MD_DIR/$MANIFEST" >&2; exit 2; }

[[ -x "$CG_MD_DIR/build/cg_md" ]] || {
  echo "build/cg_md is missing or not executable." >&2
  echo "Build it inside the container first:" >&2
  echo "  docker run --rm --platform $PLATFORM --mount type=bind,source=$CG_MD_DIR,target=/data \\" >&2
  echo "    $IMAGE bash -lc 'cd /data && cmake -S . -B build && cmake --build build -j8'" >&2
  exit 2
}
if ! file "$CG_MD_DIR/build/cg_md" | grep -q "ELF.*x86-64"; then
  echo "build/cg_md is not a Linux x86-64 binary:" >&2
  file "$CG_MD_DIR/build/cg_md" >&2
  echo "Rebuild it inside the container (see DOCKER.md)." >&2
  exit 2
fi

if docker ps --format '{{.Names}}' | grep -qx "$CONTAINER"; then
  echo "'$CONTAINER' is already running. Follow it with:  docker logs -f $CONTAINER" >&2
  exit 2
fi
if [[ "$MODE" != "check" ]] && docker ps -a --format '{{.Names}}' | grep -qx "$CONTAINER"; then
  echo "removing the previous '$CONTAINER' container (its docker logs go with it;" >&2
  echo "runs/_batch/batch.out keeps the same output across runs)." >&2
  docker rm "$CONTAINER" >/dev/null
fi

arch="$(uname -m)"
if [[ "$arch" == "arm64" && "$PLATFORM" == "linux/amd64" ]]; then
  echo "note      Apple Silicon running an amd64 image under emulation."
  echo "          Measured faster than native arm64 GROMACS on this machine; the"
  echo "          numbers are in this script's own comments."
fi

mkdir -p "$CG_MD_DIR/runs/_batch"
probe="mount_probe.$$"
printf 'host\n' > "$CG_MD_DIR/runs/_batch/$probe"
simd_file="$CG_MD_DIR/runs/_batch/$probe.simd"
if ! docker run --rm --platform "$PLATFORM" \
      --mount "type=bind,source=$CG_MD_DIR,target=/data" "$IMAGE" \
      bash -lc "test -f /data/runs/_batch/$probe \
                && printf 'container\n' > /data/runs/_batch/$probe.back \
                && gmx_mpi --version 2>/dev/null | sed -n 's/^SIMD instructions: *//p' \
                   > /data/runs/_batch/$probe.simd" 2>/dev/null; then
  rm -f "$CG_MD_DIR/runs/_batch/$probe"
  echo "The bind mount does not expose this directory to the container." >&2
  echo "  source: $CG_MD_DIR -> /data" >&2
  echo "In Docker Desktop, check Settings > Resources > File sharing, and that" >&2
  echo "Docker has permission to read this folder (System Settings > Privacy)." >&2
  exit 2
fi
if [[ ! -f "$CG_MD_DIR/runs/_batch/$probe.back" ]]; then
  rm -f "$CG_MD_DIR/runs/_batch/$probe"
  echo "The container can read the mount but its writes do not reach the Mac." >&2
  echo "Nothing the batch produces would be saved. Restart Docker Desktop and retry." >&2
  exit 2
fi
simd="$(cat "$simd_file" 2>/dev/null | tr -d '\r')"
if [[ -z "$simd" ]]; then
  echo "warning   could not read GROMACS's SIMD level from the image." >&2
elif [[ "$simd" != AVX* ]]; then
  echo >&2
  echo "  WARNING: GROMACS reports SIMD '$simd', not an AVX level." >&2
  echo "  A compute-bound run will be several times slower than it should be." >&2
  echo "  On Apple Silicon this usually means Docker Desktop is emulating with" >&2
  echo "  QEMU instead of Rosetta: enable" >&2
  echo "    Settings > General > 'Use Rosetta for x86_64/amd64 emulation'" >&2
  echo "  and re-run. Set ASSUME_YES=1 to start anyway." >&2
  echo >&2
  if [[ "${ASSUME_YES:-0}" != "1" ]]; then
    if [[ -t 0 ]]; then
      read -r -p "  Continue anyway? [y/N] " reply
      [[ "$reply" == [yY] ]] || { rm -f "$CG_MD_DIR/runs/_batch/$probe"* ; exit 1; }
    else
      rm -f "$CG_MD_DIR/runs/_batch/$probe"*
      exit 1
    fi
  fi
else
  echo "simd      GROMACS $simd"
fi
rm -f "$CG_MD_DIR/runs/_batch/$probe" "$CG_MD_DIR/runs/_batch/$probe.back" "$simd_file"

args=(
  --platform "$PLATFORM"
  --cpus "$CPUS"
  --memory "$MEMORY"
  --mount "type=bind,source=$CG_MD_DIR,target=/data"
  -e "NTOMP=${NTOMP:-0}"
  -e "JOBS=$JOBS"
  -e "PER_SIM_TIMEOUT=$PER_SIM_TIMEOUT"
  -e "PROJECT_DIR=/data"
  -e "CG_MD=/data/build/cg_md"
  -e "LOG_ROOT=/data/runs/_batch"
)
# cv_readiness.py (the ITS gate of the stop rule) looks for cg_cvgen/scripts
# two levels above itself; mounting the sibling checkout there keeps the host's
# cg_md and cg_cvgen in step. Without it the gate falls back to the clone
# inside the image ($CGMD_SRC) - or never passes.
CVGEN_DIR="$(cd "$CG_MD_DIR/.." && pwd)/cg_cvgen"
if [[ -d "$CVGEN_DIR/scripts" ]]; then
  args+=(--mount "type=bind,source=$CVGEN_DIR,target=/cg_cvgen,readonly")
else
  echo "warning   $CVGEN_DIR not found: the CV-trainability gate will use the" >&2
  echo "          cg_cvgen inside the image, if the image has one." >&2
fi
[[ "$MODE" == "check" ]] || args+=(--name "$CONTAINER")
[[ $FOREGROUND -eq 1 ]] || args+=(-d --restart=no)

echo "image     $IMAGE ($PLATFORM)"
echo "mount     $CG_MD_DIR -> /data  (verified both ways)"
[[ -d "$CVGEN_DIR/scripts" ]] && echo "mount     $CVGEN_DIR -> /cg_cvgen  (read-only, for cv_readiness.py)"
if [[ -n "$NTOMP" ]]; then
  echo "resources --cpus $CPUS --memory $MEMORY, $JOBS concurrent x --ntomp $NTOMP (NTOMP forced)"
else
  echo "resources --cpus $CPUS --memory $MEMORY, $JOBS concurrent runs, --ntomp derived as CPUS/JOBS"
fi
echo "manifest  $MANIFEST"
[[ -n "$ONLY" ]] && echo "only      $ONLY"
echo

if [[ "$MODE" == "check" ]]; then
  exec docker run --rm "${args[@]}" "$IMAGE" \
    bash -lc "cd /data && ./tools/preflight.sh '$MANIFEST'"
fi

only_arg=""
[[ -n "$ONLY" ]] && only_arg=" --only '$ONLY'"
inner="set -o pipefail; cd /data && mkdir -p /data/runs/_batch && \
./tools/run_batch.sh '$MANIFEST'$only_arg 2>&1 | tee -a /data/runs/_batch/batch.out"

if [[ $FOREGROUND -eq 1 ]]; then
  docker run "${args[@]}" "$IMAGE" bash -lc "$inner"
  exit $?
fi

docker run "${args[@]}" "$IMAGE" bash -lc "$inner" >/dev/null
echo "started in the background as '$CONTAINER'."

if [[ "${NO_CAFFEINATE:-0}" != "1" ]] && command -v caffeinate >/dev/null; then
  nohup caffeinate -dims bash -c \
    "while [ \"\$(docker inspect -f '{{.State.Running}}' $CONTAINER 2>/dev/null)\" = true ]; do sleep 60; done" \
    >/dev/null 2>&1 &
  echo "caffeinate is holding the Mac awake until the container exits (pid $!)."
else
  echo "NOTE: sleep prevention is off. Keep the Mac plugged in and set Energy Saver"
  echo "      to never sleep on power, or the container will be paused."
fi

cat <<EOF

  docker logs -f $CONTAINER
  tail -f $CG_MD_DIR/runs/_batch/batch.out
  tail -20 $CG_MD_DIR/runs/_batch/status.txt
  column -t -s\$'\t' $CG_MD_DIR/runs/_batch/summary.tsv

Stop with:  docker stop $CONTAINER   (re-launching resumes)

When it ends, check how it ended - a batch that fails still exits:
  docker inspect $CONTAINER --format '{{.State.Status}} exit={{.State.ExitCode}}'
EOF
