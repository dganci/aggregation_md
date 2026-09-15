#!/usr/bin/env bash
#
# Mac-side: build cg_md and run every test, all inside the container.
#

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CG_MD_DIR="$(cd "$HERE/.." && pwd)"

IMAGE="${IMAGE:-dganci/aggregation-md:alma9-metatomic}"
PLATFORM="${PLATFORM:-linux/amd64}"
JOBS="${JOBS:-8}"

CLEAN=0; TESTS=1
while [[ $# -gt 0 ]]; do
  case "$1" in
    --clean)    CLEAN=1; shift ;;
    --no-tests) TESTS=0; shift ;;
    *) echo "usage: $0 [--clean] [--no-tests]" >&2; exit 2 ;;
  esac
done

command -v docker >/dev/null || { echo "docker not found; start Docker Desktop" >&2; exit 2; }
docker info >/dev/null 2>&1 || { echo "Docker is not running; start Docker Desktop" >&2; exit 2; }

[[ $CLEAN -eq 1 ]] && { echo "removing $CG_MD_DIR/build"; rm -rf "$CG_MD_DIR/build"; }

inner="set -e
cd /data
cmake -S . -B build
cmake --build build -j$JOBS"

if [[ $TESTS -eq 1 ]]; then
  inner="$inner
echo
echo '=== unit tests (C++) ==='
./build/cg_md_tests
echo
echo '=== tool tests (shell) ==='
./tests/test_tools.sh"
fi

docker run --rm --platform "$PLATFORM" \
  --mount "type=bind,source=$CG_MD_DIR,target=/data" \
  "$IMAGE" bash -lc "$inner"

echo
echo "build/cg_md: $(file "$CG_MD_DIR/build/cg_md" | cut -d: -f2-)"
