#!/usr/bin/env bash
#
# Tests for the batch tooling in tools/ - the shell half of the pipeline, which
# had no tests at all when a 15-entry batch silently produced nothing.
#

set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"

if [[ -z "${BASH_VERSINFO:-}" || "${BASH_VERSINFO[0]}" -lt 4 ]]; then
  echo "these tests need bash 4+ (macOS ships 3.2); run them inside the container" >&2
  exit 2
fi
for part in "$HERE"/tools/[0-9]*.sh; do
  # shellcheck source=/dev/null
  source "$part"
done

echo
echo "$((PASS))/$((PASS + FAIL)) tool tests passed"
[[ $FAIL -eq 0 ]]
