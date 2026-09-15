# shellcheck shell=bash
source "$ROOT/tools/lib/manifest.sh"
source "$ROOT/tools/lib/logdir.sh"

PASS=0; FAIL=0
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

ok()   { PASS=$((PASS + 1)); echo "[PASS] $1"; }
bad()  { FAIL=$((FAIL + 1)); echo "[FAIL] $1: $2"; }
check()      { if [[ "$2" == "$3" ]]; then ok "$1"; else bad "$1" "expected '$3', got '$2'"; fi; }
check_true() { if [[ "$2" -eq 0 ]]; then ok "$1"; else bad "$1" "expected success, got status $2"; fi; }
check_fail() { if [[ "$2" -ne 0 ]]; then ok "$1"; else bad "$1" "expected failure, got success"; fi; }
