#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools.
# Usage: ./tests/run_golden.sh
set -uo pipefail

MYTOOLS=${MYTOOLS:-"$(cd "$(dirname "$0")/.." && pwd)/mytools"}
DATA=$(dirname "$0")/../data
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0

# check <name> -- <args...>
#   runs "$MYTOOLS <args>" and "bedtools <args>", diffs them
check() {
  local name=$1; shift; shift        # drop the literal --
  "$MYTOOLS" "$@" > "$tmp/got"  2>"$tmp/got.err"
  local got_rc=$?
  bedtools   "$@" > "$tmp/want" 2>/dev/null
  local want_rc=$?

  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return
  fi
  if diff -q "$tmp/want" "$tmp/got" >/dev/null; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name"
    diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
    (( fail++ ))
  fi
}

# checkstdin <name> <stdin-file> -- <args...>
#   like check(), but feeds stdin-file on stdin to both tools
checkstdin() {
  local name=$1 stdinfile=$2; shift; shift; shift   # drop the literal --
  "$MYTOOLS" "$@" < "$stdinfile" > "$tmp/got"  2>"$tmp/got.err"
  local got_rc=$?
  bedtools   "$@" < "$stdinfile" > "$tmp/want" 2>/dev/null
  local want_rc=$?

  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return
  fi
  if diff -q "$tmp/want" "$tmp/got" >/dev/null; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name"
    diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
    (( fail++ ))
  fi
}

# --- intersect (issue #14) -------------------------------------------------
# `sort`/`merge`/`subtract`/`closest` aren't implemented yet (separate
# issues) so this suite only covers `intersect` for now -- add their cases
# back in alongside their implementations.
# a.bed/b.bed are deliberately unsorted; intersect (unlike merge/closest)
# doesn't require sorted input in bedtools' default (non -sorted) mode, so
# these run straight against the fixtures as committed.

check "intersect default"        -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed"
check "intersect -wa"            -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -wa
check "intersect -wb"            -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -wb
check "intersect -wa -wb"        -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -wa -wb
check "intersect -v"             -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -v
check "intersect -u"             -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -u
check "intersect -s"             -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -s
check "intersect -S"             -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -S
check "intersect -f 0.5"         -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -f 0.5
check "intersect -f 0.5 -r"      -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -f 0.5 -r
check "intersect -u -s"          -- intersect -a "$DATA/a.bed" -b "$DATA/b.bed" -u -s
check "intersect a vs genes.bed" -- intersect -a "$DATA/a.bed" -b "$DATA/genes.bed"
checkstdin "intersect -a stdin" "$DATA/a.bed" -- intersect -a - -b "$DATA/b.bed"

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
