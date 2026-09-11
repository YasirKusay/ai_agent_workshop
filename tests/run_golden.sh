#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools.
# Usage: ./tests/run_golden.sh
set -uo pipefail

# Resolve to *this* checkout's mytools, not whatever "mytools" happens to be
# on PATH -- multiple checkouts of this repo can coexist on one box.
MYTOOLS=${MYTOOLS:-"$(cd "$(dirname "$0")/.." && pwd)/mytools"}
DATA=$(dirname "$0")/../data
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0

report() {
  local name=$1 got_rc=$2 want_rc=$3
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

# check_diff_args <name> -- <mytools args...> -- <bedtools args...>
#   Use when mytools' own CLI for a subcommand differs from bedtools' (e.g.
#   sort takes a plain file arg / "-" for stdin, where bedtools requires -i).
check_diff_args() {
  local name=$1; shift; shift
  local mt_args=()
  while [[ "$1" != "--" ]]; do mt_args+=("$1"); shift; done
  shift
  local bt_args=("$@")

  "$MYTOOLS" "${mt_args[@]}" > "$tmp/got" 2>"$tmp/got.err"
  local got_rc=$?
  bedtools "${bt_args[@]}" > "$tmp/want" 2>/dev/null
  local want_rc=$?
  report "$name" "$got_rc" "$want_rc"
}

# check_diff_args_stdin <name> <input_file> -- <mytools args...> -- <bedtools args...>
#   Same as check_diff_args, but feeds <input_file> to both tools' stdin.
check_diff_args_stdin() {
  local name=$1 input=$2; shift; shift; shift
  local mt_args=()
  while [[ "$1" != "--" ]]; do mt_args+=("$1"); shift; done
  shift
  local bt_args=("$@")

  "$MYTOOLS" "${mt_args[@]}" < "$input" > "$tmp/got" 2>"$tmp/got.err"
  local got_rc=$?
  bedtools "${bt_args[@]}" < "$input" > "$tmp/want" 2>/dev/null
  local want_rc=$?
  report "$name" "$got_rc" "$want_rc"
}

# check <name> -- <args...>
#   Same args to both -- for subcommands whose flags mirror bedtools exactly
#   (merge, intersect, subtract, closest).
check() {
  local name=$1; shift; shift
  check_diff_args "$name" -- "$@" -- "$@"
}

check_diff_args      "sort a.bed (file arg)" -- sort "$DATA/a.bed" -- sort -i "$DATA/a.bed"
check_diff_args      "sort b.bed (file arg)" -- sort "$DATA/b.bed" -- sort -i "$DATA/b.bed"
check_diff_args_stdin "sort a.bed (stdin)" "$DATA/a.bed" -- sort - -- sort -i -

# Add the rest here as later subcommands land. merge and closest need sorted
# input -- sort into $tmp first.

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
