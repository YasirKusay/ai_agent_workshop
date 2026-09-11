#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools.
# Usage: ./tests/run_golden.sh
set -uo pipefail

MYTOOLS=${MYTOOLS:-$(cd "$(dirname "$0")/.." && pwd)/mytools}     # override to test a different build
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

# `sort` (#12) isn't implemented yet -- its golden cases land in that issue's commit.
# merge needs sorted input -- sort a.bed/b.bed into $tmp once, up front (via unix
# sort, standing in for `mytools sort` until it exists).
sort -t $'\t' -k1,1 -k2,2n "$DATA/a.bed" > "$tmp/a.sorted.bed"
sort -t $'\t' -k1,1 -k2,2n "$DATA/b.bed" > "$tmp/b.sorted.bed"

check "merge a.bed default"            -- merge -i "$tmp/a.sorted.bed"
check "merge b.bed default"            -- merge -i "$tmp/b.sorted.bed"
check "merge a.bed -d 10"              -- merge -i "$tmp/a.sorted.bed" -d 10
check "merge a.bed -d 0 (explicit)"    -- merge -i "$tmp/a.sorted.bed" -d 0
check "merge a.bed -d -5 (min overlap)" -- merge -i "$tmp/a.sorted.bed" -d -5
check "merge a.bed -c 4 -o collapse"   -- merge -i "$tmp/a.sorted.bed" -c 4 -o collapse
check "merge a.bed -c 4,5 -o collapse,sum" -- merge -i "$tmp/a.sorted.bed" -c 4,5 -o collapse,sum
check "merge a.bed -c 5 -o mean,median,min,max,absmin,absmax,stdev,sstdev,mode,antimode,count,count_distinct,first,last" \
  -- merge -i "$tmp/a.sorted.bed" -c 5,5,5,5,5,5,5,5,5,5,5,5,5,5 \
     -o mean,median,min,max,absmin,absmax,stdev,sstdev,mode,antimode,count,count_distinct,first,last
check "merge a.bed -c 4,6 -o distinct" -- merge -i "$tmp/a.sorted.bed" -c 4,6 -o distinct
check "merge a.bed -c 5 -o distinct_sort_num,distinct_sort_num_desc" \
  -- merge -i "$tmp/a.sorted.bed" -c 5,5 -o distinct_sort_num,distinct_sort_num_desc
check "merge a.bed -s"                 -- merge -i "$tmp/a.sorted.bed" -s
check "merge a.bed -S +"               -- merge -i "$tmp/a.sorted.bed" -S +
check "merge a.bed -S -"               -- merge -i "$tmp/a.sorted.bed" -S -
check "merge a.bed -s -c 4,5 -o collapse,sum" -- merge -i "$tmp/a.sorted.bed" -s -c 4,5 -o collapse,sum
check "merge b.bed -d 20"              -- merge -i "$tmp/b.sorted.bed" -d 20
check "merge b.bed -s -c 4,5,6 -o collapse,sum,distinct" \
  -- merge -i "$tmp/b.sorted.bed" -s -c 4,5,6 -o collapse,sum,distinct

# stdin needs its own invocation shape (no -i arg passed to the binaries).
check_stdin() {
  local name=$1
  "$MYTOOLS" merge < "$tmp/a.sorted.bed" > "$tmp/got" 2>"$tmp/got.err"
  local got_rc=$?
  bedtools merge < "$tmp/a.sorted.bed" > "$tmp/want" 2>/dev/null
  local want_rc=$?
  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"; (( fail++ )); return
  fi
  if diff -q "$tmp/want" "$tmp/got" >/dev/null; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name"; diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20; (( fail++ ))
  fi
}
check_stdin "merge from stdin"

check "merge unsorted input errors"    -- merge -i "$DATA/a.bed"

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
