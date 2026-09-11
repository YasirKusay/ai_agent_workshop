"""`mytools merge` -- merge overlapping (and, at -d 0, bookended) BED intervals.

Interval semantics (see SPEC.md #3 and CLAUDE.md):

  - BED is 0-based, half-open. Two intervals overlap iff
    `a.start < b.end AND b.start < a.end` (strict `<`).
  - Bookended intervals (`a.end == b.start`) do not overlap, but DO merge
    under the default `-d 0`.
  - Zero-length intervals (`start == end`) are legal. Real bedtools pads a
    zero-length feature's coordinates by 1bp on each side -- but only for the
    purpose of deciding adjacency and computing a merged cluster's bounds
    once it actually merges with something. A zero-length feature that ends
    up alone in its own cluster is reported with its original, unpadded
    coordinates. This is empirically verified against real bedtools (not
    guessed): e.g. `chr1 500 500` merged with a neighbouring `chr1 500 600`
    is reported as `chr1 499 600`, but `chr1 500 500` on its own is reported
    unchanged.

Input must be pre-sorted the same way bedtools requires it: chromosomes may
appear in any order, but a chromosome's records must be contiguous and
non-decreasing by start.
"""

import heapq
import itertools
import sys

from . import bed
from .bed import BedDataError, UsageError

_STRAND_COL = 6  # 1-based column for BED6 strand


def _pad(rec):
    """Padded (start, end) used for adjacency tests and cluster bounds."""
    if rec.start == rec.end:
        return rec.start - 1, rec.end + 1
    return rec.start, rec.end


def _strand_of(rec):
    if len(rec.fields) < _STRAND_COL:
        raise BedDataError(
            f"mytools: -s/-S requires a strand column (BED6+), but a record "
            f"at line {rec.lineno} only has {len(rec.fields)} fields"
        )
    return rec.fields[_STRAND_COL - 1]


class _Cluster:
    __slots__ = ("chrom", "first_start", "first_end", "pstart", "pend", "n", "col_values")

    def __init__(self, rec, cols):
        p0, p1 = _pad(rec)
        self.chrom = rec.chrom
        self.first_start = rec.start
        self.first_end = rec.end
        self.pstart = p0
        self.pend = p1
        self.n = 1
        self.col_values = [[bed.field(rec, c)] for c in cols] if cols else None

    def add(self, rec, cols):
        p0, p1 = _pad(rec)
        self.pstart = min(self.pstart, p0)
        self.pend = max(self.pend, p1)
        self.n += 1
        if cols:
            for values, c in zip(self.col_values, cols):
                values.append(bed.field(rec, c))

    def bounds(self):
        if self.n == 1:
            return self.first_start, self.first_end
        return self.pstart, self.pend


def _cluster_single(records, d, cols):
    """One running cluster at a time -- default merge, or merge over a
    single strand already filtered by -S."""
    current = None
    for rec in records:
        if current is not None and _pad(rec)[0] <= current.pend + d:
            current.add(rec, cols)
        else:
            if current is not None:
                yield current
            current = _Cluster(rec, cols)
    if current is not None:
        yield current


def _cluster_forced_strand(records, d, cols):
    """-s: only merge features sharing a strand. Behaves as an independent
    merge sweep per strand value, with the results of all strands combined
    back in position order -- verified empirically against real bedtools
    (a same-strand pair merges across an intervening different-strand
    feature that sits between them in the sorted input).

    Implemented as a streaming k-way merge: one open cluster per strand
    value seen so far, plus a small heap of clusters that have closed but
    are not yet safe to emit (because another still-open cluster might still
    turn out to start earlier once padding from a zero-length feature is
    accounted for).
    """
    trackers = {}
    pending = []  # heap of (bounds_start, seq, cluster)
    seq = itertools.count()

    def push_pending(cluster):
        heapq.heappush(pending, (cluster.bounds()[0], next(seq), cluster))

    def safety_floor():
        # Conservative lower bound on any open tracker's eventual start,
        # using each tracker's padded pstart (never larger than its final
        # bounds()[0]) so we never flush a pending cluster too early.
        open_pstarts = [t.pstart for t in trackers.values() if t is not None]
        return min(open_pstarts) if open_pstarts else None

    for rec in records:
        strand = _strand_of(rec)
        current = trackers.get(strand)
        if current is not None and _pad(rec)[0] <= current.pend + d:
            current.add(rec, cols)
        else:
            if current is not None:
                push_pending(current)
            trackers[strand] = _Cluster(rec, cols)

        floor = safety_floor()
        if floor is not None:
            while pending and pending[0][0] <= floor:
                yield heapq.heappop(pending)[2]

    for t in trackers.values():
        if t is not None:
            push_pending(t)
    while pending:
        yield heapq.heappop(pending)[2]


_NUMERIC_OPS = frozenset(
    {"sum", "min", "max", "absmin", "absmax", "mean", "median", "mode", "antimode", "stdev", "sstdev"}
)
_VALID_OPS = _NUMERIC_OPS | {
    "collapse", "distinct", "distinct_sort_num", "distinct_sort_num_desc",
    "count", "count_distinct", "first", "last",
}


def _fmt_num(x):
    return f"{x:.10g}"


def _as_float(v):
    try:
        return float(v)
    except ValueError:
        raise BedDataError(
            f"mytools: non-numeric value {v!r} for a numeric merge operation"
        )


def _apply_op(op, values, delim):
    if op in _NUMERIC_OPS:
        nums = [_as_float(v) for v in values]
        if op == "sum":
            return _fmt_num(sum(nums))
        if op == "min":
            return _fmt_num(min(nums))
        if op == "max":
            return _fmt_num(max(nums))
        if op == "absmin":
            return _fmt_num(min(nums, key=abs))
        if op == "absmax":
            return _fmt_num(max(nums, key=abs))
        if op == "mean":
            return _fmt_num(sum(nums) / len(nums))
        if op == "median":
            s = sorted(nums)
            mid = len(s) // 2
            if len(s) % 2:
                return _fmt_num(s[mid])
            return _fmt_num((s[mid - 1] + s[mid]) / 2)
        if op in ("mode", "antimode"):
            counts = {}
            for v in nums:
                counts[v] = counts.get(v, 0) + 1
            best = max if op == "mode" else min
            target = best(counts.values())
            candidates = [v for v, c in counts.items() if c == target]
            return _fmt_num(min(candidates))
        if op in ("stdev", "sstdev"):
            mean = sum(nums) / len(nums)
            denom = len(nums) if op == "stdev" else len(nums) - 1
            if denom <= 0:
                # sample stdev is undefined for a single value
                return "."
            variance = sum((v - mean) ** 2 for v in nums) / denom
            return _fmt_num(variance ** 0.5)
    if op == "collapse":
        return delim.join(values)
    if op == "distinct":
        return delim.join(sorted(set(values)))
    if op in ("distinct_sort_num", "distinct_sort_num_desc"):
        seen = sorted({_as_float(v) for v in values}, reverse=(op == "distinct_sort_num_desc"))
        return delim.join(_fmt_num(v) for v in seen)
    if op == "count":
        return str(len(values))
    if op == "count_distinct":
        return str(len(set(values)))
    if op == "first":
        return values[0]
    if op == "last":
        return values[-1]
    raise AssertionError(f"unhandled op {op!r}")


def _format_columns(cluster, ops, delim):
    out = []
    for values, op in zip(cluster.col_values, ops):
        out.append(_apply_op(op, values, delim))
    return out


USAGE = "usage: mytools merge [options]"


def parse_args(argv):
    opts = {"input": None, "d": 0, "c": None, "o": None, "s": False, "S": None}
    positional = []
    i = 0
    while i < len(argv):
        a = argv[i]
        if a == "-i":
            i += 1
            if i >= len(argv):
                raise UsageError("-i requires an argument")
            opts["input"] = argv[i]
        elif a == "-d":
            i += 1
            if i >= len(argv):
                raise UsageError("-d requires an integer argument")
            try:
                opts["d"] = int(argv[i])
            except ValueError:
                raise UsageError(f"-d requires an integer argument, got {argv[i]!r}")
        elif a == "-c":
            i += 1
            if i >= len(argv):
                raise UsageError("-c requires an argument")
            opts["c"] = argv[i]
        elif a == "-o":
            i += 1
            if i >= len(argv):
                raise UsageError("-o requires an argument")
            opts["o"] = argv[i]
        elif a == "-s":
            opts["s"] = True
        elif a == "-S":
            i += 1
            if i >= len(argv):
                raise UsageError("-S requires + or -")
            if argv[i] not in ("+", "-"):
                raise UsageError(f"-S requires + or -, got {argv[i]!r}")
            opts["S"] = argv[i]
        elif a == "-":
            positional.append(a)
        elif a.startswith("-") and len(a) > 1:
            raise UsageError(f"unrecognized argument: {a}")
        else:
            positional.append(a)
        i += 1

    if opts["input"] is None and positional:
        opts["input"] = positional[0]

    if opts["s"] and opts["S"] is not None:
        raise UsageError("-s and -S are mutually exclusive")

    return opts


def _resolve_cols_ops(opts):
    c_given = opts["c"] is not None
    o_given = opts["o"] is not None
    if c_given != o_given:
        raise UsageError("-c and -o must be given together")
    if not c_given:
        return None, None, ","

    try:
        cols = [int(x) for x in opts["c"].split(",")]
    except ValueError:
        raise UsageError(f"-c expects a comma-separated list of column numbers, got {opts['c']!r}")
    if any(c < 1 for c in cols):
        raise UsageError("-c column numbers must be >= 1")

    ops = opts["o"].split(",")
    for op in ops:
        if op not in _VALID_OPS:
            raise UsageError(f"{op!r} is not a valid -o operation")

    if len(cols) == len(ops):
        pass
    elif len(cols) == 1:
        cols = cols * len(ops)
    elif len(ops) == 1:
        ops = ops * len(cols)
    else:
        raise UsageError(
            f"-c has {len(cols)} column(s) but -o has {len(ops)} operation(s); "
            "provide one operation for all columns, one column for all operations, "
            "or a matching number of each"
        )

    return cols, ops, ","


def run(argv):
    opts = parse_args(argv)
    cols, ops, delim = _resolve_cols_ops(opts)
    name = bed.display_name(opts["input"])

    records = bed.iter_records(opts["input"])
    records = bed.require_sorted(records, name)

    if opts["S"] is not None:
        records = (r for r in records if _strand_of(r) == opts["S"])

    out = sys.stdout
    for _, chrom_records in itertools.groupby(records, key=lambda r: r.chrom):
        if opts["s"]:
            clusters = _cluster_forced_strand(chrom_records, opts["d"], cols)
        else:
            clusters = _cluster_single(chrom_records, opts["d"], cols)
        for cluster in clusters:
            start, end = cluster.bounds()
            row = [cluster.chrom, str(start), str(end)]
            if cols:
                row.extend(_format_columns(cluster, ops, delim))
            out.write("\t".join(row) + "\n")
    return 0


def main(argv):
    try:
        return run(argv)
    except UsageError as e:
        print(USAGE, file=sys.stderr)
        print(f"mytools merge: {e}", file=sys.stderr)
        return 2
    except BedDataError as e:
        print(str(e), file=sys.stderr)
        return 1
