"""mytools intersect -- report overlaps between two BED files.

Flag semantics and output columns are meant to match `bedtools intersect`
exactly (issue #14 / SPEC.md sec4). See "Overlap semantics" below for the
one place this deliberately goes beyond the documented strict-inequality
predicate in SPEC.md sec3.

Two-file CLI: SPEC.md sec9 only spells out a single positional [file|-] for
one-input subcommands (sort, merge). intersect takes two files, so this
mirrors bedtools' own -a/-b flags instead -- flagged in issue #14 as a spec
gap rather than something pinned down by SPEC.md itself.
"""

import sys

from mytools_lib.bedio import iter_bed_records, strand_of
from mytools_lib.errors import BedError, UsageError

USAGE = "usage: mytools intersect -a <file|-> -b <file|-> [-wa] [-wb] [-v] [-u] [-f F] [-r] [-s] [-S]"

DEFAULT_F = 1e-9


class Args:
    def __init__(self):
        self.a = None
        self.b = None
        self.wa = False
        self.wb = False
        self.v = False
        self.u = False
        self.f = None
        self.r = False
        self.s = False
        self.S = False


def parse_args(argv):
    args = Args()
    i = 0
    while i < len(argv):
        tok = argv[i]
        if tok == "-a":
            i += 1
            if i >= len(argv):
                raise UsageError(f"{USAGE}\nmissing value for -a")
            args.a = argv[i]
        elif tok == "-b":
            i += 1
            if i >= len(argv):
                raise UsageError(f"{USAGE}\nmissing value for -b")
            args.b = argv[i]
        elif tok == "-wa":
            args.wa = True
        elif tok == "-wb":
            args.wb = True
        elif tok == "-v":
            args.v = True
        elif tok == "-u":
            args.u = True
        elif tok == "-r":
            args.r = True
        elif tok == "-s":
            args.s = True
        elif tok == "-S":
            args.S = True
        elif tok == "-f":
            i += 1
            if i >= len(argv):
                raise UsageError(f"{USAGE}\nmissing value for -f")
            try:
                args.f = float(argv[i])
            except ValueError:
                raise UsageError(f"{USAGE}\n-f requires a float, got {argv[i]!r}")
        else:
            raise UsageError(f"{USAGE}\nunrecognized option: {tok}")
        i += 1

    if args.a is None or args.b is None:
        raise UsageError(f"{USAGE}\n-a and -b are both required")
    return args


# ---------------------------------------------------------------------------
# Overlap semantics
#
# SPEC.md sec3's strict predicate (a.start < b.end AND b.start < a.end) is
# what real bedtools uses for two normal-width intervals. It can't be the
# whole story for zero-length features though: a zero-length interval never
# satisfies start < end against anything, yet real bedtools does report
# overlaps for them. Reverse-engineered against real bedtools (v2.31.1) by
# probing a.bed/b.bed plus hand-built cases (a zero-length point at a
# normal interval's start/end/interior, two zero-length features at the
# same position and at varying offsets, and -f/-r fraction thresholds
# swept across a zero-length pair): bedtools appears to treat *every*
# zero-length feature [s, s) as if it were the 2bp window [s-1, s+1) for
# the purposes of (a) deciding whether an overlap exists and (b) how much
# overlap there is (both plain and -f/-r fraction math), while continuing
# to *report* a zero-length A's own true, unexpanded coordinates in
# default-mode output. That single expansion rule reproduces every case
# tried, including the otherwise-inexplicable one where two zero-length
# features exactly 1bp apart are reported as overlapping while ones 2bp+
# apart are not.
# ---------------------------------------------------------------------------


def _expand_if_zero(start, end):
    if start == end:
        return start - 1, end + 1
    return start, end


def overlaps(a_start, a_end, b_start, b_end):
    ea_start, ea_end = _expand_if_zero(a_start, a_end)
    eb_start, eb_end = _expand_if_zero(b_start, b_end)
    return ea_start < eb_end and eb_start < ea_end


def clipped_coords(a_start, a_end, b_start, b_end):
    """The reported intersection region for the *default* output mode.

    Assumes overlaps(a_start, a_end, b_start, b_end) is already True.

    A zero-length A always reports its own coordinates unchanged -- it's
    trivially contained in any interval it overlaps, and bedtools doesn't
    expand it for display purposes. Otherwise, clip A's true coordinates
    against B's coordinates -- expanded per the zero-length rule above if
    B is zero-length, which is what turns an otherwise-degenerate
    (zero-width) clip into the single adjacent base bedtools reports.
    """
    if a_start == a_end:
        return a_start, a_end
    eb_start, eb_end = _expand_if_zero(b_start, b_end)
    return max(a_start, eb_start), min(a_end, eb_end)


def fraction_ok(a_start, a_end, b_start, b_end, f, reciprocal):
    """-f / -r: minimum overlap as a fraction of A (and of B, if reciprocal).

    Both overlap length and each feature's own length use the zero-length
    expansion rule above (a zero-length feature's effective length is 2),
    confirmed by sweeping -f across a zero-length pair against real
    bedtools and finding the pass/fail threshold sits at exactly 1/2.
    """
    ea_start, ea_end = _expand_if_zero(a_start, a_end)
    eb_start, eb_end = _expand_if_zero(b_start, b_end)
    ov = max(0, min(ea_end, eb_end) - max(ea_start, eb_start))
    a_len = ea_end - ea_start
    if (ov / a_len) < f:
        return False
    if not reciprocal:
        return True
    b_len = eb_end - eb_start
    return (ov / b_len) >= f


def _load_b(path):
    # B is loaded fully into memory; only A is streamed line-by-line.
    # SPEC.md sec6 says intersect streams over "assumed-sorted" input, but
    # real bedtools' *default* (non -sorted) mode -- which is what -a/-b
    # without -sorted actually runs, and what byte-for-byte matching the
    # oracle on the deliberately-unsorted a.bed/b.bed requires -- doesn't
    # assume sorted input and needs all of B available to find every A
    # record's matches in B's own file order. -sorted (chromsweep, fully
    # streaming both sides) is out of scope for issue #14.
    by_chrom = {}
    for rec in iter_bed_records(path):
        by_chrom.setdefault(rec.chrom, []).append(rec)
    return by_chrom


def _find_matches(a_rec, b_by_chrom, args):
    f = args.f if args.f is not None else DEFAULT_F
    a_strand = strand_of(a_rec.fields) if (args.s or args.S) else None
    matches = []
    for b_rec in b_by_chrom.get(a_rec.chrom, ()):
        if args.s or args.S:
            b_strand = strand_of(b_rec.fields)
            same = a_strand == b_strand
            if args.s and not same:
                continue
            if args.S and same:
                continue
        if not overlaps(a_rec.start, a_rec.end, b_rec.start, b_rec.end):
            continue
        if args.f is not None or args.r:
            if not fraction_ok(a_rec.start, a_rec.end, b_rec.start, b_rec.end, f, args.r):
                continue
        matches.append(b_rec)
    return matches


def run(argv, out=None):
    out = out or sys.stdout
    args = parse_args(argv)

    b_by_chrom = _load_b(args.b)

    for a_rec in iter_bed_records(args.a):
        matches = _find_matches(a_rec, b_by_chrom, args)

        if args.v:
            if not matches:
                out.write(a_rec.raw + "\n")
            continue

        if args.u:
            if matches:
                out.write(a_rec.raw + "\n")
            continue

        for b_rec in matches:
            if args.wa and args.wb:
                out.write(a_rec.raw + "\t" + b_rec.raw + "\n")
            elif args.wa:
                out.write(a_rec.raw + "\n")
            else:
                lo, hi = clipped_coords(a_rec.start, a_rec.end, b_rec.start, b_rec.end)
                clipped_fields = list(a_rec.fields)
                clipped_fields[1] = str(lo)
                clipped_fields[2] = str(hi)
                clipped = "\t".join(clipped_fields)
                if args.wb:
                    out.write(clipped + "\t" + b_rec.raw + "\n")
                else:
                    out.write(clipped + "\n")

    return 0


def main(argv):
    try:
        return run(argv)
    except UsageError as e:
        print(str(e), file=sys.stderr)
        return 2
    except BedError as e:
        print(str(e), file=sys.stderr)
        return 1
