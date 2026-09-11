"""Shared BED I/O helpers for mytools subcommands.

BED is 0-based, half-open: `chrom  start  end` covers bases start..end-1.
"""

import collections
import gzip
import io
import sys

_SKIP_PREFIXES = ("track", "browser", "#")


class BedDataError(Exception):
    """Bad input data: malformed line, missing file, unsorted input, etc.

    Maps to exit code 1.
    """


class UsageError(Exception):
    """Bad command-line usage: unknown flag, bad flag value/combination.

    Maps to exit code 2.
    """


BedRecord = collections.namedtuple("BedRecord", "chrom start end fields lineno")


def display_name(path):
    return "stdin" if path in (None, "-") else path


def _is_gzip_magic(head):
    return head[:2] == b"\x1f\x8b"


def _open_binary(path):
    if path in (None, "-"):
        return sys.stdin.buffer
    try:
        return open(path, "rb")
    except OSError:
        raise BedDataError(
            f"mytools: cannot open {path}: No such file or directory"
        )


def open_text(path):
    """Open a BED source as a decoded, line-iterable text stream.

    `path` of None or "-" reads stdin. Gzip input is detected by magic bytes
    (not just the `.gz` extension), matching bedtools' own auto-detection.
    """
    raw = _open_binary(path)
    if path in (None, "-"):
        try:
            head = raw.peek(2)[:2]
        except (AttributeError, OSError):
            head = b""
        if _is_gzip_magic(head):
            return io.TextIOWrapper(gzip.GzipFile(fileobj=raw))
        return io.TextIOWrapper(raw)
    else:
        head = raw.read(2)
        raw.seek(0)
        if _is_gzip_magic(head):
            return io.TextIOWrapper(gzip.GzipFile(fileobj=raw))
        return io.TextIOWrapper(raw)


def iter_records(path):
    """Yield BedRecord for every data line of a BED source.

    Blank lines and `track`/`browser`/`#` lines are skipped, matching
    bedtools. A malformed line (too few columns, non-numeric start/end,
    negative start, or start > end) raises BedDataError naming the file and
    line number.
    """
    name = display_name(path)
    fh = open_text(path)
    is_stdin = path in (None, "-")
    lineno = 0
    try:
        for raw_line in fh:
            lineno += 1
            line = raw_line.rstrip("\n").rstrip("\r")
            if line == "" or line.startswith(_SKIP_PREFIXES):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise BedDataError(
                    f"mytools: malformed BED line in {name} at line {lineno}: "
                    f"expected at least 3 tab-separated fields, got {len(fields)}"
                )
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                raise BedDataError(
                    f"mytools: malformed BED line in {name} at line {lineno}: "
                    "start and end must be integers"
                )
            if start < 0:
                raise BedDataError(
                    f"mytools: malformed BED line in {name} at line {lineno}: "
                    f"start ({start}) must be >= 0"
                )
            if start > end:
                raise BedDataError(
                    f"mytools: malformed BED line in {name} at line {lineno}: "
                    f"start ({start}) > end ({end})"
                )
            yield BedRecord(chrom, start, end, fields, lineno)
    finally:
        if not is_stdin:
            fh.close()


def require_sorted(records, name):
    """Wrap a BedRecord iterator, raising BedDataError on out-of-order input.

    Matches bedtools' own sortedness requirement: chromosomes may appear in
    any order, but records within a chromosome must be non-decreasing by
    start, and a chromosome's records must not be interleaved with another
    chromosome's (i.e. a chromosome may not reappear once we've moved on).
    """
    seen_chroms = set()
    last_chrom = None
    last_start = None
    for rec in records:
        if rec.chrom != last_chrom:
            if rec.chrom in seen_chroms:
                raise BedDataError(
                    f"mytools: sorted input specified, but {name} has an "
                    f"out-of-order record at line {rec.lineno}: "
                    f"{rec.chrom}\t{rec.start}\t{rec.end} "
                    f"(chromosome {rec.chrom!r} appeared earlier and is not contiguous)"
                )
            seen_chroms.add(rec.chrom)
            last_chrom = rec.chrom
            last_start = rec.start
        else:
            if rec.start < last_start:
                raise BedDataError(
                    f"mytools: sorted input specified, but {name} has an "
                    f"out-of-order record at line {rec.lineno}: "
                    f"{rec.chrom}\t{rec.start}\t{rec.end}"
                )
            last_start = rec.start
        yield rec


def overlaps(a_start, a_end, b_start, b_end):
    """The core BED overlap predicate (SPEC.md #3): strict `<` on both sides.

    Bookended intervals (one's end equals the other's start) do NOT overlap
    by this predicate -- they only merge under `merge -d 0` as a separate,
    explicit bookend allowance, not because this function says they overlap.
    """
    return a_start < b_end and b_start < a_end


def field(rec, col):
    """1-based column lookup into a record's raw fields, as used by -c/-o."""
    if col < 1 or col > len(rec.fields):
        raise BedDataError(
            f"mytools: requested column {col}, but a record at line {rec.lineno} "
            f"only has fields 1-{len(rec.fields)}"
        )
    return rec.fields[col - 1]
