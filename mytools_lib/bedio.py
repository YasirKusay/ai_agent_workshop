"""Minimal streaming BED reader shared by mytools subcommands.

BED is 0-based, half-open (SPEC.md sec3). This module only parses chrom/
start/end plus whatever other columns are present; it never interprets or
reformats columns beyond the three required ones, so name/score/strand/etc.
round-trip byte-for-byte.
"""

import sys
from collections import namedtuple

from mytools_lib.errors import BedError

BedRecord = namedtuple("BedRecord", ["chrom", "start", "end", "fields", "raw"])

_SKIP_PREFIXES = ("#", "track", "browser")


def open_input(path):
    """Open path for reading text, or stdin if path is '-'.

    Raises BedError with the SPEC.md sec7 missing-file message.
    """
    if path == "-":
        return sys.stdin
    try:
        return open(path, "r")
    except OSError:
        raise BedError(f"mytools: cannot open {path}: No such file or directory")


def _label(path):
    return "<stdin>" if path == "-" else path


def iter_bed_records(path):
    """Yield BedRecord for each data line in path (or stdin for '-').

    Blank lines and track/browser/# lines are skipped, matching bedtools
    (verified against real bedtools on fixtures with such lines).
    Raises BedError, naming the file and line number, on a malformed line:
    fewer than 3 columns, non-integer start/end, or start > end.
    """
    fh = open_input(path)
    label = _label(path)
    try:
        for lineno, line in enumerate(fh, start=1):
            raw = line.rstrip("\n")
            if raw == "" or raw.startswith(_SKIP_PREFIXES):
                continue
            fields = raw.split("\t")
            if len(fields) < 3:
                raise BedError(
                    f"mytools: {label}:{lineno}: malformed BED line: "
                    f"expected at least 3 columns, got {len(fields)}"
                )
            chrom = fields[0]
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                raise BedError(
                    f"mytools: {label}:{lineno}: malformed BED line: "
                    "start/end must be integers"
                )
            if start > end:
                raise BedError(
                    f"mytools: {label}:{lineno}: malformed BED line: "
                    f"start ({start}) > end ({end})"
                )
            yield BedRecord(chrom, start, end, fields, raw)
    finally:
        if fh is not sys.stdin:
            fh.close()


def strand_of(fields):
    """Column 6 (1-based) if present, else None."""
    return fields[5] if len(fields) > 5 else None
