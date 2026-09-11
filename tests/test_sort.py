#!/usr/bin/env python3
"""Unit tests for `mytools sort`. Run without bedtools: python3 -m unittest tests/test_sort.py"""

import importlib.util
import pathlib
import unittest
from importlib.machinery import SourceFileLoader

ROOT = pathlib.Path(__file__).resolve().parent.parent
_loader = SourceFileLoader("mytools_main", str(ROOT / "mytools"))
_spec = importlib.util.spec_from_loader(_loader.name, _loader)
mytools = importlib.util.module_from_spec(_spec)
_loader.exec_module(mytools)


def sorted_lines(raw_lines):
    records = mytools.parse_bed_lines(raw_lines, "test")
    return [line for _, _, line in mytools.sort_records(records)]


class ChromOrderTests(unittest.TestCase):
    def test_chrom17_sorts_before_chrom7(self):
        # Lexicographic chrom order, not numeric: "17" < "7" as strings.
        lines = ["chr7\t10\t20\n", "chr17\t10\t20\n", "chr2\t10\t20\n"]
        self.assertEqual(
            sorted_lines(lines),
            ["chr17\t10\t20", "chr2\t10\t20", "chr7\t10\t20"],
        )


class StartTiebreakTests(unittest.TestCase):
    def test_numeric_tiebreak_on_start_within_chrom(self):
        lines = ["chr1\t200\t300\n", "chr1\t50\t60\n", "chr1\t100\t150\n"]
        self.assertEqual(
            sorted_lines(lines),
            ["chr1\t50\t60", "chr1\t100\t150", "chr1\t200\t300"],
        )


class BookendedTests(unittest.TestCase):
    def test_bookended_intervals_both_kept_in_start_order(self):
        # a.end == b.start: not an overlap, but sort just orders by start --
        # both records survive, in start order.
        lines = ["chr1\t100\t200\n", "chr1\t0\t100\n"]
        self.assertEqual(sorted_lines(lines), ["chr1\t0\t100", "chr1\t100\t200"])


class ZeroLengthTests(unittest.TestCase):
    def test_zero_length_interval_is_legal(self):
        lines = ["chr1\t500\t600\n", "chr1\t500\t500\n"]
        # Same start (500): stable sort keeps input order for the tie.
        self.assertEqual(sorted_lines(lines), ["chr1\t500\t600", "chr1\t500\t500"])


class NestedTests(unittest.TestCase):
    def test_nested_interval_sorts_by_start_not_end(self):
        lines = ["chr1\t100\t500\n", "chr1\t200\t300\n"]
        self.assertEqual(sorted_lines(lines), ["chr1\t100\t500", "chr1\t200\t300"])


class PositionZeroTests(unittest.TestCase):
    def test_interval_at_position_zero_sorts_first(self):
        lines = ["chr1\t50\t60\n", "chr1\t0\t10\n"]
        self.assertEqual(sorted_lines(lines), ["chr1\t0\t10", "chr1\t50\t60"])


class MalformedInputTests(unittest.TestCase):
    def test_too_few_columns_raises(self):
        with self.assertRaises(mytools.BedError):
            mytools.parse_bed_lines(["chr1\t100\n"], "test")

    def test_non_integer_start_raises(self):
        with self.assertRaises(mytools.BedError):
            mytools.parse_bed_lines(["chr1\tfoo\t100\n"], "test")

    def test_start_greater_than_end_raises(self):
        with self.assertRaises(mytools.BedError):
            mytools.parse_bed_lines(["chr1\t100\t50\n"], "test")

    def test_blank_and_comment_lines_are_skipped(self):
        lines = ["#comment\n", "\n", "track name=foo\n", "browser position chr1\n", "chr1\t0\t10\n"]
        self.assertEqual(sorted_lines(lines), ["chr1\t0\t10"])


if __name__ == "__main__":
    unittest.main()
