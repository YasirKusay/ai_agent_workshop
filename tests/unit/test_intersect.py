"""Unit tests for mytools_lib.intersect. Run without bedtools installed.

    python3 -m unittest discover -s tests/unit
"""

import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from mytools_lib.intersect import clipped_coords, fraction_ok, overlaps  # noqa: E402

MYTOOLS = os.path.join(ROOT, "mytools")


def run_mytools(args, stdin_text=None):
    return subprocess.run(
        [sys.executable, MYTOOLS, *args],
        input=stdin_text,
        capture_output=True,
        text=True,
    )


class TestOverlapPredicate(unittest.TestCase):
    """SPEC.md sec3: a.start < b.end AND b.start < a.end, strict, for two
    normal-width (non-zero-length) intervals."""

    def test_clear_overlap(self):
        self.assertTrue(overlaps(100, 200, 150, 250))

    def test_no_overlap_disjoint(self):
        self.assertFalse(overlaps(100, 200, 300, 400))

    def test_nested_is_overlap(self):
        # a06 fully inside a05, per data/a.bed
        self.assertTrue(overlaps(300, 400, 320, 350))

    def test_position_zero(self):
        self.assertTrue(overlaps(0, 100, 0, 50))
        self.assertFalse(overlaps(0, 100, 100, 200))


class TestBookended(unittest.TestCase):
    """Bookended (a.end == b.start) must NOT count as overlap: SPEC.md sec3."""

    def test_bookend_a_before_b(self):
        self.assertFalse(overlaps(0, 100, 100, 200))

    def test_bookend_b_before_a(self):
        self.assertFalse(overlaps(100, 200, 0, 100))


class TestZeroLength(unittest.TestCase):
    """Zero-length (start == end) features: legal, and bedtools treats them
    as if expanded to [s-1, s+1) for overlap/fraction purposes -- see the
    "Overlap semantics" note in mytools_lib/intersect.py for how this was
    derived against real bedtools v2.31.1. These pin the cases that were
    checked against the oracle."""

    def test_zero_length_matches_at_normal_start(self):
        # point at a normal interval's own start counts as overlapping
        self.assertTrue(overlaps(500, 600, 500, 500))

    def test_zero_length_matches_at_normal_end(self):
        # bookended from the point's side still counts -- unlike two
        # normal-width bookended intervals
        self.assertTrue(overlaps(0, 100, 100, 100))

    def test_zero_length_matches_interior(self):
        self.assertTrue(overlaps(31, 43, 41, 41))

    def test_zero_length_no_match_far_away(self):
        self.assertFalse(overlaps(100, 200, 500, 500))

    def test_zero_length_vs_zero_length_same_point(self):
        self.assertTrue(overlaps(500, 500, 500, 500))

    def test_zero_length_vs_zero_length_1bp_apart(self):
        # a genuine bedtools quirk (both treated as 2bp windows one base
        # apart still touch) -- confirmed against the oracle, not guessed
        self.assertTrue(overlaps(500, 500, 501, 501))

    def test_zero_length_vs_zero_length_2bp_apart_no_match(self):
        self.assertFalse(overlaps(500, 500, 502, 502))

    def test_clipped_coords_zero_length_a_keeps_own_coords(self):
        self.assertEqual(clipped_coords(500, 500, 400, 600), (500, 500))

    def test_clipped_coords_zero_length_b_interior(self):
        self.assertEqual(clipped_coords(31, 43, 41, 41), (40, 42))

    def test_clipped_coords_zero_length_b_at_a_start(self):
        self.assertEqual(clipped_coords(31, 43, 31, 31), (31, 32))

    def test_clipped_coords_zero_length_b_at_a_end(self):
        self.assertEqual(clipped_coords(31, 43, 43, 43), (42, 43))

    def test_fraction_zero_length_effective_length_is_two(self):
        # a07-style pair: zero-length A's fraction denominator behaves as
        # if its length were 2, so a 1bp overlap passes -f 0.5 but not
        # -f 0.6 (crossover confirmed against real bedtools)
        self.assertTrue(fraction_ok(7, 7, 6, 7, 0.5, reciprocal=False))
        self.assertFalse(fraction_ok(7, 7, 6, 7, 0.6, reciprocal=False))


class TestCLI(unittest.TestCase):
    def test_default_streams_from_file(self):
        with tempfile.TemporaryDirectory() as d:
            a = os.path.join(d, "a.bed")
            b = os.path.join(d, "b.bed")
            with open(a, "w") as f:
                f.write("chr1\t0\t100\ta1\t0\t+\n")
            with open(b, "w") as f:
                f.write("chr1\t50\t150\tb1\t0\t+\n")
            result = run_mytools(["intersect", "-a", a, "-b", b])
            self.assertEqual(result.returncode, 0)
            self.assertEqual(result.stdout, "chr1\t50\t100\ta1\t0\t+\n")

    def test_stdin_a(self):
        with tempfile.TemporaryDirectory() as d:
            b = os.path.join(d, "b.bed")
            with open(b, "w") as f:
                f.write("chr1\t50\t150\tb1\t0\t+\n")
            result = run_mytools(
                ["intersect", "-a", "-", "-b", b],
                stdin_text="chr1\t0\t100\ta1\t0\t+\n",
            )
            self.assertEqual(result.returncode, 0)
            self.assertEqual(result.stdout, "chr1\t50\t100\ta1\t0\t+\n")

    def test_missing_file_exit_1(self):
        result = run_mytools(["intersect", "-a", "/no/such/file.bed", "-b", "/no/such/file.bed"])
        self.assertEqual(result.returncode, 1)
        self.assertEqual(result.stdout, "")

    def test_unknown_flag_exit_2(self):
        result = run_mytools(["intersect", "-a", "-", "-b", "-", "--bogus"])
        self.assertEqual(result.returncode, 2)

    def test_start_gt_end_exit_1(self):
        with tempfile.TemporaryDirectory() as d:
            a = os.path.join(d, "a.bed")
            with open(a, "w") as f:
                f.write("chr1\t100\t50\tbad\t0\t+\n")
            result = run_mytools(["intersect", "-a", a, "-b", a])
            self.assertEqual(result.returncode, 1)
            self.assertIn(":1:", result.stderr)


if __name__ == "__main__":
    unittest.main()
