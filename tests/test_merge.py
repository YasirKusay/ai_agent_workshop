"""Unit tests for `mytools merge`. No bedtools required -- see tests/README.md
for why these exist alongside (not instead of) the golden tests in
run_golden.sh.
"""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from mytools_lib import bed, merge  # noqa: E402


def rec(chrom, start, end, name="x", score="0", strand="+", lineno=1):
    fields = [chrom, str(start), str(end), name, score, strand]
    return bed.BedRecord(chrom, start, end, fields, lineno)


class OverlapPredicateTest(unittest.TestCase):
    """BED is 0-based half-open: a.start < b.end AND b.start < a.end, strict."""

    def test_real_overlap(self):
        self.assertTrue(bed.overlaps(100, 200, 150, 250))
        self.assertTrue(bed.overlaps(150, 250, 100, 200))

    def test_nested(self):
        self.assertTrue(bed.overlaps(100, 500, 200, 300))

    def test_bookended_does_not_overlap(self):
        # a.end == b.start: touching, not overlapping.
        self.assertFalse(bed.overlaps(100, 200, 200, 300))
        self.assertFalse(bed.overlaps(200, 300, 100, 200))

    def test_gap_does_not_overlap(self):
        self.assertFalse(bed.overlaps(100, 200, 250, 300))

    def test_zero_length_never_satisfies_strict_predicate(self):
        # A zero-length feature can never have start < end on either side.
        self.assertFalse(bed.overlaps(500, 500, 500, 600))
        self.assertFalse(bed.overlaps(500, 600, 500, 500))
        self.assertFalse(bed.overlaps(0, 0, 0, 0))

    def test_position_zero(self):
        self.assertTrue(bed.overlaps(0, 100, 50, 150))
        self.assertFalse(bed.overlaps(0, 100, 100, 200))


class BookendMergeTest(unittest.TestCase):
    """Bookended intervals merge under the default -d 0, but not once a real
    overlap (negative -d) is required."""

    def test_bookended_merges_at_d0(self):
        records = [rec("chr1", 100, 200), rec("chr1", 200, 300)]
        clusters = list(merge._cluster_single(records, 0, None))
        self.assertEqual([c.bounds() for c in clusters], [(100, 300)])

    def test_bookended_does_not_merge_when_overlap_required(self):
        records = [rec("chr1", 100, 200), rec("chr1", 200, 300)]
        clusters = list(merge._cluster_single(records, -1, None))
        self.assertEqual([c.bounds() for c in clusters], [(100, 200), (200, 300)])

    def test_gap_of_one_does_not_merge_at_d0(self):
        records = [rec("chr1", 100, 200), rec("chr1", 201, 300)]
        clusters = list(merge._cluster_single(records, 0, None))
        self.assertEqual([c.bounds() for c in clusters], [(100, 200), (201, 300)])

    def test_gap_of_one_merges_at_d1(self):
        records = [rec("chr1", 100, 200), rec("chr1", 201, 300)]
        clusters = list(merge._cluster_single(records, 1, None))
        self.assertEqual([c.bounds() for c in clusters], [(100, 300)])


class ZeroLengthMergeTest(unittest.TestCase):
    """Verified against real bedtools: a zero-length feature that merges with
    a neighbour is padded 1bp on each side for the purposes of the merged
    boundary; a zero-length feature that ends up alone keeps its original
    coordinates. `bedtools merge` on sorted data/a.bed reports the
    chr1:500-500 feature merged with chr1:500-600 as chr1:499-600."""

    def test_zero_length_alone_keeps_original_coords(self):
        clusters = list(merge._cluster_single([rec("chr1", 500, 500)], 0, None))
        self.assertEqual([c.bounds() for c in clusters], [(500, 500)])

    def test_zero_length_bookended_with_neighbour_expands_left(self):
        records = [rec("chr1", 500, 500), rec("chr1", 500, 600)]
        clusters = list(merge._cluster_single(records, 0, None))
        self.assertEqual([c.bounds() for c in clusters], [(499, 600)])

    def test_zero_length_after_normal_feature_expands_right(self):
        records = [rec("chr1", 400, 500), rec("chr1", 500, 500)]
        clusters = list(merge._cluster_single(records, 0, None))
        self.assertEqual([c.bounds() for c in clusters], [(400, 501)])

    def test_zero_length_between_two_normals_merges_all_three(self):
        records = [rec("chr1", 400, 500), rec("chr1", 500, 500), rec("chr1", 500, 600)]
        clusters = list(merge._cluster_single(records, 0, None))
        self.assertEqual([c.bounds() for c in clusters], [(400, 600)])

    def test_zero_length_with_true_gap_does_not_merge(self):
        # A gap of 2 defeats even the 1bp padding.
        records = [rec("chr1", 300, 400), rec("chr1", 402, 402)]
        clusters = list(merge._cluster_single(records, 0, None))
        self.assertEqual([c.bounds() for c in clusters], [(300, 400), (402, 402)])

    def test_zero_length_with_gap_of_one_merges_via_padding(self):
        records = [rec("chr1", 300, 400), rec("chr1", 401, 401)]
        clusters = list(merge._cluster_single(records, 0, None))
        self.assertEqual([c.bounds() for c in clusters], [(300, 402)])


class NestedAndPositionZeroTest(unittest.TestCase):
    def test_fully_nested_interval_absorbed(self):
        records = [rec("chr1", 100, 500), rec("chr1", 200, 300)]
        clusters = list(merge._cluster_single(records, 0, None))
        self.assertEqual([c.bounds() for c in clusters], [(100, 500)])

    def test_interval_at_position_zero(self):
        records = [rec("chr1", 0, 100)]
        clusters = list(merge._cluster_single(records, 0, None))
        self.assertEqual([c.bounds() for c in clusters], [(0, 100)])


class ForcedStrandMergeTest(unittest.TestCase):
    """-s: same-strand features merge across an intervening different-strand
    feature that sits between them positionally (verified against real
    bedtools -- this is not the same as a single strand-blind sweep)."""

    def test_same_strand_merges_across_other_strand_feature(self):
        records = [
            rec("chr1", 100, 200, strand="+"),
            rec("chr1", 150, 250, strand="-"),
            rec("chr1", 200, 300, strand="+"),
        ]
        clusters = list(merge._cluster_forced_strand(records, 0, None))
        self.assertEqual(
            sorted(c.bounds() for c in clusters),
            [(100, 300), (150, 250)],
        )


class ColumnOpsTest(unittest.TestCase):
    def test_sum(self):
        self.assertEqual(merge._apply_op("sum", ["10", "20", "30"], ","), "60")
        self.assertEqual(merge._apply_op("sum", ["1", "2", "3"], ","), "6")

    def test_collapse_keeps_order_and_duplicates(self):
        self.assertEqual(merge._apply_op("collapse", ["b", "a", "b"], ","), "b,a,b")

    def test_distinct_dedups_and_sorts(self):
        self.assertEqual(merge._apply_op("distinct", ["b", "a", "b"], ","), "a,b")

    def test_count_and_count_distinct(self):
        self.assertEqual(merge._apply_op("count", ["a", "a", "b"], ","), "3")
        self.assertEqual(merge._apply_op("count_distinct", ["a", "a", "b"], ","), "2")

    def test_first_last(self):
        self.assertEqual(merge._apply_op("first", ["a", "b", "c"], ","), "a")
        self.assertEqual(merge._apply_op("last", ["a", "b", "c"], ","), "c")


class SortednessCheckTest(unittest.TestCase):
    def test_sorted_within_chrom_passes_through(self):
        records = [rec("chr1", 0, 10), rec("chr1", 10, 20), rec("chr2", 0, 5)]
        self.assertEqual(len(list(bed.require_sorted(records, "test"))), 3)

    def test_decreasing_start_within_chrom_errors(self):
        records = [rec("chr1", 10, 20), rec("chr1", 5, 15)]
        with self.assertRaises(bed.BedDataError):
            list(bed.require_sorted(records, "test"))

    def test_interleaved_chrom_errors(self):
        records = [rec("chr1", 0, 10), rec("chr2", 0, 10), rec("chr1", 20, 30)]
        with self.assertRaises(bed.BedDataError):
            list(bed.require_sorted(records, "test"))

    def test_any_chrom_order_is_fine_as_long_as_contiguous(self):
        records = [rec("chrX", 0, 10), rec("chr1", 0, 10)]
        self.assertEqual(len(list(bed.require_sorted(records, "test"))), 2)


class MalformedLineTest(unittest.TestCase):
    def _run(self, content):
        with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as f:
            f.write(content)
            path = f.name
        try:
            return list(bed.iter_records(path))
        finally:
            os.unlink(path)

    def test_too_few_columns_errors(self):
        with self.assertRaises(bed.BedDataError):
            self._run("chr1\t100\n")

    def test_non_numeric_start_errors(self):
        with self.assertRaises(bed.BedDataError):
            self._run("chr1\tfoo\t200\n")

    def test_start_greater_than_end_errors(self):
        with self.assertRaises(bed.BedDataError):
            self._run("chr1\t200\t100\n")

    def test_negative_start_errors(self):
        with self.assertRaises(bed.BedDataError):
            self._run("chr1\t-5\t100\n")

    def test_track_browser_comment_and_blank_lines_are_skipped(self):
        records = self._run("track name=foo\nbrowser position chr1\n#c\n\nchr1\t0\t100\n")
        self.assertEqual(len(records), 1)
        self.assertEqual((records[0].chrom, records[0].start, records[0].end), ("chr1", 0, 100))


class UsageErrorTest(unittest.TestCase):
    def test_unknown_flag_is_usage_error(self):
        with self.assertRaises(bed.UsageError):
            merge.parse_args(["-i", "x.bed", "--bogus"])

    def test_c_without_o_is_usage_error(self):
        opts = merge.parse_args(["-i", "x.bed", "-c", "4"])
        with self.assertRaises(bed.UsageError):
            merge._resolve_cols_ops(opts)

    def test_mismatched_c_and_o_counts_is_usage_error(self):
        opts = merge.parse_args(["-i", "x.bed", "-c", "4,5,6", "-o", "sum,mean"])
        with self.assertRaises(bed.UsageError):
            merge._resolve_cols_ops(opts)

    def test_invalid_operation_name_is_usage_error(self):
        opts = merge.parse_args(["-i", "x.bed", "-c", "4", "-o", "bogus"])
        with self.assertRaises(bed.UsageError):
            merge._resolve_cols_ops(opts)

    def test_invalid_strand_value_is_usage_error(self):
        with self.assertRaises(bed.UsageError):
            merge.parse_args(["-i", "x.bed", "-S", "x"])


if __name__ == "__main__":
    unittest.main()
