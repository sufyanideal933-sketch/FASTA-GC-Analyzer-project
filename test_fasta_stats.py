import math
import tempfile
import unittest
from pathlib import Path

from fasta_stats import read_fasta, gc_content, build_table


class FastaAnalyzerTests(unittest.TestCase):
    def write_temp(self, text):
        handle = tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8")
        handle.write(text)
        handle.close()
        self.addCleanup(lambda: Path(handle.name).unlink(missing_ok=True))
        return handle.name

    def test_multiline_and_lowercase(self):
        path = self.write_temp(">seq1 description\nacg\nTTA\n")
        self.assertEqual(read_fasta(path), [("seq1", "ACGTTA")])

    def test_gc_ignores_ambiguous_bases(self):
        gc, denominator = gc_content("ACGTRYN-")
        self.assertAlmostEqual(gc, 50.0)
        self.assertEqual(denominator, 4)

    def test_only_ambiguous_bases_has_undefined_gc(self):
        gc, denominator = gc_content("NNRY-.")
        self.assertTrue(math.isnan(gc))
        self.assertEqual(denominator, 0)

    def test_rejects_invalid_character(self):
        path = self.write_temp(">seq1\nACGTX\n")
        with self.assertRaisesRegex(ValueError, "unsupported character"):
            read_fasta(path)

    def test_rejects_duplicate_ids(self):
        path = self.write_temp(">seq1\nACGT\n>seq1 duplicate\nGGCC\n")
        with self.assertRaisesRegex(ValueError, "Duplicate FASTA"):
            read_fasta(path)

    def test_rejects_empty_sequence(self):
        path = self.write_temp(">seq1\n>seq2\nACGT\n")
        with self.assertRaisesRegex(ValueError, "empty sequence"):
            read_fasta(path)

    def test_table_has_counts(self):
        df = build_table([("seq1", "ACGTRN-.")])
        row = df.iloc[0]
        self.assertEqual(row["length"], 8)
        self.assertEqual(row["A_count"], 1)
        self.assertEqual(row["C_count"], 1)
        self.assertEqual(row["G_count"], 1)
        self.assertEqual(row["T_count"], 1)
        self.assertEqual(row["ambiguous_count"], 2)
        self.assertEqual(row["gap_count"], 2)


if __name__ == "__main__":
    unittest.main()
