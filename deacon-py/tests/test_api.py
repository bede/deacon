import base64
import inspect
import json
import tempfile
import unittest
from pathlib import Path

from deacon import Index


ROOT = Path(__file__).resolve().parents[2]
READS = ROOT / "tests/data/test_small_1.fastq.gz"
# Exact index built from READS with the default k=31, w=15 parameters.
INDEX_BYTES = base64.b64decode(
    "Ax8PQQAAAAAAAACve33W9VXXBESYNVc0k/ANrPWbvtZvcxLn2gd9//W+G89JPdM1xNcVtW1Qym/XCivXZ11fdU0UFa1UtVGcmykSrVPjVL3GpxEJZtxZ51QZB0QNr3t91vUVfff33d+Xax+t1xAXVm6mCL0XtVHXF1ED1V1Xd13ddS2x4wlm3FnnFFB1vRe1UdcXbdbLWX/+vRdUbqZIhFlzBUXkvsPlUP8G1fVe1EZdXwR91TVRVG6mCH0RNbzu9VkXlPLbtcLaMw1d70Vt1PVFFNT3r+8bfb8CpNA3lfA6TQ/1VddEUbmZIm3WkELfVMIrReVmikSYNRe9bdaQQt9UAjmjHEaev7UttVLVRnFupgg94y7CTJzvNBSVmykSYdYcvdP3Xd8zfR/QtJiVr9nKGm3U9UXU8LoXkedvbRuU8htRw+ten3V9FdS/YW/QtJgVKb9dK6w90xR0DfF1ddfVHd3fd39frn0bLWbla7ay1i9Z11ddE0V9FY997/R91/cMXddGXWN1rRTx+27ff33PNb1B02JWvmYrfde1UddYXSs1VtdKVRvFOdC9bdaQQt8U1bVTFFRd7wVEUbmZIhFmDdU1UdRXXes1FVGjUF07RQHXeg1xYeVmCn3X90zfd38fTfVQMNVDQjPVXVd3XZt1PT3T1K1F5L4Dd13ddW3W9R3dWkTuO1wONVHfv75v9P0K"
)
TEST_DIRECTORY = None
INDEX = None


def setUpModule():
    global TEST_DIRECTORY, INDEX
    TEST_DIRECTORY = tempfile.TemporaryDirectory()
    INDEX = Path(TEST_DIRECTORY.name) / "test.idx"
    INDEX.write_bytes(INDEX_BYTES)


def tearDownModule():
    TEST_DIRECTORY.cleanup()


class SignatureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index = Index(INDEX)

    def test_signatures_are_future_proof(self):
        self.assertEqual(
            str(inspect.signature(Index)),
            "(path, /, *, complexity_threshold=None)",
        )
        self.assertEqual(
            str(inspect.signature(Index.fetch)),
            "(*, name='panhuman-1', k=31, w=15, output=None, complexity_threshold=None)",
        )

        parameters = inspect.signature(Index.filter).parameters
        self.assertEqual(parameters["input"].kind, inspect.Parameter.POSITIONAL_ONLY)
        self.assertEqual(parameters["cbq_block_size"].default, 16)
        self.assertIsNone(parameters["inverse_output"].default)
        self.assertIsNone(parameters["inverse_output2"].default)
        for name, parameter in parameters.items():
            if name not in {"self", "input"}:
                self.assertEqual(parameter.kind, inspect.Parameter.KEYWORD_ONLY, name)

    def test_optional_positional_arguments_are_rejected(self):
        with self.assertRaises(TypeError):
            Index(INDEX, None)
        with self.assertRaises(TypeError):
            Index.fetch("panhuman-1")
        with self.assertRaises(TypeError):
            self.index.filter(READS, None)

    def test_primary_paths_are_positional_only(self):
        with self.assertRaises(TypeError):
            Index(path=INDEX)
        with self.assertRaises(TypeError):
            self.index.filter(input=READS)


class FilteringTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.index = Index(INDEX)

    def test_pathlike_output_and_json_summary(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            output = directory / "output.fastq"
            summary_path = directory / "summary.json"

            summary = self.index.filter(
                READS,
                output=output,
                summary=summary_path,
                threads=1,
            )

            self.assertTrue(output.is_file())
            self.assertEqual(json.loads(summary_path.read_text()), summary)
            self.assertEqual(summary["check_pairs"], False)
            self.assertGreater(summary["seqs_in"], 0)

    def test_inverse_output(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            primary = directory / "primary.fastq"
            inverse = directory / "inverse.fastq"
            summary = self.index.filter(
                READS,
                output=primary,
                inverse_output=inverse,
                threads=1,
            )

            self.assertTrue(primary.is_file())
            self.assertTrue(inverse.is_file())
            self.assertEqual(summary["inverse_output"], str(inverse))
            self.assertEqual(
                summary["seqs_out"] + summary["seqs_removed"],
                summary["seqs_in"],
            )
            records = len(primary.read_text().splitlines()) // 4
            inverse_records = len(inverse.read_text().splitlines()) // 4
            self.assertEqual(records + inverse_records, summary["seqs_in"])

    def test_check_pairs(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mate1 = directory / "mate1.fastq"
            mate2 = directory / "mate2.fastq"
            mate1.write_text("@read/1\nACGT\n+\nIIII\n")
            mate2.write_text("@read/2\nACGT\n+\nIIII\n")
            summary = self.index.filter(
                mate1,
                input2=mate2,
                output=directory / "output1.fastq",
                output2=directory / "output2.fastq",
                check_pairs=True,
                threads=1,
            )
            self.assertTrue(summary["check_pairs"])

            mate1.write_text("@read-a/1\nACGT\n+\nIIII\n")
            mate2.write_text("@read-b/2\nACGT\n+\nIIII\n")
            with self.assertRaisesRegex(RuntimeError, "Paired record name mismatch"):
                self.index.filter(
                    mate1,
                    input2=mate2,
                    output=directory / "bad1.fastq",
                    output2=directory / "bad2.fastq",
                    check_pairs=True,
                    threads=1,
                )

    def test_cbq_options_and_round_trip(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            cbq = directory / "output.cbq"
            first = self.index.filter(
                READS,
                output=cbq,
                cbq_block_size=1,
                threads=1,
            )
            second = self.index.filter(
                cbq,
                output=directory / "roundtrip.fastq",
                threads=1,
            )
            self.assertEqual(second["seqs_in"], first["seqs_out"])

            cba = directory / "output.cba"
            self.index.filter(READS, output=cba, cbq_block_size=1, threads=1)
            fasta = directory / "quality-free.fasta"
            self.index.filter(cba, output=fasta, threads=1)
            if fasta.stat().st_size:
                self.assertEqual(fasta.read_bytes()[:1], b">")

            with self.assertRaisesRegex(RuntimeError, "CBQ input does not support INPUT2"):
                self.index.filter(cbq, input2=READS, threads=1)
            with self.assertRaisesRegex(RuntimeError, "CBQ output does not support OUTPUT2"):
                self.index.filter(
                    READS,
                    output=directory / "invalid.cbq",
                    output2=directory / "invalid.fastq",
                    threads=1,
                )

    def test_argument_validation(self):
        with self.assertRaisesRegex(ValueError, "abs_threshold"):
            self.index.filter(READS, abs_threshold=0)
        with self.assertRaisesRegex(ValueError, "cbq_block_size"):
            self.index.filter(READS, cbq_block_size=0)
        with self.assertRaisesRegex(ValueError, "cbq_block_size"):
            self.index.filter(READS, cbq_block_size=1025)
        with self.assertRaisesRegex(ValueError, "input2"):
            self.index.filter(READS, input2=READS, interleaved=True)
        with self.assertRaisesRegex(RuntimeError, "check-pairs requires paired input"):
            self.index.filter(READS, check_pairs=True)


if __name__ == "__main__":
    unittest.main()
