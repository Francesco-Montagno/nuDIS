from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from smearing.event_file import add_standard_normal_columns


PREAMBLE = (
    "# =============================================================\n"
    "# RUN CONFIGURATION\n"
    "# process = d_p\n"
    "# PID Parton\tPID Nucleon\tx\tQ2\tE\tWeight\n"
    "\n"
)
EVENTS = (
    "-1\t2212\t0.018658759744039518\t10.161244539533463\t303.37342114147646\t3.6390910414969932e-06\n"
    "-1\t2212\t0.01801517833736843\t10.128450995459293\t302.0176785777517\t3.862170820234393e-06\n"
)


class EventFileTests(unittest.TestCase):
    def test_preserves_input_and_appends_reproducible_random_columns(self) -> None:
        with TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "events.dat"
            first = root / "first.dat"
            second = root / "second.dat"
            original = PREAMBLE + EVENTS
            source.write_text(original, encoding="utf-8")

            columns = ("z_x", "z_Q2", "z_E")
            count = add_standard_normal_columns(source, first, columns, seed=12345)
            add_standard_normal_columns(source, second, columns, seed=12345)

            self.assertEqual(count, 2)
            self.assertEqual(source.read_text(encoding="utf-8"), original)
            self.assertEqual(first.read_bytes(), second.read_bytes())

            output_lines = first.read_text(encoding="utf-8").splitlines()
            original_comments = [line for line in original.splitlines() if line.startswith("#")]
            output_comments = [line for line in output_lines if line.startswith("#")]
            self.assertEqual(output_comments[:-1], original_comments)
            self.assertEqual(
                output_comments[-1],
                "# Added post-processing columns\tz_x\tz_Q2\tz_E",
            )

            output_events = [line for line in output_lines if line and not line.startswith("#")]
            original_events = EVENTS.splitlines()
            expected = np.random.default_rng(12345).normal(size=(2, 3))
            for index, output_event in enumerate(output_events):
                fields = output_event.split("\t")
                prefix = "\t".join(fields[:-3])
                self.assertEqual(prefix, original_events[index])
                np.testing.assert_array_equal(
                    np.asarray(fields[-3:], dtype=float),
                    expected[index],
                )

    def test_rejects_in_place_processing(self) -> None:
        with TemporaryDirectory() as directory:
            source = Path(directory) / "events.dat"
            source.write_text(PREAMBLE + EVENTS, encoding="utf-8")
            with self.assertRaises(ValueError):
                add_standard_normal_columns(source, source, ("z_x",), seed=1)


if __name__ == "__main__":
    unittest.main()
