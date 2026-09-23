"""Post-process generator event files without changing the generator itself."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path
from typing import Sequence

import numpy as np


def _event_line_indices(lines: Sequence[str]) -> list[int]:
    """Return indices of non-empty, non-comment event rows."""
    return [
        index
        for index, line in enumerate(lines)
        if line.strip() and not line.lstrip().startswith("#")
    ]


def add_standard_normal_columns(
    input_path: str | Path,
    output_path: str | Path,
    column_names: Sequence[str],
    seed: int | None = None,
) -> int:
    """Append independent standard-normal columns to a generator event file.

    Existing lines and values are copied verbatim. A tab and the generated
    value(s) are appended only to event rows. The input and output paths must
    be different.

    Returns the number of processed events.
    """
    source = Path(input_path)
    destination = Path(output_path)

    if source.resolve() == destination.resolve():
        raise ValueError("Input and output must be different files")
    if not column_names:
        raise ValueError("At least one output column is required")
    if len(set(column_names)) != len(column_names):
        raise ValueError("Output column names must be unique")

    lines = source.read_text(encoding="utf-8").splitlines(keepends=True)
    event_indices = _event_line_indices(lines)

    rng = np.random.default_rng(seed)
    random_values = rng.normal(
        loc=0.0,
        scale=1.0,
        size=(len(event_indices), len(column_names)),
    )

    # Keep every original comment unchanged and add a separate description of
    # the post-processing columns immediately before the first event.
    column_note = "# Added post-processing columns\t" + "\t".join(column_names) + "\n"
    insertion_index = event_indices[0] if event_indices else len(lines)
    lines.insert(insertion_index, column_note)

    # Account for the inserted metadata line when updating event rows.
    for row_number, original_index in enumerate(event_indices):
        index = original_index + 1
        original = lines[index]
        newline = "\r\n" if original.endswith("\r\n") else "\n" if original.endswith("\n") else ""
        content = original[: -len(newline)] if newline else original
        generated = "\t".join(repr(float(value)) for value in random_values[row_number])
        lines[index] = f"{content}\t{generated}{newline}"

    destination.parent.mkdir(parents=True, exist_ok=True)
    file_descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        text=True,
    )
    try:
        with os.fdopen(file_descriptor, "w", encoding="utf-8", newline="") as stream:
            stream.writelines(lines)
        os.replace(temporary_name, destination)
    except BaseException:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise

    return len(event_indices)
