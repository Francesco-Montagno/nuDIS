#!/usr/bin/env python3
"""Add reusable random variables to an existing generator event file."""

from __future__ import annotations

import argparse

from smearing.event_file import add_standard_normal_columns

SMEARING_COLUMNS = ("z_x", "z_Q2", "z_E")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Add detector-smearing random variables to an event file"
    )
    parser.add_argument("--input", required=True, help="Input generator event file")
    parser.add_argument("--output", required=True, help="New output event file")
    parser.add_argument("--seed", type=int, default=None, help="Optional random seed")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    event_count = add_standard_normal_columns(
        input_path=args.input,
        output_path=args.output,
        column_names=SMEARING_COLUMNS,
        seed=args.seed,
    )
    print(
        f"Wrote {event_count} events with "
        f"{', '.join(SMEARING_COLUMNS)} to: {args.output}"
    )


if __name__ == "__main__":
    main()
