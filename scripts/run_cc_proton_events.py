#!/usr/bin/env python3
"""
Run events.py for the four CC proton processes
(d_p, s_p, ubar_p, cbar_p) in parallel, scanning three x bins per flavor:
  bin 1: x in [1e-3, 1e-2]
  bin 2: x in [1e-2, 1e-1]
  bin 3: x in [1e-1, 1.0]

Each bin is saved in its own subfolder:  events/x_<xmin>_<xmax>/

Usage:
    python scripts/run_cc_proton_events.py
    python scripts/run_cc_proton_events.py --replica 3
    python scripts/run_cc_proton_events.py --replica 3 --card card/events_card.dat
    python scripts/run_cc_proton_events.py --max_workers 6
"""
import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PROCESSES = ["d_p", "s_p", "ubar_p", "cbar_p"]
X_BINS    = [
    (0.001000, 0.001995),
    (0.001995, 0.003981),
    (0.003981, 0.007943),
    (0.007943, 0.015849),
    (0.015849, 0.031623),
    (0.031623, 0.063096),
    (0.063096, 0.125893),
    (0.125893, 0.251189),
    (0.251189, 0.501187),
    (0.501187, 1.000000)
]

REPO_ROOT     = Path(__file__).resolve().parent.parent
EVENTS_SCRIPT = REPO_ROOT / "events.py"

# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Run CC proton events in parallel")
    parser.add_argument('--replica', type=int, default=0,
                        help='PDF replica ID (default: 0)')
    parser.add_argument('--card', type=str, default='card/events_card.dat',
                        help='Run card path relative to repo root (default: card/events_card.dat)')
    parser.add_argument('--max_workers', type=int, default=2,
                        help='Max parallel jobs (default: 2)')
    return parser.parse_args()


def _fmt(v: float) -> str:
    """Format an arbitrary float into a filesystem-safe string.
    e.g. 0.001 -> '0p001', 0.001995 -> '0p001995', 1.0 -> '1'
    """
    return f"{v:.6g}".replace('.', 'p')


def run_process(process: str, replica: int, card: Path, base_folder: Path,
                x_min: float, x_max: float) -> tuple[str, str, int]:
    bin_label = f"x_{_fmt(x_min)}_{_fmt(x_max)}"
    folder    = base_folder / bin_label
    folder.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(EVENTS_SCRIPT),
        "--card",       str(card),
        "--process",    process,
        "--replica",    str(replica),
        "--folder",     str(folder),
        "--x_min_bin",  str(x_min),
        "--x_max_bin",  str(x_max),
    ]
    tag = f"{process}|{bin_label}"
    print(f"[{tag}] starting  (replica={replica})")
    result = subprocess.run(cmd, cwd=str(REPO_ROOT))
    status = "done" if result.returncode == 0 else f"FAILED (exit {result.returncode})"
    print(f"[{tag}] {status}")
    return tag, bin_label, result.returncode


def main():
    args        = parse_args()
    card        = REPO_ROOT / args.card
    base_folder = REPO_ROOT / "data" / "output" / f"replica_{args.replica}" / "events"
    base_folder.mkdir(parents=True, exist_ok=True)
    print(f"Base output   : {base_folder}")
    print(f"Card          : {card}")
    print(f"Replica       : {args.replica}")
    print(f"Processes     : {PROCESSES}")
    print(f"x bins        : {X_BINS}")
    print(f"Total jobs    : {len(PROCESSES) * len(X_BINS)}\n")

    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {
            executor.submit(run_process, p, args.replica, card, base_folder, x_min, x_max): (p, x_min, x_max)
            for p in PROCESSES
            for x_min, x_max in X_BINS
        }
        failed = []
        for future in as_completed(futures):
            tag, bin_label, code = future.result()
            if code != 0:
                failed.append(tag)

    print("\n--- Summary ---")
    if failed:
        print(f"Failed: {failed}")
        sys.exit(1)
    else:
        print("All processes completed successfully.")
        print(f"Output: {base_folder}")


if __name__ == '__main__':
    main()
