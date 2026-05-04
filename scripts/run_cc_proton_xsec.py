#!/usr/bin/env python3
"""
Run cross_section.py for the four CC proton processes
(d_p, s_p, ubar_p, cbar_p) in parallel over a given x range.

Usage:
    python scripts/run_cc_proton_xsec.py --x_min 1e-3 --x_max 1e-2 --num_bins_x 50
    python scripts/run_cc_proton_xsec.py --x_min 1e-2 --x_max 1e-1 --num_bins_x 50 --replica 3
    python scripts/run_cc_proton_xsec.py --x_min 1e-1 --x_max 1   --num_bins_x 50 --max_workers 4
"""
import argparse
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

PROCESSES   = ["d_p", "s_p", "ubar_p", "cbar_p"]
REPO_ROOT   = Path(__file__).resolve().parent.parent
XSEC_SCRIPT = REPO_ROOT / "cross_section.py"

# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(description="Run CC proton cross sections in parallel")
    parser.add_argument('--replica', type=int, default=0,
                        help='PDF replica ID (default: 0)')
    parser.add_argument('--card', type=str, default='card/cross_section_card.dat',
                        help='Run card path relative to repo root (default: card/cross_section_card.dat)')
    parser.add_argument('--max_workers', type=int, default=4,
                        help='Max parallel jobs (default: 4)')
    parser.add_argument('--x_min_bin', type=float, default=None,
                        help='Minimum x for binning (overrides x_min_bin from card)')
    parser.add_argument('--x_max_bin', type=float, default=None,
                        help='Maximum x for binning (overrides x_max_bin from card)')
    parser.add_argument('--num_bins_x', type=int, default=None,
                        help='Number of x bins (overrides num_bins_x from card)')
    parser.add_argument('--Q2_min_bin', type=float, default=None,
                        help='Minimum Q2 for binning (overrides Q2_min_bin from card)')
    parser.add_argument('--Q2_max_bin', type=float, default=None,
                        help='Maximum Q2 for binning (overrides Q2_max_bin from card)')
    parser.add_argument('--num_bins_Q2', type=int, default=None,
                        help='Number of Q2 bins (overrides num_bins_Q2 from card)')
    return parser.parse_args()


def run_process(process: str, replica: int, card: Path, folder: Path,
                x_min: float | None, x_max: float | None, num_bins_x: int | None,
                Q2_min: float | None, Q2_max: float | None, num_bins_Q2: int | None,
                ) -> tuple[str, int]:
    cmd = [
        sys.executable, str(XSEC_SCRIPT),
        "--card",    str(card),
        "--process", process,
        "--replica", str(replica),
        "--folder",  str(folder),
    ]
    if x_min is not None:
        cmd += ["--x_min_bin", str(x_min)]
    if x_max is not None:
        cmd += ["--x_max_bin", str(x_max)]
    if num_bins_x is not None:
        cmd += ["--num_bins_x", str(num_bins_x)]
    if Q2_min is not None:
        cmd += ["--Q2_min_bin", str(Q2_min)]
    if Q2_max is not None:
        cmd += ["--Q2_max_bin", str(Q2_max)]
    if num_bins_Q2 is not None:
        cmd += ["--num_bins_Q2", str(num_bins_Q2)]
    print(f"[{process}] starting  (replica={replica})")
    print(f"[{process}] cmd: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(REPO_ROOT))
    status = "done" if result.returncode == 0 else f"FAILED (exit {result.returncode})"
    print(f"[{process}] {status}")
    return process, result.returncode


def main():
    args   = parse_args()
    card   = REPO_ROOT / args.card
    folder = REPO_ROOT / "data" / "output" / f"replica_{args.replica}" / "xsec"
    folder.mkdir(parents=True, exist_ok=True)
    print(f"Output folder : {folder}")
    print(f"Card          : {card}")
    print(f"Replica       : {args.replica}")
    print(f"Processes     : {PROCESSES}")
    if args.x_min_bin is not None:
        print(f"x_min_bin     : {args.x_min_bin}")
    if args.x_max_bin is not None:
        print(f"x_max_bin     : {args.x_max_bin}")
    if args.num_bins_x is not None:
        print(f"num_bins_x    : {args.num_bins_x}")
    if args.Q2_min_bin is not None:
        print(f"Q2_min_bin    : {args.Q2_min_bin}")
    if args.Q2_max_bin is not None:
        print(f"Q2_max_bin    : {args.Q2_max_bin}")
    if args.num_bins_Q2 is not None:
        print(f"num_bins_Q2   : {args.num_bins_Q2}")
    print()

    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {
            executor.submit(run_process, p, args.replica, card, folder,
                            args.x_min_bin, args.x_max_bin, args.num_bins_x,
                            args.Q2_min_bin, args.Q2_max_bin, args.num_bins_Q2): p
            for p in PROCESSES
        }
        failed = []
        for future in as_completed(futures):
            process, code = future.result()
            if code != 0:
                failed.append(process)

    print("\n--- Summary ---")
    if failed:
        print(f"Failed: {failed}")
        sys.exit(1)
    else:
        print("All processes completed successfully.")
        print(f"Output: {folder}")


if __name__ == '__main__':
    main()
