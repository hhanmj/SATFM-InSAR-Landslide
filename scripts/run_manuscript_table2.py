from __future__ import annotations

import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from satfm.manuscript_table2 import MC_SEEDS, run_manuscript_table2


def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce manuscript Table 2 synthetic benchmark.")
    parser.add_argument("--output", default="outputs/manuscript_table2", help="Output directory.")
    parser.add_argument("--seeds", type=int, default=MC_SEEDS, help="Number of Monte Carlo seeds.")
    args = parser.parse_args()

    _, summary = run_manuscript_table2(args.output, n_seeds=args.seeds)
    print(summary[["Method", "East RMSE (mm/yr)", "North RMSE (mm/yr)", "Up RMSE (mm/yr)"]].to_string(index=False))


if __name__ == "__main__":
    main()

