"""Replay formal three-pool dynamics from a complete frozen yearly input pack."""
from pathlib import Path
import argparse
import json

from lake_mc.toxin_pools.replay import VARIANTS, run_replay


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs', type=Path, required=True, help='Year directory extracted from the frozen lake input archive')
    parser.add_argument('--output', type=Path, required=True, help='Output directory')
    parser.add_argument('--variant', choices=list(VARIANTS), default='dynamic')
    parser.add_argument('--max-hours', type=int, default=None, help='Optional initial diagnostic slice; omitted for the full frozen window')
    args = parser.parse_args()
    audit = run_replay(args.inputs, args.output, args.variant, max_hours=args.max_hours)
    print(json.dumps(audit, indent=2))


if __name__ == '__main__':
    main()
