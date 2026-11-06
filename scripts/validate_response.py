"""Validate the public response/diagnostic inputs without any model solver."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from lake_mc.response import surface


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, help='Optional JSON validation report')
    parser.add_argument('--allow-numpy-fallback', action='store_true',
                        help='Validate the fallback only; not a historical cKDTree replay')
    args = parser.parse_args()
    if surface.cKDTree is None and not args.allow_numpy_fallback:
        parser.error('SciPy cKDTree is required for the archived response regression; install the verified requirements.')
    suite = unittest.defaultTestLoader.discover(str(ROOT / 'tests'), pattern='test_response.py')
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    constraint_audit = json.loads((ROOT / 'data/response/source_constraint_audit.json').read_text())
    known_count = len(constraint_audit['exception_rows'])
    status = ('PASS_FROZEN_ARCHIVE_WITH_KNOWN_SOURCE_EXCEPTION' if known_count else 'PASS') if result.wasSuccessful() else 'FAIL'
    report = {
        'status': status,
        'tests_run': result.testsRun,
        'failures': len(result.failures),
        'errors': len(result.errors),
        'known_strict_scientific_expected_failures': len(result.expectedFailures),
        'unexpected_successes': len(result.unexpectedSuccesses),
        'scipy_backend_available': surface.cKDTree is not None,
        'source_scientific_status': constraint_audit['status'],
        'source_constraint_exceptions': known_count,
        'strict_equality_certified': False if known_count else True,
        'response_sha256': hashlib.sha256((ROOT / 'data/response/metabolic_response.csv.gz').read_bytes()).hexdigest(),
        'new_lp_runs': 0, 'new_lake_model_runs': 0, 'parameter_changes': 0,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))
    return 0 if result.wasSuccessful() else 1


if __name__ == '__main__':
    raise SystemExit(main())
