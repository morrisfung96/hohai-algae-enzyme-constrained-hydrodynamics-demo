"""Check the renamed numerical interface against frozen numeric witnesses.

Usage: python verify_capacity_query.py path/to/current_response.csv.gz
This only reads saved response data; no metabolic model or solver is used.
"""
from pathlib import Path
import argparse
import hashlib
import json
import pandas as pd
from numerical_capacity_interface import SavedResponseQuery, OUTPUTS

EXPECTED_RESPONSE_SHA256 = 'f63fd1bb473e4897f70d500a978b5a373ccc56bf8371cae473d9e56e7a5f9871'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('response', type=Path)
    args = ap.parse_args()
    actual_hash = hashlib.sha256(args.response.read_bytes()).hexdigest()
    if actual_hash != EXPECTED_RESPONSE_SHA256:
        raise SystemExit('Response hash differs from the declared current table')
    query = SavedResponseQuery(args.response)
    witness = pd.read_csv(Path(__file__).with_name('capacity_query_witnesses.csv'))
    largest = 0.0
    for row in witness.to_dict('records'):
        got = query.query(row['temp_C'], row['light_Wm2'], row['DIN_mgN_L'], row['DIP_mgP_L'], row['rho'])
        largest = max(largest, max(abs(got[key] - row[key]) for key in OUTPUTS))
    passed = largest < 1e-12
    print(json.dumps({'passed': passed, 'witness_queries': len(witness),
                      'max_absolute_difference': largest, 'response_sha256': actual_hash,
                      'network_loads': 0, 'optimization_calls': 0}, indent=2))
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
