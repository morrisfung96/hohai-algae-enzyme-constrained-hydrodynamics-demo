"""Rescore released frozen predictions; no scientific model or fit is run."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from lake_mc.evaluation import reproduce


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-dir',type=Path,default=ROOT/'data'/'frozen_evaluation')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'outputs'/'evaluation')
    parser.add_argument('--skip-bootstrap',action='store_true',help='Skip paired resampling intervals only; retain all point metrics and identity checks.')
    args=parser.parse_args()
    summary=reproduce(args.data_dir,args.output_dir,bootstrap=not args.skip_bootstrap)
    print(json.dumps({k:summary[k] for k in ['status','assets_verified','checks_passed','model_runs','new_fits']},indent=2))
    for record in summary['primary_mc']:
        print(f"{record['model']}: log MAE={record['MAE_log10']:.12f}; max={record['prediction_max']:.6f} ug/L; >=5 TP={record['TP']}, FP={record['FP']}")


if __name__=='__main__':
    main()
