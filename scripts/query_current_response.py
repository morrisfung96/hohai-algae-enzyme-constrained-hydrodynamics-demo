"""Query the current precomputed response with the named numerical interface."""
from pathlib import Path
import argparse,json,sys
ROOT=Path(__file__).resolve().parents[1]
from numerical_capacity_interface import SavedResponseQuery
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--temperature',type=float,required=True)
    p.add_argument('--light',type=float,required=True)
    p.add_argument('--din',type=float,required=True)
    p.add_argument('--dip',type=float,required=True)
    p.add_argument('--rho',type=float,required=True)
    a=p.parse_args()
    r=SavedResponseQuery(ROOT/'current_results/data/response/environment_response.csv.gz')
    print(json.dumps(r.query(a.temperature,a.light,a.din,a.dip,a.rho),indent=2))
if __name__=='__main__':main()
