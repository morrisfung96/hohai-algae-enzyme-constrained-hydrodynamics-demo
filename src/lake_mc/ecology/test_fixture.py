"""Regression against a retained-driver 16-hour arithmetic fixture.

This checks the extracted update, not a fresh full-lake historical replay.
"""
from pathlib import Path
import json
import sys
import unittest

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[1]))
from lake_mc.ecology.config import V16Config
from lake_mc.ecology.state import simulate_grid_coupled_generator
from lake_mc.ecology.production import project_two_species_conservative_production


class EcologyArithmeticFixture(unittest.TestCase):
    def test_retained_environmental_driver_update(self):
        fixture=HERE/'fixtures'
        driver=pd.read_csv(fixture/'driver.csv')
        initial=json.loads((fixture/'initial.json').read_text())
        cfg=V16Config(**json.loads((HERE/'frozen_parameters.json').read_text()))
        cfg.max_grids=None;cfg.n_workers=1;cfg.progress_every_grids=0
        queries=json.loads((fixture/'lookup_queries.json').read_text())
        expected=pd.read_csv(fixture/'expected_state.csv')
        class FixedComputedResponse:
            def __init__(self):self.index=0
            def lookup(inner,**kwargs):
                record=queries[inner.index];inner.index+=1
                self.assertEqual(set(kwargs),set(record['query']))
                np.testing.assert_allclose([kwargs[k] for k in sorted(kwargs)],
                                           [record['query'][k] for k in sorted(kwargs)],rtol=0,atol=1e-13)
                return record['response']
        response=FixedComputedResponse()
        actual=[]
        for payload in simulate_grid_coupled_generator(driver,initial,response,cfg,[],source_state_scenario='OE1_envpow0p75_existing'):
            row=payload['source_row'];row.update(Bm_next=payload['Bm_next'],Bo_next=payload['Bo_next'])
            actual.append({column:float(row[column]) for column in expected})
        self.assertEqual(len(actual),16)
        self.assertEqual(response.index,len(queries))
        np.testing.assert_allclose(pd.DataFrame(actual),expected,rtol=0,atol=1e-14)
        self.assertGreater(expected.micro_growth_state_h.max(),0)
        self.assertTrue((expected.micro_growth_state_h==0).any())
        self.assertLess(expected.state_balance_residual.abs().max(),1e-14)

    def test_production_capacity_projection_exposes_infeasibility(self):
        feasible=project_two_species_conservative_production(B_M=1,B_O=2,P_total_native=3,w_M_target=.9,C_M=1,C_O=3)
        self.assertEqual(feasible.P_M_projected,1)
        self.assertEqual(feasible.P_O_projected,2)
        self.assertEqual(feasible.production_residual,0)
        impossible=project_two_species_conservative_production(B_M=1,B_O=2,P_total_native=5,w_M_target=.5,C_M=1,C_O=2)
        self.assertTrue(impossible.capacity_budget_inconsistent)
        self.assertEqual(impossible.production_residual,-2)
        self.assertEqual(impossible.unallocatable_budget,2)


if __name__=='__main__':
    unittest.main()
