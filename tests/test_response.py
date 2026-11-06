"""Scientific regression of the archived response interface and data identity."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from lake_mc.response import ResponseConfig, ResponseSurfaceV16, RecoveredCacheArchive
from lake_mc.response import surface as surface_module

DATA = ROOT / 'data/response'


class ResponseRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.table = pd.read_csv(DATA / 'metabolic_response.csv.gz')
        cls.config = ResponseConfig(require_native_v16_response=True)

    def test_formal_identity_units_and_original_acceptance_guard(self):
        self.assertEqual(hashlib.sha256((DATA / 'metabolic_response.csv.gz').read_bytes()).hexdigest(),
                         'c68f04181b4825ce6dae1781d5c7b88e426d2945f1566fdd50468ede1942178d')
        self.assertEqual(self.table.shape, (7555, 16))
        self.assertEqual(set(self.table.toxin_solution_status), {'optimal'})
        self.assertEqual(set(self.table.growth_constraint_mode), {'equality'})
        np.testing.assert_allclose(self.table.mu_toxin_solution_h,
            self.table.biomass_flux_toxin_solution_raw * 1000, rtol=1e-12, atol=1e-14)
        # Reproduce the producer's pre-existing raw-unit acceptance guard.
        # The stricter mathematical check below remains a known failure.
        report = json.loads((DATA / 'source_constraint_audit.json').read_text())
        bad = self.table[self.table.growth_constraint_residual_raw.abs() > report['raw_residual_tolerance']]
        self.assertEqual(list(bad.env_id), [r['env_id'] for r in report['exception_rows']])
        np.testing.assert_allclose(bad.growth_constraint_residual_raw,
            [r['growth_constraint_residual_raw'] for r in report['exception_rows']], rtol=1e-14, atol=1e-20)
        self.assertLess(abs(self.table.growth_constraint_residual_raw).max(), 2e-6)

    @unittest.expectedFailure
    def test_known_source_strict_growth_equality_exception(self):
        # This is intentionally an expected failure, not a looser strict test.
        # A source correction must revise the frozen identity and issue record.
        np.testing.assert_allclose(self.table.mu_committed_h,
            self.table.mu_max_h * self.table.rho_growth_commitment, rtol=1e-12, atol=1e-14)
        self.assertLess(abs(self.table.growth_constraint_residual_raw).max(), 1e-10)

    def test_night_gate_has_no_growth_capacity_or_cache_entry(self):
        rs = ResponseSurfaceV16(DATA / 'metabolic_response.csv.gz', self.config)
        for light in [-100, 0, 1e-9]:
            response = rs.lookup(temp_C=24.5, light_Wm2=light, DIN_mgN_L=.3,
                DIP_mgP_L=.003, rho_growth_commitment=.7)
            self.assertEqual(response['matched_env_ids'], 'EXACT_NIGHT_GATE')
            for key in ['mu_max_h', 'mu_committed_h', 'toxin_capacity_upper_per_gDW_h',
                        'mu_toxin_solution_h', 'growth_constraint_target_raw']:
                self.assertEqual(response[key], 0.)
        self.assertEqual(rs.lookup_cache, {})

    def test_original_class_frozen_regression(self):
        fixture = json.loads((DATA / 'query_regression_cases.json').read_text())
        for backend in ['scipy', 'numpy']:
            if backend == 'scipy' and surface_module.cKDTree is None:
                continue
            rs = ResponseSurfaceV16(DATA / 'metabolic_response.csv.gz', self.config)
            if backend == 'numpy':
                rs.tree = None
            for case in fixture['cases']:
                if case['backend'] != backend:
                    continue
                query = {k:v for k,v in case.items() if k not in ['backend','case','output']}
                actual = rs.lookup(**query)
                self.assertEqual(actual.keys(), case['output'].keys())
                for key, expected in case['output'].items():
                    with self.subTest(backend=backend, case=case['case'], field=key):
                        if isinstance(expected, str):
                            self.assertEqual(actual[key], expected)
                        else:
                            self.assertAlmostEqual(actual[key], expected,
                                delta=max(1e-14, abs(expected)*1e-12))

    def test_rho_clip_and_capacity_range(self):
        rs = ResponseSurfaceV16(DATA / 'metabolic_response.csv.gz', self.config)
        query = dict(temp_C=24.5, light_Wm2=100, DIN_mgN_L=.3, DIP_mgP_L=.003)
        self.assertEqual(rs.lookup(**query, rho_growth_commitment=-1),
                         rs.lookup(**query, rho_growth_commitment=0))
        self.assertEqual(rs.lookup(**query, rho_growth_commitment=2),
                         rs.lookup(**query, rho_growth_commitment=1))
        for rho in [0,.25,.5,.75,1]:
            response = rs.lookup(**query, rho_growth_commitment=rho)
            self.assertGreaterEqual(response['toxin_capacity_upper_per_gDW_h'], 0)
            self.assertLessEqual(response['toxin_capacity_upper_per_gDW_h'],
                                self.table.toxin_capacity_upper_per_gDW_h.max())

    def test_rounded_cache_retains_first_query(self):
        q1 = dict(temp_C=24.5, light_Wm2=100, DIN_mgN_L=.3, DIP_mgP_L=.003,
                  rho_growth_commitment=.9501)
        q2 = dict(q1, rho_growth_commitment=.9504)
        rs = ResponseSurfaceV16(DATA / 'metabolic_response.csv.gz', self.config)
        a = rs.lookup(**q1)
        self.assertIs(rs.lookup(**q2), a)
        reverse = ResponseSurfaceV16(DATA / 'metabolic_response.csv.gz', self.config)
        b = reverse.lookup(**q2)
        self.assertIs(reverse.lookup(**q1), b)
        self.assertNotEqual(a['toxin_capacity_upper_per_gDW_h'],
                            b['toxin_capacity_upper_per_gDW_h'])

    def test_metabolic_diagnostic_identity(self):
        directory = ROOT / 'data/metabolic_diagnostics'
        manifest = json.loads((directory / 'manifest.json').read_text())
        for entry in manifest['files']:
            p = directory / entry['file']
            self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), entry['sha256'])
            self.assertEqual(len(pd.read_csv(p)), entry['rows'])
            self.assertNotIn(':\\', p.read_text(encoding='utf-8-sig'))

    def test_recovered_cache_identity_coverage_and_strict_context(self):
        archive = RecoveredCacheArchive(DATA / 'recovered_worker_cache_2022.csv.gz')
        self.assertEqual(len(archive.grids), 123)
        self.assertEqual(len(archive.df), 88191)
        audit = pd.read_csv(DATA / 'recovered_worker_cache_2022_forward_audit.csv')
        self.assertTrue(audit.strict_pass_1e10.all())
        self.assertLess(audit.max_abs_state_process_difference.max(), 1e-10)
        surface = ResponseSurfaceV16(DATA / 'metabolic_response.csv.gz', self.config)
        for gid in [archive.grids[0], archive.grids[len(archive.grids)//2], archive.grids[-1]]:
            context = archive.for_grid(surface, gid)
            row = archive.df[archive.df.grid_id.eq(gid)].iloc[0]
            response = context.lookup(temp_C=row.key_temp_C, light_Wm2=row.key_light_Wm2,
                DIN_mgN_L=row.key_DIN_mgN_L, DIP_mgP_L=row.key_DIP_mgP_L,
                rho_growth_commitment=row.key_rho, k_env=int(row.key_k_env))
            self.assertEqual(response['toxin_capacity_upper_per_gDW_h'], row.toxin_capacity_upper_per_gDW_h)
            night = context.lookup(temp_C=24.5, light_Wm2=0, DIN_mgN_L=.3,
                DIP_mgP_L=.003, rho_growth_commitment=.7)
            self.assertEqual(night['matched_env_ids'], 'EXACT_NIGHT_GATE')
            with self.assertRaisesRegex(ValueError, 'left the retained'):
                context.lookup(temp_C=-999, light_Wm2=100, DIN_mgN_L=.3,
                    DIP_mgP_L=.003, rho_growth_commitment=.7)
        with self.assertRaisesRegex(ValueError, 'outside'):
            archive.for_grid(surface, 'NOT_RETAINED')


if __name__ == '__main__':
    unittest.main()
