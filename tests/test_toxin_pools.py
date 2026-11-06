"""Numerical regressions from retained states; these are not whole-lake tests."""
from pathlib import Path
import json
import unittest
import numpy as np

from lake_mc.toxin_pools import (
    BIOMASS_SCALE_GDW_M3_PER_STATE,
    InventoryContractError,
    InventoryHandoff,
    release_loss_gdw,
    update_variant,
)

FIXTURE = Path(__file__).resolve().parents[1] / 'data/replay_example'


class ToxinPoolRegression(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.meta = json.loads((FIXTURE / 'real_source_pool_steps.json').read_text('utf-8'))
        cls.data = np.load(FIXTURE / 'real_source_pool_steps.npz', allow_pickle=False)
        cls.fields = {str(n): j for j, n in enumerate(cls.data['state_fields'])}

    @classmethod
    def tearDownClass(cls):
        cls.data.close()

    def test_repaired_original_arithmetic_on_real_frozen_inputs(self):
        """All variants, four years, and both release paths retain their values."""
        self.assertEqual(set(c['year'] for c in self.meta['cases']), {2019, 2020, 2021, 2022})
        self.assertEqual(len(set(c['variant'] for c in self.meta['cases'])), 3)
        self.assertGreater(sum(c['fallback_positive_cells'] for c in self.meta['cases']), 0)
        scale = float(self.data['biomass_scale'])
        self.assertEqual(scale, BIOMASS_SCALE_GDW_M3_PER_STATE)
        for case in self.meta['cases']:
            with self.subTest(case=case['key']):
                key = case['key']
                pools = tuple(self.data[key + '_pools'])
                nxt, diag = update_variant(
                    self.data[key + '_state'], self.fields, 0, pools,
                    case['dt_h'], case['variant'], case['constants'],
                    biomass_scale_gDW_m3_per_state=scale,
                )
                np.testing.assert_array_equal(np.stack(nxt), self.data[key + '_next'])
                for name, values in diag.items():
                    np.testing.assert_array_equal(values, self.data[key + '_diag_' + name])
                for values in nxt:
                    self.assertTrue(np.isfinite(values).all())
                    self.assertTrue((values >= 0).all())

    def test_total_pool_balance_on_real_steps(self):
        for case in self.meta['cases']:
            key = case['key']
            before = self.data[key + '_pools'].sum(axis=0)
            after = self.data[key + '_next'].sum(axis=0)
            inputs = self.data[key + '_diag_total_input_flux'] * case['dt_h']
            # Independently sum the three states and the reported residual.
            residual = self.data[key + '_diag_pool_mass_balance_residual']
            self.assertLess(float(np.max(np.abs(residual))), 1e-10)
            self.assertTrue((after <= before + inputs + 1e-12).all())

    def test_loss_fallback_uses_dry_weight_units(self):
        scale = BIOMASS_SCALE_GDW_M3_PER_STATE
        actual = release_loss_gdw(
            np.array([2.0, 0.0, 0.0]), np.array([10.0, 0.2, 0.3]),
            np.array([0.1, 0.25, 0.0]), scale,
        )
        np.testing.assert_array_equal(actual, [2.0, 0.2 * 0.25 * scale, 0.0])
        for invalid in [0, -1, np.nan, np.inf]:
            with self.assertRaises(ValueError):
                release_loss_gdw(np.zeros(1), np.ones(1), np.ones(1), invalid)

    def test_inventory_handoff_real_distinct_volumes_and_order(self):
        ids = tuple(map(str, self.data['inventory_ids']))
        src = self.data['source_volume_m3']
        receiver = self.data['receiver_volume_m3']
        self.assertGreater(float(np.max(np.abs(src - receiver))), 0.0)
        # Reordering receiver IDs checks that inventory follows grid identity.
        contract = InventoryHandoff(ids, ids[::-1], src, receiver[::-1])
        concentration = np.arange(1.0, len(ids) + 1.0)
        inventory = contract.to_inventory(concentration)
        np.testing.assert_array_equal(inventory, (concentration * src)[::-1])
        np.testing.assert_allclose(
            contract.to_concentration(concentration) * receiver[::-1],
            inventory, rtol=2e-16, atol=0.0,
        )
        with self.assertRaises(InventoryContractError):
            InventoryHandoff(ids, ids, np.zeros_like(src), receiver)


if __name__ == '__main__':
    unittest.main()
