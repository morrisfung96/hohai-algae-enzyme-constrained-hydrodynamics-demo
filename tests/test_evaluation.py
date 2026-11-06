"""Scientific regression checks against the frozen, common-sample release."""
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from lake_mc.evaluation import log_errors, reproduce, threshold_metrics
from lake_mc.evaluation.metrics import rank_auc, spearman
from lake_mc.evaluation.reproduction import verify_assets


class FrozenEvaluationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary=tempfile.TemporaryDirectory()
        cls.output=Path(cls.temporary.name)/'scores'
        cls.result=reproduce(ROOT/'data'/'frozen_evaluation',cls.output)

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    def test_exact_current_mc_metrics(self):
        result={row['source_model']:row for row in self.result['primary_mc']}
        self.assertAlmostEqual(result['Dynamic RevB']['MAE_log10'],.366237359737432,places=13)
        self.assertAlmostEqual(result['Fixed capacity RevB']['MAE_log10'],.3708119251159395,places=13)
        self.assertAlmostEqual(result['Biomass source RevB']['MAE_log10'],.36082840564210417,places=13)
        self.assertAlmostEqual(result['Dynamic q80 RevB']['MAE_log10'],.429533890252874,places=13)
        self.assertAlmostEqual(result['Training median']['MAE_log10'],.47335155422729386,places=13)
        self.assertEqual((result['Dynamic RevB']['TP'],result['Dynamic RevB']['observed_positive']),(0,36))
        self.assertEqual((result['Dynamic q80 RevB']['TP'],result['Dynamic q80 RevB']['FP']),(19,30))
        self.assertAlmostEqual(result['Dynamic RevB']['prediction_max'],4.6541237288,places=9)
        self.assertEqual(len(pd.read_csv(self.output/'archived_180_metrics_checks.csv')),180)

    def test_sample_identity_and_no_new_fitting(self):
        self.assertEqual((self.result['mc_rows'],self.result['biomass_surface_rows']),(164,64))
        self.assertEqual((self.result['biomass_development_rows'],self.result['biomass_2019_holdout_rows']),(141,40))
        self.assertEqual(self.result['stratification_rows'],124)
        self.assertEqual(self.result['mapping_interval_rows'],12)
        self.assertEqual(len(pd.read_csv(self.output/'biomass_adjacent_changes_149.csv')),149)
        self.assertEqual((self.result['model_runs'],self.result['new_fits'],self.result['new_parameter_selection']),(0,0,0))
        self.assertTrue(pd.read_csv(self.output/'identity_checks.csv').passed.all())

    def test_biomass_surface_and_transfer_values(self):
        rows=self.result['biomass']
        surface=[row for row in rows if row['n']==64]
        self.assertAlmostEqual(surface[0]['MAE_log10'],.3161447273828053,places=13)
        self.assertAlmostEqual(surface[1]['MAE_log10'],.2764397276034949,places=13)
        self.assertAlmostEqual(surface[0]['spearman'],.5640567765567764,places=13)
        self.assertAlmostEqual(surface[1]['spearman'],.721565934065934,places=13)
        metabolic=[row for row in rows if row['model'].startswith('Metabolic')]
        by_count={row['n']:row for row in metabolic}
        self.assertAlmostEqual(by_count[141]['MAE_log10'],.2511367371,places=9)
        self.assertAlmostEqual(by_count[40]['MAE_log10'],.2614030490,places=9)

    def test_fpr10_is_a_separate_frozen_risk_branch(self):
        all_scope={row['threshold_ugL']:row for row in self.result['primary_risk']}
        self.assertEqual((all_scope[5]['TP'],all_scope[5]['FP']),(11,20))
        self.assertEqual((all_scope[8]['TP'],all_scope[8]['FP']),(5,1))
        self.assertEqual((all_scope[10]['TP'],all_scope[10]['FP']),(2,1))

    def test_input_corruption_stops_scoring(self):
        with tempfile.TemporaryDirectory() as directory:
            damaged=Path(directory)/'data'
            shutil.copytree(ROOT/'data'/'frozen_evaluation',damaged)
            with (damaged/'formal_predictions_164.csv').open('ab') as stream:
                stream.write(b'\n')
            with self.assertRaisesRegex(ValueError,'hash mismatch'):
                verify_assets(damaged)

    def test_metric_definitions_ties_and_missing_rows(self):
        self.assertAlmostEqual(rank_auc([True,False],[1,1]),.5)
        self.assertAlmostEqual(spearman([1,1,2],[1,1,2]),1.)
        threshold=threshold_metrics([5,4,10],[5,5,0],5)
        self.assertEqual((threshold['TP'],threshold['FP'],threshold['FN'],threshold['TN']),(1,1,1,0))
        with self.assertRaises(ValueError):
            log_errors([1,np.nan],[1,2])
        with self.assertRaises(ValueError):
            log_errors([-1,2],[1,2])
        with self.assertRaises(ValueError):
            log_errors([0,1],[1,2],offset=0)


if __name__=='__main__':
    unittest.main()
