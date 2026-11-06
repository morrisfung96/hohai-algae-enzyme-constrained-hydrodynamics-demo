"""Reproduce archived statistics from released, frozen input tables.

The module evaluates fixed predictions and checks fixed mapping coefficients.
It does not fit coefficients, regenerate their bootstrap intervals, run the
lake model, execute a metabolic LP, or select a risk branch.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .metrics import concentration_metrics, log_errors, spearman, threshold_metrics

KEYS = ['sample_id', 'year', 'window_id', 'station_name', 'date']
OFFSET = .01
THRESHOLDS = [5, 8, 10, 12]
READABLE_MODELS = {
    'Dynamic RevB': 'Dynamic metabolic capacity + depth-retention surface expression (LOYO OLS)',
    'Fixed capacity RevB': 'Fixed capacity + depth-retention surface expression (LOYO OLS)',
    'Biomass source RevB': 'Biomass-proportional source + depth-retention surface expression (LOYO OLS)',
    'Dynamic q80 RevB': 'Dynamic metabolic capacity + depth-retention surface expression (LOYO q0.80)',
    'Model biomass': 'Model biomass regression (LOYO)',
    'Model biomass and composition': 'Model biomass and composition regression (LOYO)',
    'Training median': 'Training-year median baseline (LOYO)',
}


def _clean_json(value):
    if isinstance(value, dict):
        return {k: _clean_json(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_clean_json(v) for v in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (float, np.floating)) and not np.isfinite(value):
        return None
    if isinstance(value, np.floating):
        return float(value)
    return value


def _check_equal(a, b, label, checks, atol=1e-12):
    a = np.asarray(a); b = np.asarray(b)
    if a.shape != b.shape or not np.allclose(a, b, rtol=0, atol=atol, equal_nan=True):
        raise ValueError('Evaluation identity failed: ' + label)
    checks.append(dict(check=label, passed=True))


def _scopes(frame):
    yield 'all', frame
    for year, group in frame.groupby('year', sort=True):
        yield str(year), group


def _write(frame, path):
    frame.to_csv(path, index=False, encoding='utf-8')


def verify_assets(data_dir):
    """Verify every published frozen CSV hash and row count before scoring."""
    data_dir = Path(data_dir)
    provenance = json.loads((data_dir / 'provenance.json').read_text(encoding='utf-8'))
    records = []
    for entry in provenance['records']:
        path = data_dir / entry['file']
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if digest != entry['sha256']:
            raise ValueError('Frozen input hash mismatch: ' + entry['file'])
        frame = pd.read_csv(path)
        if len(frame) != entry['rows'] or frame.columns.tolist() != entry['columns']:
            raise ValueError('Frozen input schema or coverage changed: ' + entry['file'])
        records.append(dict(file=entry['file'], sha256=digest, rows=len(frame)))
    return records


def _mapping_checks(formal, parameters, intervals, checks):
    rows = []
    for row in parameters.itertuples(index=False):
        stem = row.scenario.replace(' ', '_')
        test = formal[formal.year.eq(row.holdout_year)]
        training = formal[~formal.year.eq(row.holdout_year)]
        if len(test) != row.n_test or len(training) != row.n_train:
            raise ValueError('LOYO train/test coverage mismatch')
        predictor = test[stem + '_process_ugL']
        predicted = np.maximum(10 ** (row.intercept + row.slope * np.log10(predictor + OFFSET)) - OFFSET, 0)
        archived = test[f'{row.scenario}_{row.mapping}_prediction_ugL']
        _check_equal(predicted, archived, f'Fixed coefficients: {row.scenario}, {row.mapping}, {row.holdout_year}', checks)
        rows.append(dict(model=READABLE_MODELS.get(row.scenario, row.scenario), source_scenario=row.scenario,
                         mapping=row.mapping, holdout_year=int(row.holdout_year), n_train=len(training), n_test=len(test),
                         train_key_hash=hashlib.sha256('\n'.join(sorted(training.index)).encode()).hexdigest(),
                         test_key_hash=hashlib.sha256('\n'.join(sorted(test.index)).encode()).hexdigest(),
                         max_prediction_difference=float(np.max(abs(predicted-archived)))))
    if len(intervals) != 12:
        raise ValueError('Current coefficient interval archive must have exactly 12 rows.')
    for row in intervals.itertuples(index=False):
        p = parameters[(parameters.scenario == row.input) & (parameters.mapping == 'refit_ols')
                       & (parameters.holdout_year == row.holdout_year)]
        if len(p) != 1 or int(p.iloc[0].n_train) != row.n_train:
            raise ValueError('S24 point coefficient identity mismatch')
        _check_equal([row.intercept, row.slope], p.iloc[0][['intercept','slope']].to_numpy(float),
                     f'S24 point coefficients: {row.input}, {row.holdout_year}', checks)
    return pd.DataFrame(rows)


def _formal_reference(formal, reference, checks):
    rows=[]
    for record in reference.to_dict('records'):
        group = formal if str(record['scope']) == 'all' else formal[formal.year.eq(int(record['scope']))]
        observed = group.surface_LODhalf
        predicted = group[f"{record['scenario']}_{record['mapping']}_prediction_ugL"]
        c = concentration_metrics(observed, predicted)
        t = threshold_metrics(observed, predicted, record['threshold_ugL'])
        translated = dict(MAE_log10_offset001=c['MAE_log10'], prediction_max_ugL=c['prediction_max'],
                          prediction_p95_ugL=c['prediction_p95'], prediction_median_ugL=c['prediction_median'], n=c['n'])
        translated.update({k:t[k] for k in ['observed_positive','predicted_positive','TP','FP','FN','TN']})
        _check_equal(list(translated.values()), [record[k] for k in translated],
                     f"Archived statistics: {record['scenario']}, {record['mapping']}, {record['scope']}, {record['threshold_ugL']}", checks)
        rows.append(dict(scenario=record['scenario'], mapping=record['mapping'], scope=record['scope'],
                         threshold_ugL=record['threshold_ugL'], max_absolute_difference=float(max(abs(translated[k]-record[k]) for k in translated))))
    return pd.DataFrame(rows)


def _paired_differences(formal, common):
    """Original 2,000 paired station-year bootstrap draws; no predictions refit."""
    errors = common.apply(lambda predicted: abs(log_errors(formal.surface_LODhalf, predicted)))
    rng = np.random.default_rng(290222)
    blocks = [group.index.to_numpy() for _,group in formal.groupby(['year','station_name'])]
    if len(blocks) != 32:
        raise ValueError('Expected 32 complete station-year blocks.')
    draws = [np.concatenate([blocks[i] for i in rng.integers(0,len(blocks),len(blocks))]) for _ in range(2000)]
    annual = [np.concatenate([formal[formal.year.eq(y)].index.to_numpy() for y in years])
              for years in itertools.product(sorted(formal.year.unique()), repeat=4)]
    rows=[]
    for model in ['Biomass source RevB','Dynamic q80 RevB','Fixed capacity RevB','Model biomass',
                  'Model biomass and composition','Training median']:
        delta = errors[model] - errors['Dynamic RevB']
        bootstrap = np.array([delta.loc[index].mean() for index in draws])
        by_year = np.array([delta.loc[index].mean() for index in annual])
        rows.append(dict(model=READABLE_MODELS[model], source_model=model, n=len(formal),
                         mean_difference=float(delta.mean()), block95_low=float(np.quantile(bootstrap,.025)),
                         block95_high=float(np.quantile(bootstrap,.975)), year95_low=float(np.quantile(by_year,.025)),
                         year95_high=float(np.quantile(by_year,.975))))
    return pd.DataFrame(rows)


def _biomass(data_dir, checks):
    b64 = pd.read_csv(data_dir / 'biomass_surface_2022_64.csv')
    if len(b64) != 64 or b64[['year','station_name','date']].duplicated().any():
        raise ValueError('64-row biomass coverage mismatch.')
    if b64.uses_model_fraction_in_target.any():
        raise ValueError('Biomass proxy must not contain a model-derived fraction.')
    _check_equal(b64.field_chla_ugL/10, b64.proxy_total_algae_gDW_m3, 'Independent Chl-a/10 biomass proxy, 64 rows', checks)
    rows=[]
    for column,label in [('aggregation_offswitch_biomass_gDW_m3','Surface expression off'),
                         ('shallow_retention_max3_biomass_gDW_m3','Depth-retention surface expression')]:
        metrics=concentration_metrics(b64.proxy_total_algae_gDW_m3,b64[column],offset=0)
        rows.append(dict(dataset='2022 surface-expression comparison; 64 common rows',model=label,**metrics))
    b181 = pd.read_csv(data_dir / 'biomass_development_and_holdout_181.csv')
    if len(b181)!=181 or b181[['year','station_name','date']].duplicated().any():
        raise ValueError('181-row biomass coverage mismatch.')
    if not (b181[b181.year.eq(2019)].evaluation_role.eq('independent_validation').all()
            and len(b181[b181.year.eq(2019)])==40):
        raise ValueError('2019 is a frozen 40-row conditional holdout.')
    for role,g in b181.groupby('evaluation_role'):
        for col,label in [('R18_surface_x145','Metabolic capacity, fixed x1.45 surface scale'),
                          ('R18E_surface_x145','Empirical capacity, fixed x1.45 surface scale')]:
            rows.append(dict(dataset=f'181-row set: {role}',model=label,**concentration_metrics(g.obs,g[col],offset=0)))
    changes=[]
    for (year,station), group in b181.groupby(['year','station_name'],sort=False):
        g=group.sort_values('date')
        observed=np.diff(np.log10(g.obs));predicted=np.diff(np.log10(g.R18_surface_x145))
        for i,(o,p) in enumerate(zip(observed,predicted),1):
            changes.append(dict(year=int(year),station_name=station,evaluation_role=g.evaluation_role.iloc[0],
                                date_from=g.date.iloc[i-1],date_to=g.date.iloc[i],
                                delta_log10_obs=float(o),delta_log10_model=float(p),direction_match=bool(np.sign(o)==np.sign(p))))
    if len(changes)!=149:
        raise ValueError('Expected 149 adjacent dates within station-years.')
    development=pd.read_csv(data_dir/'biomass_development_metrics_reference.csv')
    holdout=pd.read_csv(data_dir/'biomass_2019_metrics_reference.csv')
    for table,role,cols in [(development,'development',{'R18_surface':'R18_surface_x145','R18E_surface':'R18E_surface_x145'}),
                            (holdout,'independent_validation',{'R18_surface_x145':'R18_surface_x145','R18E_surface_x145':'R18E_surface_x145'})]:
        g=b181[b181.evaluation_role.eq(role)]
        for source_model,column in cols.items():
            reference=table[table.model.eq(source_model)].iloc[0]
            metrics=concentration_metrics(g.obs,g[column],offset=0)
            _check_equal([metrics['n'],metrics['MAE_log10'],metrics['spearman']],
                         [reference.n,reference.log10_MAE,reference.spearman],f'181 biomass reference: {role}, {source_model}',checks)
    return pd.DataFrame(rows),pd.DataFrame(changes)


def reproduce(data_dir, output_dir, bootstrap=True):
    """Write the released evaluation outputs and return the audit summary.

    `bootstrap=False` skips only paired resampling intervals. It does not change
    predictions or concentration/threshold statistics. Defaults reproduce the
    original paired resampling protocol.
    """
    data_dir=Path(data_dir);output_dir=Path(output_dir)
    if data_dir.resolve()==output_dir.resolve():
        raise ValueError('Output directory must differ from frozen inputs.')
    output_dir.mkdir(parents=True,exist_ok=True)
    assets=verify_assets(data_dir)
    read=lambda filename:pd.read_csv(data_dir/filename)
    formal=read('formal_predictions_164.csv').set_index('sample_id',drop=False)
    features=read('target_and_process_features_164.csv').set_index('sample_id',drop=False)
    common=read('common_predictions_164.csv').set_index('sample_id')
    checks=[]
    if len(formal)!=164 or formal.index.duplicated().any():
        raise ValueError('Expected 164 unique official MC station-date keys.')
    if set(formal.index)!=set(features.index) or set(formal.index)!=set(common.index):
        raise ValueError('Prediction/target common-sample identity failed.')
    features=features.loc[formal.index];common=common.loc[formal.index]
    if not formal[KEYS].equals(features[KEYS]):
        raise ValueError('Official station, date, year or window identity mismatch.')
    _check_equal(formal.surface_LODhalf,features.surface_LODhalf,'Exact surface total-MC observation target',checks,atol=0)
    _check_equal(features.surface_depth_m,np.repeat(.75,len(formal)),'0.75-m surface observation support',checks,atol=0)
    _check_equal(formal.revb_factor,1+2*formal.depth_gate*formal.retention_gate,'Depth-retention factor 1+2*depth_gate*retention_gate',checks)
    for name,variable in [('Dynamic RevB','DYNAMIC_METABOLIC_CAPACITY'),('Fixed capacity RevB','FIXED_CAPACITY_DEVELOPMENT_MATCHED'),
                           ('Biomass source RevB','BIOMASS_PROPORTIONAL_DEVELOPMENT_MATCHED')]:
        stem=name.replace(' ','_')
        _check_equal(formal[[stem+'_I_ugL',stem+'_P_ugL',stem+'_D_ugL']].sum(axis=1),features[variable],
                     name+' water-column three-phase sum',checks,atol=1e-10)
        _check_equal(formal[stem+'_D_ugL']+formal.revb_factor*(formal[stem+'_I_ugL']+formal[stem+'_P_ugL']),
                     formal[stem+'_process_ugL'],name+' particulate-only surface expression',checks,atol=1e-10)
        _check_equal(common[name],formal[stem+'_loyo_prediction_ugL'],name+' formal prediction identity',checks)
    _check_equal(common['Dynamic q80 RevB'],formal['Dynamic RevB_refit_q80_prediction_ugL'],'Dynamic q0.80 prediction identity',checks)
    mapping=_mapping_checks(formal,read('formal_loyo_parameters.csv'),read('mapping_intervals_12.csv'),checks)
    _write(mapping,output_dir/'mapping_fixed_parameter_checks.csv')
    _write(_formal_reference(formal,read('formal_metrics_reference.csv'),checks),output_dir/'archived_180_metrics_checks.csv')
    mc=[]
    for model in common.columns:
        for scope,g in _scopes(formal):
            o=g.surface_LODhalf;p=common.loc[g.index,model]
            for threshold in THRESHOLDS:
                mc.append(dict(source_model=model,model=READABLE_MODELS[model],scope=scope,
                               **concentration_metrics(o,p),**threshold_metrics(o,p,threshold)))
    mc=pd.DataFrame(mc);_write(mc,output_dir/'mc_metrics.csv')
    phase=[]
    _check_equal(formal.dynamic_revb_dissolved_fraction,formal.Dynamic_RevB_D_ugL/formal.Dynamic_RevB_process_ugL,'Surface-expression modeled dissolved fraction',checks)
    _check_equal(formal.observed_dissolved_fraction,features.surface_dissolved/(features.surface_dissolved+features.surface_particulate),'Observed dissolved fraction',checks)
    for scope,g in _scopes(formal):
        phase.append(dict(scope=scope,n=len(g),observed_median=float(g.observed_dissolved_fraction.median()),
                          model_median=float(g.dynamic_revb_dissolved_fraction.median()),
                          MAE_fraction=float(abs(g.dynamic_revb_dissolved_fraction-g.observed_dissolved_fraction).mean()),
                          spearman=spearman(g.observed_dissolved_fraction,g.dynamic_revb_dissolved_fraction)))
    _write(pd.DataFrame(phase),output_dir/'phase_fraction_summary.csv')
    risk=read('risk_fpr10_predictions_164.csv').set_index('sample_id',drop=False).loc[formal.index]
    if not risk[KEYS].equals(formal[KEYS]):
        raise ValueError('S27 risk branch sample identity mismatch.')
    _check_equal(risk.surface_LODhalf,formal.surface_LODhalf,'Risk branch independent target identity',checks,atol=0)
    _check_equal(risk.background_loyo_prediction_ugL,common['Dynamic RevB'],'Risk branch background equals formal OLS',checks)
    _check_equal(risk.q80_loyo_prediction_ugL,common['Dynamic q80 RevB'],'Risk branch q0.80 equals formal q0.80',checks)
    riskrows=[]
    for scope,g in _scopes(risk):
        for threshold in THRESHOLDS:
            riskrows.append(dict(model='Two-stage risk/tail branch (fpr10; frozen nested LOYO selection)',scope=scope,
                                 **concentration_metrics(g.surface_LODhalf,g.two_stage_prediction_ugL),
                                 **threshold_metrics(g.surface_LODhalf,g.two_stage_prediction_ugL,threshold)))
    riskrows=pd.DataFrame(riskrows);_write(riskrows,output_dir/'risk_fpr10_metrics.csv')
    riskreference=read('risk_fpr10_metrics_reference.csv')
    reference_models=riskreference.model.unique().tolist()
    candidates=[name for name in reference_models if 'two_stage' in name]
    if len(candidates)!=1:
        raise ValueError('Cannot uniquely identify archived two-stage model.')
    for row in riskrows.to_dict('records'):
        ref=riskreference[(riskreference.model.eq(candidates[0]))&(riskreference.scope.astype(str).eq(row['scope']))
                          &(riskreference.threshold_ugL.eq(row['threshold_ugL']))].iloc[0]
        fields=['n','observed_positive','predicted_positive','TP','FP','FN','TN']
        _check_equal([row[k] for k in fields],[ref[k] for k in fields],f"S27: {row['scope']}, {row['threshold_ugL']}",checks,atol=0)
    biomass,changes=_biomass(data_dir,checks)
    _write(biomass,output_dir/'biomass_metrics.csv');_write(changes,output_dir/'biomass_adjacent_changes_149.csv')
    groups=read('stratification_variables_124.csv').set_index('sample_id')
    if set(groups.index)!=set(formal[formal.year.isin([2020,2021,2022])].index):
        raise ValueError('S23 must retain its exact 124 rows, separately from 164 MC evaluation rows.')
    _check_equal(groups.observed,formal.loc[groups.index,'surface_LODhalf'],'S23 current target identity',checks,atol=0)
    strata=[]
    for row in read('stratified_metrics_reference.csv').itertuples(index=False):
        if row.dimension=='observed_change':
            index=groups.index[groups.delta_log_chla>0] if row.stratum=='rising' else groups.index[groups.delta_log_chla<0]
        else:
            column='chla_ugL' if row.dimension=='observed_chla' else row.dimension
            v=groups[column];index=v.index[v<=v.median()] if row.stratum=='low' else v.index[v>v.median()]
        e=log_errors(formal.loc[index,'surface_LODhalf'],common.loc[index,row.model])
        _check_equal([len(index),abs(e).mean(),e.mean()],[row.n,row.MAE_log_offset,row.bias_log_offset],
                     f'S23: {row.dimension}, {row.stratum}, {row.model}',checks)
        strata.append(dict(dimension=row.dimension,stratum=row.stratum,source_model=row.model,model=READABLE_MODELS[row.model],
                           n=len(index),MAE_log10=float(abs(e).mean()),bias_log10=float(e.mean())))
    _write(pd.DataFrame(strata),output_dir/'stratified_metrics_124.csv')
    if bootstrap:
        differences=_paired_differences(formal,common)
        reference=read('paired_error_differences_reference.csv').set_index('model')
        for row in differences.itertuples(index=False):
            fields=['mean_difference','block95_low','block95_high','year95_low','year95_high','n']
            _check_equal([getattr(row,k) for k in fields],reference.loc[row.source_model,fields].to_numpy(float),
                         'Paired frozen-prediction intervals: '+row.source_model,checks)
        _write(differences,output_dir/'paired_error_differences.csv')
    # Current S24 intervals are archived coefficients, not newly estimated uncertainty.
    _write(read('mapping_intervals_12.csv'),output_dir/'archived_mapping_intervals_12.csv')
    _write(pd.DataFrame(checks),output_dir/'identity_checks.csv')
    summary=dict(status='PASS_FROZEN_PREDICTION_RESCORING',assets_verified=len(assets),checks_passed=len(checks),
                 model_runs=0,new_fits=0,new_parameter_selection=0,
                 mc_rows=164,biomass_surface_rows=64,biomass_development_rows=141,biomass_2019_holdout_rows=40,
                 stratification_rows=124,mapping_interval_rows=12,
                 protocol=dict(mc_target='Surface particulate+dissolved total MC, ELISA, LOD/2; ug/L',surface_depth_m=.75,
                               log10_mc_offset_ugL=.01,biomass_log10_offset=0,
                               paired_bootstrap_reproduced=bool(bootstrap),paired_station_year_blocks=32,
                               bootstrap_draws=2000 if bootstrap else 0,bootstrap_seed=290222,
                               annual_combinations=256 if bootstrap else 0,
                               year_2019='Original LOYO scoring retained; frozen conditional biomass transfer holdout; no new selection.'),
                 primary_mc=mc[(mc.scope=='all')&(mc.threshold_ugL==5)].to_dict('records'),
                 primary_risk=riskrows[riskrows.scope=='all'].to_dict('records'),
                 biomass=biomass.to_dict('records'),phase_fraction=phase,
                 upstream_model='FACHB905 model and upstream LP generation are outside this release scope.')
    summary=_clean_json(summary)
    (output_dir/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False),encoding='utf-8')
    return summary
