"""Replay the formal three-pool source/fate chain from frozen ecological states.

The input directory contains a full 5000-grid state bundle, its identity and
constants, corrected initial pools, and the exported 4723-grid wet geometry.
This module performs no ecological rerun, parameter fitting, or upstream LP.
"""
from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from .config import BIOMASS_SCALE_GDW_M3_PER_STATE, MC_LR_MW, T3S_DEFAULTS
from .inventory import InventoryHandoff
from .source import update_variant
from .transport import transport_loss_step


VARIANTS = {
    'dynamic': 'DYNAMIC_METABOLIC_CAPACITY',
    'fixed': 'FIXED_CAPACITY_DEVELOPMENT_MATCHED',
    'proportional': 'BIOMASS_PROPORTIONAL_DEVELOPMENT_MATCHED',
}


class ReplayInputError(ValueError):
    """A frozen replay asset does not satisfy its declared identity."""


def load_inputs(directory: str | Path, variant: str) -> dict:
    """Read only explicitly supplied portable assets, never source-project paths."""
    directory = Path(directory)
    if variant not in VARIANTS:
        raise ReplayInputError(f'Unknown variant {variant!r}; choose {list(VARIANTS)}')
    metadata = json.loads((directory / 'source_identity.json').read_text('utf-8-sig'))
    fate = json.loads((directory / 'fate_config.json').read_text('utf-8-sig'))
    constants = json.loads((directory / 'source_constants.json').read_text('utf-8-sig'))
    if 'constants' in constants:
        constants = constants['constants']
    bundle = np.load(directory / 'source_states.npy', mmap_mode='r', allow_pickle=False)
    names = metadata.get('field_names', metadata.get('fields'))
    fields = names if isinstance(names, dict) else {str(n): i for i, n in enumerate(names)}
    times = pd.DatetimeIndex(pd.to_datetime(metadata['times'], utc=True)).to_numpy()
    source_ids = np.asarray(metadata['grid_ids']).astype(str)
    if bundle.ndim != 3 or bundle.shape[:2] != (len(times), len(source_ids)):
        raise ReplayInputError('State shape, timestamps, and source grid IDs differ')
    if 'shape' in metadata and list(bundle.shape) != list(metadata['shape']):
        raise ReplayInputError('State array differs from the declared shape')
    if metadata.get('status', 'COMPLETE') != 'COMPLETE':
        raise ReplayInputError('Source bundle is not marked complete')
    if fate['source_pool_defaults'] != asdict(T3S_DEFAULTS):
        raise ReplayInputError('Source defaults differ from the frozen formal defaults')
    if fate['MC_LR_MW'] != MC_LR_MW or fate['biomass_scale_gDW_m3_per_state'] != BIOMASS_SCALE_GDW_M3_PER_STATE:
        raise ReplayInputError('Molecular weight or dry-weight conversion identity differs')
    with np.load(directory / 'frozen_geometry.npz', allow_pickle=False) as z:
        take = z['source_indices'].copy()
        ids = z['grid_ids'].astype(str).tolist()
        if not np.array_equal(source_ids[take], np.asarray(ids)):
            raise ReplayInputError('Wet receiver IDs do not match selected source columns')
        if not np.array_equal(pd.to_datetime(z['times_utc'], utc=True), pd.to_datetime(times, utc=True)):
            raise ReplayInputError('Velocity timestamps differ from ecological timestamps')
        geom = {'grid_ids': ids, 'volume': z['receiver_volume_m3'].copy(),
                'depth': z['depth'].copy(), 'face_width': z['face_width'].copy(),
                'edges': {n: z['edge_' + n].copy() for n in ['src', 'dst', 'ex', 'ey', 'dist_m']}}
        for name in ['u_particle', 'v_particle', 'u_dissolved', 'v_dissolved']:
            geom[name] = z[name].copy()
            if geom[name].shape != (len(times), len(ids)):
                raise ReplayInputError(f'Invalid velocity shape: {name}')
        contract = InventoryHandoff(tuple(ids), tuple(ids), z['source_volume_m3'], z['receiver_volume_m3'])
        W = csr_matrix((z['W_data'], z['W_indices'], z['W_indptr']), shape=tuple(z['W_shape']))
        stations = z['station_names'].astype(str).tolist()
    if W.shape != (len(stations), len(ids)) or not np.allclose(np.asarray(W.sum(axis=1)).ravel(), 1., atol=1e-14, rtol=0.):
        raise ReplayInputError('Station weights have invalid dimensions or normalization')
    with np.load(directory / f'initial_{variant}.npz', allow_pickle=False) as z:
        pools = tuple(z[n].copy() for n in ['intra', 'dead', 'dissolved'])
    if any(p.shape != (len(source_ids),) or not np.isfinite(p).all() or (p < 0).any() for p in pools):
        raise ReplayInputError('Corrected initial pools must cover all source columns')
    return {'bundle': bundle, 'fields': fields, 'times': times, 'geom': geom,
            'take': take, 'contract': contract, 'W': W, 'stations': stations,
            'initial_pools': pools, 'constants': constants, 'fate': fate,
            'metadata': metadata, 'variant': VARIANTS[variant]}


def replay(data: dict, *, max_hours: int | None = None) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    """Run the unchanged formal stepping order from the retained initial pools.

    ``max_hours`` verifies an initial slice and is explicitly labelled partial.
    It never uses a later timestamp with the first-timestamp initial pools.
    """
    bundle, fields, times = data['bundle'], data['fields'], data['times']
    geom, take, contract = data['geom'], data['take'], data['contract']
    W, stations = data['W'], data['stations']
    constants, variant, fate = data['constants'], data['variant'], data['fate']
    cfg = SimpleNamespace(**fate['transport_config'])
    pools = tuple(z.copy() for z in data['initial_pools'])
    V = geom['volume']
    Vs = contract.source_volume_m3
    factor = Vs / V
    dead, diss = contract.to_inventory(pools[1][take]), contract.to_inventory(pools[2][take])
    n = len(times) if max_hours is None else min(int(max_hours), len(times))
    if n <= 0:
        raise ReplayInputError('At least one hour must be requested')
    m = len(stations)
    names = ['I_ugL', 'P_ugL', 'D_ugL', 'C_ugL', 'J_mmol_m3_h',
             'loading_mmol_m3_h', 'release_mmol_m3_h', 'I_source_ugL',
             'J_receiver_mmol_m3_h', 'release_source_mmol_m3_h']
    values = np.zeros((n, m, len(names)))
    ledger = []
    pool_error = mass_error = inventory_error = 0.
    for ti in range(n):
        timestamp = times[ti]
        local_i = pools[0][take]
        nxt, diag = update_variant(
            bundle, fields, ti, pools, 1., variant, constants,
            biomass_scale_gDW_m3_per_state=fate['biomass_scale_gDW_m3_per_state'],
        )
        ci = contract.to_inventory(local_i) / V
        pool_error = max(pool_error, float(np.max(abs(diag['pool_mass_balance_residual'][take]))))
        inventory_error = max(inventory_error, float(np.max(abs(ci * V - local_i * Vs) / np.maximum(local_i * Vs, 1e-30))))
        J, load, rel = [diag[name][take] for name in ('synthesis_flux', 'standing_loading_flux', 'release_flux')]
        fields_out = [ci * MC_LR_MW, dead / V * MC_LR_MW, diss / V * MC_LR_MW,
                      (ci + (dead + diss) / V) * MC_LR_MW, J, load, rel * factor,
                      local_i * MC_LR_MW, J * factor, rel]
        for j, value in enumerate(fields_out):
            values[ti, :, j] = W @ value
        # The original full replay omits a mobile update after its final hour.
        # An initial partial slice can also omit that unobserved final update.
        if ti + 1 < n:
            if (pd.Timestamp(times[ti + 1]) - pd.Timestamp(timestamp)).total_seconds() != 3600:
                raise ReplayInputError('Frozen replay requires consecutive hourly states')
            before = float((dead + diss).sum())
            released = contract.to_inventory(rel)
            transfer = np.minimum(np.maximum(T3S_DEFAULTS.particulate_release_rate_h * dead, 0), dead)
            db = np.maximum(dead + (1 - T3S_DEFAULTS.release_to_dissolved_fraction) * released - transfer, 0)
            ds = np.maximum(diss + T3S_DEFAULTS.release_to_dissolved_fraction * released + transfer, 0)
            dead, diss, td = transport_loss_step(db, ds, ti, 1., geom, cfg)
            settling = float(td['particulate_settling_loss'].sum())
            decay = float(td['dissolved_degradation_loss'].sum())
            after = float((dead + diss).sum())
            residual = (after - before - released.sum() + settling + decay) / max(before, after, 1.)
            mass_error = max(mass_error, abs(residual))
            ledger.append(dict(time=str(timestamp), source_release_mmol=float(rel @ Vs),
                receiver_release_mmol=float(released.sum()), mobile_before_mmol=before,
                mobile_after_mmol=after, settling_mmol=settling, degradation_mmol=decay,
                relative_residual=residual, source_synthesis_mmol=float(J @ Vs),
                source_loading_mmol=float(load @ Vs), source_intra_mmol=float(local_i @ Vs),
                receiver_intra_mmol=float(ci @ V)))
            if not (np.isfinite(dead).all() and np.isfinite(diss).all() and (dead >= 0).all() and (diss >= 0).all()):
                raise ArithmeticError('Mobile inventories lost finiteness or nonnegativity')
        pools = nxt
    used_times = times[:n]
    hour = pd.DataFrame(values.reshape(n * m, len(names)), columns=names)
    hour.insert(0, 'station_name', np.tile(stations, n))
    hour.insert(0, 'time', np.repeat(pd.to_datetime(used_times, utc=True).astype(str), m))
    hour['date'] = pd.to_datetime(hour.time, utc=True).dt.tz_convert('America/Detroit').dt.strftime('%Y-%m-%d')
    hour['elapsed_day'] = (pd.to_datetime(hour.time, utc=True) - pd.to_datetime(times[0], utc=True)).dt.total_seconds() / 86400
    ag = {name: (name, 'mean') for name in names}
    ag.update(n_hours=('time', 'size'), elapsed_day=('elapsed_day', 'mean'))
    daily = hour.groupby(['station_name', 'date'], as_index=False).agg(**ag)
    year = int(fate['year'])
    daily['year'], daily['variant'], daily['scenario'] = year, variant, 'retained'
    audit = dict(year=year, variant=variant, scenario='retained', pool_residual_max=pool_error,
                 mobile_relative_residual_max=mass_error, inventory_relative_error_max=inventory_error,
                 n_station_hours=len(hour), n_station_days=len(daily), wet_grids=len(V),
                 source_hours_available=len(times), hours_replayed=n, full_window=(n == len(times)),
                 scope='Frozen ecological-state source/fate/transport replay; no ecology recomputation',
                 parameters_fitted=False, ecology_rerun=False, historical_native_flux_validated=False)
    if pool_error >= 1e-9 or mass_error >= 1e-12 or inventory_error >= 1e-12:
        raise ArithmeticError(f'Formal conservation tolerances exceeded: {audit}')
    return hour, daily, pd.DataFrame(ledger), audit


def run_replay(directory: str | Path, output: str | Path, variant: str, *, max_hours: int | None = None) -> dict:
    """Write station outputs and a conservation audit from the supplied asset pack."""
    data = load_inputs(directory, variant)
    hour, daily, ledger, audit = replay(data, max_hours=max_hours)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    hour.to_csv(output / 'station_hourly.csv.gz', index=False, encoding='utf-8-sig')
    daily.to_csv(output / 'station_daily.csv', index=False, encoding='utf-8-sig')
    ledger.to_csv(output / 'mobile_mass_ledger.csv.gz', index=False, encoding='utf-8-sig')
    (output / 'audit.json').write_text(json.dumps(audit, indent=2), encoding='utf-8')
    return audit
