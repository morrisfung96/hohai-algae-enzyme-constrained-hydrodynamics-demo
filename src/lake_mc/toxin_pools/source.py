"""Formal source-pool update with the repaired dry-weight loss fallback.

Extracted without changing equations from the frozen project implementation.
All arrays and configuration are supplied by the caller; no upstream solver
or project path is accessed. Fluxes use mmol m-3 h-1 and pools mmol m-3.
"""
from __future__ import annotations
import numpy as np
from .config import T3S_DEFAULTS, EPS

def vec(bundle: np.ndarray, fields: dict[str, int], field: str, ti: int, default: float) -> np.ndarray:
    if field not in fields:
        return np.full(bundle.shape[1], float(default), dtype=float)
    out = np.asarray(bundle[ti, :, fields[field]], dtype=float)
    return np.where(np.isfinite(out), out, float(default))

def dynamic_source_vectors(
    bundle: np.ndarray, fields: dict[str, int], ti: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    active = np.clip(vec(bundle, fields, "active_state", ti, 1.0), 0.0, 1.0)
    live = np.maximum(vec(bundle, fields, "B_tox_live_gDW_m3", ti, 0.0), 0.0)
    q_real = vec(bundle, fields, "toxin_capacity_realized_per_gDW_h", ti, np.nan)
    q_upper = vec(bundle, fields, "toxin_capacity_upper_per_gDW_h", ti, 0.0)
    q_use = np.where(np.isfinite(q_real) & (q_real > 0.0), q_real, q_upper)
    candidate = np.maximum(q_use * live * active * T3S_DEFAULTS.source_scale, 0.0)
    cap_default = q_use * np.maximum(live, EPS)
    cap = vec(bundle, fields, "MCLR_capacity_ceiling_flux_mmol_m3_h", ti, np.nan)
    cap = np.maximum(np.where(np.isfinite(cap), cap, cap_default), 0.0)
    synth = np.where(cap > 0.0, np.minimum(candidate, cap), candidate)
    return synth, q_use, live, active, cap

def release_loss_gdw(direct_loss, micro_loss_state_h, ftox, biomass_scale_gDW_m3_per_state):
    scale=float(biomass_scale_gDW_m3_per_state)
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError('A positive, finite upstream state-to-dry-weight scale is required')
    fallback_gdw_m3_h=micro_loss_state_h*ftox*scale
    return np.where(direct_loss>0.0,direct_loss,fallback_gdw_m3_h)

def update_variant(
    bundle: np.ndarray,
    fields: dict[str, int],
    ti: int,
    pools: tuple[np.ndarray, np.ndarray, np.ndarray],
    dt_h: float,
    variant: str,
    constants: dict[str, float],
    *, biomass_scale_gDW_m3_per_state: float,
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], dict[str, np.ndarray]]:
    cfg = T3S_DEFAULTS
    intra, dead, dissolved = pools
    dt_h = max(float(dt_h), 0.0)
    active = np.clip(vec(bundle, fields, "active_state", ti, 1.0), 0.0, 1.0)
    bloom = np.clip(vec(bundle, fields, "bloom_stage_gate_state", ti, 0.0), 0.0, 1.0)
    maturity_raw = vec(bundle, fields, "release_maturity_gate_state", ti, np.nan)
    maturity = np.clip(np.where(np.isfinite(maturity_raw), maturity_raw, bloom), 0.0, 1.0)
    pressure = np.clip(vec(bundle, fields, "release_pressure_state", ti, 0.0), 0.0, 1.0)
    live = np.maximum(vec(bundle, fields, "B_tox_live_gDW_m3", ti, 0.0), 0.0)
    inactive = np.maximum(vec(bundle, fields, "B_tox_inactive_gDW_m3", ti, 0.0), 0.0)
    tox_pool = live + inactive

    if variant == "DYNAMIC_METABOLIC_CAPACITY":
        synthesis, q_unit, _, _, capacity = dynamic_source_vectors(bundle, fields, ti)
    elif variant == "FIXED_CAPACITY_DEVELOPMENT_MATCHED":
        q_unit = np.full_like(live, constants["q_fixed_active_mmol_gDW_h"])
        synthesis = np.maximum(q_unit * live * active * cfg.source_scale, 0.0)
        capacity = np.maximum(q_unit * live * cfg.source_scale, 0.0)
    elif variant == "BIOMASS_PROPORTIONAL_DEVELOPMENT_MATCHED":
        q_unit = np.full_like(live, constants["q_fixed_biomass_mmol_gDW_h"])
        synthesis = np.maximum(q_unit * live * cfg.source_scale, 0.0)
        capacity = synthesis.copy()
    else:
        raise ValueError(variant)

    active_decline = np.clip((cfg.active_decline_ref - active) / max(cfg.active_decline_width, EPS), 0.0, 1.0)
    release_gate = np.clip(
        pressure * (0.45 + 0.55 * maturity) * (0.55 + 0.45 * active_decline) * 0.70,
        0.0,
        1.0,
    )
    baseline_content = np.clip(
        cfg.baseline_pool_content_base_mmol_gDW
        + cfg.baseline_pool_content_maturity_boost_mmol_gDW * maturity,
        0.0,
        cfg.pool_content_max_mmol_gDW,
    )
    event_content = np.clip(
        cfg.event_pool_content_release_maturity_mmol_gDW * pressure * maturity
        + cfg.event_pool_content_active_decline_mmol_gDW * release_gate * active_decline,
        0.0,
        cfg.pool_content_max_mmol_gDW,
    )
    target_content = np.clip(baseline_content + event_content, 0.0, cfg.pool_content_max_mmol_gDW)
    pool_target = np.maximum(tox_pool * target_content, 0.0)
    deficit = np.maximum(pool_target - np.maximum(intra, 0.0), 0.0)
    raw_loading = deficit / max(cfg.toxic_pool_adjust_tau_h, EPS)
    cap_multiplier = np.maximum(
        cfg.toxic_pool_loading_capacity_multiplier
        + cfg.event_pool_loading_capacity_extra_multiplier * release_gate,
        0.0,
    )
    loading_cap = np.maximum(cap_multiplier * q_unit * tox_pool * cfg.source_scale, 0.0)
    standing_loading = np.minimum(raw_loading, loading_cap)
    total_input = synthesis + standing_loading

    direct_loss = (
        np.maximum(vec(bundle, fields, "F_tox_biomass_L_to_I_gDW_m3_h", ti, 0.0), 0.0)
        + np.maximum(vec(bundle, fields, "F_tox_biomass_turnover_gDW_m3_h", ti, 0.0), 0.0)
        + np.maximum(vec(bundle, fields, "F_tox_biomass_clear_gDW_m3_h", ti, 0.0), 0.0)
    )
    micro_loss = np.maximum(vec(bundle, fields, "micro_loss_state_h", ti, 0.0), 0.0)
    ftox_primary = vec(bundle, fields, "f_toxigenic_for_synthesis_state", ti, np.nan)
    ftox_fallback = vec(bundle, fields, "effective_toxic_biomass_fraction_state", ti, 0.0)
    ftox = np.maximum(np.where(np.isfinite(ftox_primary), ftox_primary, ftox_fallback), 0.0)
    direct_loss = release_loss_gdw(direct_loss, micro_loss, ftox, biomass_scale_gDW_m3_per_state)
    quota = np.maximum(intra, 0.0) / np.maximum(tox_pool, EPS)
    mortality_release = cfg.mortality_release_fraction * direct_loss * quota
    turnover_release = cfg.intracellular_turnover_rate_h * np.maximum(intra, 0.0)
    event_release = cfg.event_pool_release_rate_h * np.maximum(intra, 0.0) * release_gate
    release_candidate = np.maximum(mortality_release + turnover_release + event_release, 0.0)
    if dt_h > 0.0:
        available = np.maximum(intra + dt_h * total_input, 0.0) / max(dt_h, EPS)
        release = np.minimum(release_candidate, available)
    else:
        release = np.zeros_like(intra)
    part_to_dissolved = cfg.particulate_release_rate_h * np.maximum(dead, 0.0)
    dissolved_degradation = cfg.dissolved_degradation_rate_h * np.maximum(dissolved, 0.0)
    intra_next = np.maximum(intra + dt_h * (total_input - release), 0.0)
    dead_next = np.maximum(dead + dt_h * ((1.0 - cfg.release_to_dissolved_fraction) * release - part_to_dissolved), 0.0)
    dissolved_next = np.maximum(
        dissolved + dt_h * (cfg.release_to_dissolved_fraction * release + part_to_dissolved - dissolved_degradation),
        0.0,
    )
    expected = intra + dead + dissolved + dt_h * (total_input - dissolved_degradation)
    residual = intra_next + dead_next + dissolved_next - expected
    return (intra_next, dead_next, dissolved_next), {
        "release_flux": release,
        "synthesis_flux": synthesis,
        "standing_loading_flux": standing_loading,
        "total_input_flux": total_input,
        "capacity_flux": capacity,
        "pool_mass_balance_residual": residual,
    }
