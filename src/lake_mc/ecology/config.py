"""Exact source extraction from V16_3_pipeline/v16_3_algal_state_model.py; algorithm retained.
No FACHB905/SBML/LP imports. See ecology extraction provenance.
"""
from __future__ import annotations
from dataclasses import dataclass, asdict
from pathlib import Path
import json
from typing import Optional



MC_LR_MW_G_MOL = 995.17

@dataclass
class V16Config:
    # Identity and response mode.
    biomass_scale_gDW_m3_per_state: float = 1.0
    mc_lr_mw_g_mol: float = MC_LR_MW_G_MOL
    exact_night_light_threshold_Wm2: float = 1e-9
    response_k_env: int = 4
    require_native_v16_response: bool = False  # set True for final manuscript run

    # Nutrient/environment state indices. These do not create qN/qP state.
    K_N_mgN_L: float = 0.35
    K_P_mgP_L: float = 0.030
    light_stress_half_Wm2: float = 80.0
    temp_opt_C: float = 26.0
    temp_width_C: float = 7.5

    # Biomass ecology, intentionally simple and non-MC driven.
    # M4_current keeps legacy V16_3 behavior.  M0-M3 are structural ablation modes
    # for biomass closure only; they do not alter toxin/source equations directly.
    biomass_closure_mode: str = "M4_current"  # M0-M4, plus L0-L4 additive-loss modes
    carrying_capacity_state: float = 0.30
    alpha_other_on_micro: float = 0.80
    alpha_micro_on_other: float = 0.70
    micro_base_loss_rate_h: float = 0.00035
    flush_loss_max_rate_h: float = 0.0
    flush_current_half_mps: float = 0.04
    succession_loss_max_rate_h: float = 0.0
    succession_tau_h: float = 168.0
    succession_advantage_ref_h: float = 0.0
    succession_advantage_width_h: float = 0.003
    succession_pressure_mode: str = "potential_h"  # potential_h | relative_suitability
    succession_suitability_metric: str = "raw_delta"  # raw_delta | relative_advantage | log_advantage
    succession_trigger_function: str = "linear_clip"  # linear_clip | positive_saturation
    succession_suitability_delta_ref: float = 0.05
    succession_suitability_width: float = 0.20
    succession_use_competitor_fraction: bool = False
    succession_effect_mode: str = "loss"  # loss | composition_turnover
    composition_turnover_rate_h: float = 0.0
    micro_suitability_ref_h: float = 0.029147432795074

    # D3 biomass-closure experiment: live -> inactive -> lysis/settling ledger.
    # Defaults are off so historical V16_3/C2 runs are unchanged.  The inactive
    # pool is biomass-only here; toxin release is not coupled during B-stage
    # calibration.
    enable_late_inactivation_pool: bool = False
    late_inactivation_signal_mode: str = "metabolic_decline_saturation"  # metabolic_decline_saturation | off
    late_mu_memory_tau_h: float = 72.0
    late_inactivation_memory_tau_h: float = 96.0
    late_decline_half_per_day: float = 0.025
    late_decline_light_min_Wm2: float = 20.0
    late_saturation_power: float = 1.0
    late_saturation_ref: float = 0.85
    late_saturation_width: float = 0.05
    late_exhaustion_rcrit_per_day: float = 0.015
    late_exhaustion_width_per_day: float = 0.010
    late_burden_threshold_h: float = 72.0
    late_burden_half_h: float = 48.0
    late_burden_decay_tau_h: float = 96.0
    late_damage_accum_rate_h: float = 1.0 / 96.0
    late_damage_repair_tau_h: float = 336.0
    late_damage_half_saturation: float = 0.25
    micro_inactivation_max_rate_h: float = 0.0
    other_inactivation_max_rate_h: float = 0.0
    micro_inactive_lysis_rate_h: float = 0.00125
    other_inactive_lysis_rate_h: float = 0.00125
    micro_inactive_settling_rate_h: float = 0.00125
    other_inactive_settling_rate_h: float = 0.00125
    inactive_observation_equiv_fraction: float = 0.5
    inactive_density_contribution_chi: float = 0.0

    other_mu_max_h: float = 0.014
    other_loss_rate_h: float = 0.00035
    other_light_half_Wm2: float = 80.0

    # Evolving algal state.
    active_min_fraction: float = 0.02
    active_tau_h: float = 96.0
    active_mu_ref_h: float = 0.020
    active_release_pressure_weight: float = 0.45
    active_low_light_weight: float = 0.25
    alloc_tau_h: float = 96.0
    low_light_memory_tau_h: float = 48.0
    release_pressure_tau_h: float = 96.0
    release_pressure_base_loss_rate_h: float = 0.00010
    release_pressure_max_extra_loss_rate_h: float = 0.0012

    # State target weights. These are global, low-dimensional, and use non-MC drivers.
    alloc_stress_weight: float = 0.35
    alloc_maturity_weight: float = 0.35
    alloc_competition_weight: float = 0.15
    alloc_lowlight_weight: float = 0.15

    # Mapping allocation state to ecModel growth-commitment rho.
    # High allocation/stress/maturity -> lower growth commitment, leaving capacity for non-growth programs.
    rho_growth_high: float = 0.80
    rho_growth_low: float = 0.25

    # MC capacity realization. Default 1.0 means use the state-conditioned upper envelope.
    # This is deliberately kept separate and audited. Do not tune per window.
    capacity_use_fraction: float = 1.0

    # f_toxigenic dynamics. Initialization can come from mcyE/total-cyano profile.
    # Effective toxigenic biomass fraction used in source calculations.
    # Raw mcyE/total-cyano qPCR ratios are kept separately as a gene proxy;
    # they should not be used literally as biomass fractions.
    ftox_default: float = 0.020
    ftox_min: float = 1e-5
    ftox_max: float = 0.20
    ftox_tau_h: float = 240.0
    ftox_growth_competition_sensitivity: float = 0.10

    # V16_2 toxic biomass closure switches.
    # mcyE/total-cyano is diagnostic by default. It is not required and should not
    # hard-suppress source strength in high-biomass/high-release-pressure regions.
    v16_variant: str = "V16_3_full2022_initial_fidelity"
    mcye_source_mode: str = "diagnostic"  # off | diagnostic | soft_constraint
    toxic_producer_prior_fraction: float = 0.060
    toxic_producer_prior_min: float = 0.020
    toxic_producer_prior_max: float = 0.120
    toxic_prior_maturity_boost_fraction: float = 0.015
    toxic_prior_release_boost_fraction: float = 0.030
    bloom_stage_gate_min: float = 0.75
    bloom_stage_gate_max: float = 1.15
    release_maturity_gate_min: float = 0.85
    release_maturity_gate_max: float = 1.20
    mcye_soft_min_factor: float = 0.70
    mcye_soft_max_factor: float = 1.50
    mcye_soft_ref_effective_fraction: float = 0.020

    toxic_pool_fraction_max: float = 0.160
    toxic_pool_adjust_tau_h: float = 72.0
    toxic_pool_release_active_suppression: bool = False

    # V16_2C: experimental release-partition pool closure.
    # This is not qMC demand. It defines a bounded intracellular/particulate MC pool
    # target from non-MC states, then limits loading by native ecModel capacity.
    # A separate event pool target activates during high release/maturity/active-decline
    # conditions and is partitioned using the experiment-measured 35% dissolved / 65% retained rule.
    toxic_pool_loading_mode: str = "experimental_release_partition"
    baseline_pool_content_base_mmol_gDW: float = 0.0008
    baseline_pool_content_maturity_boost_mmol_gDW: float = 0.0010
    baseline_pool_content_lowlight_boost_mmol_gDW: float = 0.0005
    event_pool_content_release_maturity_mmol_gDW: float = 0.0060
    event_pool_content_active_decline_mmol_gDW: float = 0.0020
    event_pool_content_lowlight_release_mmol_gDW: float = 0.0015
    pool_content_max_mmol_gDW: float = 0.0120
    toxic_pool_loading_capacity_multiplier: float = 3.0
    event_pool_loading_capacity_extra_multiplier: float = 3.0
    active_decline_ref: float = 0.45
    active_decline_width: float = 0.25

    # V16_3: high-biomass / dominance / scum memory closure.
    # This is the cleanest V14-derived mechanism to restore first. It does not
    # introduce qN/qP or qMC/LR demand. The memory states only affect toxic-pool
    # eligibility, event-pool target, release/scum expression, and are fully audited.
    enable_patch_memory: bool = True
    high_biomass_memory_tau_h: float = 72.0
    dominance_memory_tau_h: float = 72.0
    bloom_patch_memory_tau_h: float = 96.0
    scum_memory_tau_h: float = 48.0
    patch_B_micro_ref_state: float = 0.018
    patch_B_micro_width_state: float = 0.008
    patch_dominance_ref: float = 0.12
    patch_dominance_width: float = 0.06
    patch_memory_toxic_prior_boost_fraction: float = 0.020
    patch_memory_baseline_pool_boost_mmol_gDW: float = 0.0008
    patch_memory_event_pool_boost_mmol_gDW: float = 0.0035
    patch_memory_release_gate_boost: float = 0.35
    patch_memory_surface_extra: float = 3.0

    # MC pool/release dynamics. Synthesis/loading, release, and partition are separated.
    intracellular_turnover_rate_h: float = 0.0010
    mortality_release_fraction: float = 0.80
    release_to_dissolved_fraction: float = 0.35  # user experimental constraint: 35% dissolved
    retained_particulate_fraction: float = 0.65  # user experimental constraint: 65% retained particulate/cell-bound
    event_pool_release_rate_h: float = 0.0060
    particulate_release_rate_h: float = 0.0015
    dissolved_degradation_rate_h: float = 0.0008

    # Toxin-biomass interface audit. Defaults are off/diagnostic-safe so
    # existing V16.3 source runs are unchanged unless explicitly enabled.
    enable_toxin_interface_ledger: bool = False
    toxin_production_mode: str = "off"  # off | capacity_ceiling | growth_associated | quota_homeostasis | capacity_fraction
    toxin_initial_live_quota_mmol_gDW: float = 0.0008
    toxin_growth_yield_mmol_gDW: float = 0.0008
    toxin_quota_target_mmol_gDW: float = 0.0008
    toxin_quota_homeostasis_rate_h: float = 1.0 / 72.0
    toxin_capacity_fraction_theta: float = 0.10
    clearance_toxin_fate_mode: str = "unresolved"  # unresolved | lysis | settling | partitioned

    # Aggregation/surface-expression state: observation layer only.
    # It converts bulk/surface pools into a station surface-sample proxy and does
    # not increase transported lake mass.
    enable_aggregation_state: bool = True
    aggregation_tau_h: float = 48.0
    aggregation_surface_amplification: float = 0.25
    aggregation_wind_half_mps: float = 4.0
    surface_sample_enrichment_max_extra: float = 8.0
    surface_sample_factor_max: float = 10.0
    surface_sample_release_weight: float = 0.60
    surface_sample_maturity_weight: float = 0.40

    # Runtime and speed.
    local_timezone: str = "America/Detroit"
    max_grids: int | None = None
    # For small smoke runs, do not simply take the first N grid IDs.
    # station_coverage keeps grids nearest to 2022 mcyE stations, then fills the rest with spatially spread grids.
    max_grids_strategy: str = "station_coverage"  # station_coverage | spatial | first
    n_workers: int = 1
    parallel_chunk_grids: int = 120
    progress_every_grids: int = 1000
    keep_parallel_parts: bool = False
    fast_response_audit_from_chunks: bool = True
    write_daily_summary: bool = True

    @classmethod
    def from_json(cls, path: Path | str) -> "V16Config":
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        # Backward compatibility for V16 prototypes that used the stronger word "lysis".
        # V16 now uses "release_pressure" because this state is a pressure/proxy,
        # not a direct fraction of lysed cells.
        legacy_map = {
            "lysis_tau_h": "release_pressure_tau_h",
            "lysis_base_rate_h": "release_pressure_base_loss_rate_h",
            "lysis_max_extra_rate_h": "release_pressure_max_extra_loss_rate_h",
        }
        for old_key, new_key in legacy_map.items():
            if old_key in data and new_key not in data:
                data[new_key] = data[old_key]
        allowed = set(cls.__dataclass_fields__.keys())
        data = {k: v for k, v in data.items() if k in allowed}
        return cls(**data)

    def to_json(self, path: Path | str) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2, ensure_ascii=False), encoding="utf-8")
