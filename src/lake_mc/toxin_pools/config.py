"""Frozen source-pool defaults used by the formal repaired replay."""
from __future__ import annotations
from dataclasses import dataclass

SOURCE_SCALE = 1.0
EPS = 1e-18
MC_LR_MW = 995.17
BIOMASS_SCALE_GDW_M3_PER_STATE = 15.300047

@dataclass(frozen=True)
class T3SProducerDefaults:
    source_scale: float = SOURCE_SCALE
    toxin_initial_live_quota_mmol_gDW: float = 0.0008
    toxin_growth_yield_mmol_gDW: float = 0.0008
    toxin_quota_target_mmol_gDW: float = 0.0008
    toxin_quota_homeostasis_rate_h: float = 1.0 / 72.0
    toxin_capacity_fraction_theta: float = 0.10
    inactive_clearance_lysis_fraction: float = 0.5

    baseline_pool_content_base_mmol_gDW: float = 0.0008
    baseline_pool_content_maturity_boost_mmol_gDW: float = 0.0010
    baseline_pool_content_lowlight_boost_mmol_gDW: float = 0.0005
    event_pool_content_release_maturity_mmol_gDW: float = 0.0060
    event_pool_content_active_decline_mmol_gDW: float = 0.0020
    event_pool_content_lowlight_release_mmol_gDW: float = 0.0015
    pool_content_max_mmol_gDW: float = 0.0120
    toxic_pool_adjust_tau_h: float = 72.0
    toxic_pool_loading_capacity_multiplier: float = 3.0
    event_pool_loading_capacity_extra_multiplier: float = 3.0
    active_decline_ref: float = 0.45
    active_decline_width: float = 0.25

    intracellular_turnover_rate_h: float = 0.0010
    mortality_release_fraction: float = 0.80
    release_to_dissolved_fraction: float = 0.35
    particulate_release_rate_h: float = 0.0015
    dissolved_degradation_rate_h: float = 0.0008
    event_pool_release_rate_h: float = 0.0060

    spinup_hours: int = 72
    spinup_max_cycles: int = 30
    spinup_convergence_rel_change: float = 1e-6
    negative_zero_tolerance: float = 1e-12

T3S_DEFAULTS = T3SProducerDefaults()
