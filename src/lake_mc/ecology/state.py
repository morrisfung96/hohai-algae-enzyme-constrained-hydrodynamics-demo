"""Exact source extraction from scripts/engine_R18_stateful_capacity_candidate.py; algorithm retained.
No FACHB905/SBML/LP imports. See ecology extraction provenance.
"""
from __future__ import annotations
from typing import Any
import numpy as np
import pandas as pd
from .production import project_two_species_conservative_production

SCENARIOS = {'OE1_envpow0p75_existing': {'allocation_mode': 'conservative_projected_phi_allocation', 'phi_M_mode': 'current_times_active', 'phi_O_mode': 'bottleneck', 'gamma_inventory': 1.0, 'other_capacity_mode': 'env_power_existing_biomass', 'other_capacity_env_power': 0.75, 'other_capacity_recruitable_fraction': 0.0, 'status': 'valid', 'complexity_rank': 1, 'description': 'Partial release of other hard-capacity environmental double gating: env_single^0.75.', 'dynamic_capacity_enabled': True, 'dynamic_capacity_floor_fraction': 0.5, 'dynamic_capacity_tau_h': 72.0, 'dynamic_excess_tau_h': 72.0}}

EPS = 1e-30

def safe_div(a: float, b: float, default: float = np.nan) -> float:
    if not np.isfinite(a) or not np.isfinite(b) or abs(float(b)) < EPS:
        return float(default)
    return float(a) / float(b)

def sat(x: float, k: float) -> float:
    xx = max(float(x), 0.0)
    kk = max(float(k), 1e-12)
    return float(xx / (xx + kk))

def sigmoid(x: float, midpoint: float = 0.5, width: float = 0.1) -> float:
    ww = max(float(width), 1e-12)
    z = np.clip((float(x) - float(midpoint)) / ww, -60.0, 60.0)
    return float(1.0 / (1.0 + np.exp(-z)))

def temp_suitability(temp_C: float, opt_C: float, width_C: float) -> float:
    ww = max(float(width_C), 1e-6)
    return float(np.exp(-0.5 * ((float(temp_C) - float(opt_C)) / ww) ** 2))

def phase_from_date(date: str) -> str:
    d = pd.Timestamp(date)
    if d <= pd.Timestamp("2022-07-11"):
        return "initial_to_20220711"
    if d <= pd.Timestamp("2022-07-18"):
        return "fraction_lock_20220712_20220718"
    if d <= pd.Timestamp("2022-07-25"):
        return "post_lock_to_20220725"
    if d <= pd.Timestamp("2022-08-02"):
        return "peak_approach_20220726_20220802"
    if d <= pd.Timestamp("2022-08-09"):
        return "peak_hold_20220803_20220809"
    if d <= pd.Timestamp("2022-08-17"):
        return "early_decline_20220810_20220817"
    return "terminal_late_20220818_20220823"

def basin_from_lon(lon: float) -> str:
    if not np.isfinite(lon):
        return "unknown_basin"
    if lon <= -82.5:
        return "western_basin_lon_le_-82.5"
    if lon <= -80.5:
        return "central_basin_-82.5_to_-80.5"
    return "eastern_basin_lon_gt_-80.5"

ACCUM_COLS = [
    "n_grid_hours",
    "n_grid_days",
    "sum_BM",
    "sum_BO",
    "sum_BM_inactive",
    "sum_BO_inactive",
    "sum_BL",
    "sum_BI",
    "sum_BT_intact",
    "first_BM",
    "first_BO",
    "first_BM_inactive",
    "first_BO_inactive",
    "first_BL",
    "first_BI",
    "first_BT_intact",
    "last_BM",
    "last_BO",
    "last_BM_inactive",
    "last_BO_inactive",
    "last_BL",
    "last_BI",
    "last_BT_intact",
    "sum_micro_mu_max",
    "sum_micro_mu_committed",
    "sum_micro_density_gate",
    "sum_micro_active_gate",
    "sum_micro_rho_gate",
    "sum_micro_realized",
    "sum_other_mu_max",
    "sum_other_temp_gate",
    "sum_other_light_gate",
    "sum_other_N_gate",
    "sum_other_P_gate",
    "sum_other_nutrient_gate",
    "sum_other_env_product",
    "sum_other_env_single_limiter",
    "sum_other_density_gate",
    "sum_other_realized",
    "sum_BT_phi_micro",
    "sum_BT_phi_other",
    "sum_BT_phi_ratio_micro_over_other",
    "sum_BT_ln_phi_micro_over_other",
    "sum_BT_phi_micro_lt_other",
    "sum_BT_phi_ratio_lt_0p5",
    "sum_BT_phi_ratio_gt_2",
    "sum_BT_w_micro_allocation",
    "sum_BT_w_micro_phi_pure",
    "sum_BT_w_minus_p_micro",
    "sum_BT_logit_w_minus_logit_p",
    "sum_BT_L_B_other_over_micro",
    "sum_BT_L_phi_other_over_micro",
    "sum_BT_L_w_other_over_micro_gamma",
    "sum_BT_logit_w_micro_gamma",
    "sum_BT_succession_signal",
    "sum_BT_succession_pressure",
    "sum_BT_suitability_relative_advantage",
    "sum_GT_ref_flux",
    "sum_micro_growth_ref_flux",
    "sum_other_growth_ref_flux",
    "sum_micro_growth_allocation_delta_flux",
    "sum_other_growth_allocation_delta_flux",
    "sum_allocation_residual_flux",
    "sum_micro_capacity_flux",
    "sum_other_capacity_flux",
    "sum_other_capacity_basis",
    "sum_other_capacity_basis_env_multiplier",
    "sum_other_recruitable_capacity_basis",
    "sum_micro_capacity_slack_flux",
    "sum_other_capacity_slack_flux",
    "sum_projection_distance_flux",
    "sum_unallocatable_budget_flux",
    "sum_target_infeasible_flux",
    "sum_capacity_budget_inconsistent_flux",
    "sum_micro_capacity_binding_hours",
    "sum_other_capacity_binding_hours",
    "sum_target_feasible_hours",
    "sum_micro_growth_flux",
    "sum_other_growth_flux",
    "sum_micro_loss_flux",
    "sum_other_loss_flux",
    "sum_composition_turnover_flux",
    "sum_micro_inactivation_flux",
    "sum_other_inactivation_flux",
    "sum_micro_inactive_lysis_flux",
    "sum_other_inactive_lysis_flux",
    "sum_micro_inactive_settling_flux",
    "sum_other_inactive_settling_flux",
    "sum_inactive_module_residual_flux",
    "sum_BM_late_memory",
    "sum_BM_late_damage",
    "sum_BM_late_raw_pressure",
    "sum_BM_late_decline_signal",
    "sum_BM_late_saturation_gate",
    "sum_BM_late_exhaustion_gate",
    "sum_BM_micro_inactivation_rate_h",
    "sum_BM_surface_factor",
    "sum_BM_surface_equiv",
    "sum_BT_surface_equiv",
    "sum_negative_correction_mass",
    "sum_state_balance_residual",
]

def zero_accum() -> dict[str, float]:
    return {c: 0.0 for c in ACCUM_COLS}

def state_from_initial(init_row: dict[str, Any], cfg: Any) -> dict[str, float]:
    raw_ftox_gene = init_row.get("f_toxigenic_gene_proxy_raw_state", np.nan)
    try:
        ftox_gene = float(raw_ftox_gene)
    except Exception:
        ftox_gene = np.nan
    return {
        "Bm": max(float(init_row.get("B_micro_state", 0.001) or 0.001), 0.0),
        "Bo": max(float(init_row.get("B_other_state", 0.001) or 0.001), 0.0),
        "Bm_inactive": 0.0,
        "Bo_inactive": 0.0,
        "ftox": float(np.clip(init_row.get("f_toxigenic_state", cfg.ftox_default), cfg.ftox_min, cfg.ftox_max)),
        "ftox_gene_proxy_raw": ftox_gene if np.isfinite(ftox_gene) else np.nan,
        "active_state": float(np.clip(init_row.get("active_state", 0.8), cfg.active_min_fraction, 1.0)),
        "allocation_state": float(np.clip(init_row.get("allocation_state", 0.2), 0.0, 1.0)),
        "release_pressure_state": float(np.clip(init_row.get("release_pressure_state", 0.0), 0.0, 1.0)),
        "low_light_memory": 0.0,
        "aggregation_state": 0.0,
        "high_biomass_memory_state": 0.0,
        "dominance_memory_state": 0.0,
        "bloom_patch_memory_state": 0.0,
        "scum_accumulation_memory_state": 0.0,
        "succession_pressure_state": 0.0,
        "late_inactivation_memory_state": 0.0,
        "late_viability_burden_h": 0.0,
        "late_damage_burden_state": 0.0,
        "mu_capacity_memory_h": np.nan,
        "dynamic_capacity_memory_fraction": 1.0,
    }

def new_daily_accum() -> dict[str, Any]:
    rec = zero_accum()
    rec["_started"] = False
    return rec

def update_daily_accum(rec: dict[str, Any], *, Bm: float, Bo: float,
                       Bm_inactive: float, Bo_inactive: float,
                       Bm_next: float, Bo_next: float,
                       Bm_inactive_next: float, Bo_inactive_next: float,
                       mu_max: float, mu_comm: float, density_micro: float, active: float,
                       rho_gate: float, mu_micro: float, other_mu_max: float, temp_gate: float,
                       light_gate: float, n_gate: float, p_gate: float, nutrient_gate: float,
                       env_product: float, env_single: float, density_other: float, mu_other: float,
                       phi_micro: float, phi_other: float, phi_ratio_micro_over_other: float,
                       ln_phi_micro_over_other: float,
                       w_micro_allocation: float, w_micro_phi_pure: float,
                       w_minus_p_micro: float, logit_w_minus_logit_p: float,
                       L_B_other_over_micro: float,
                       L_phi_other_over_micro: float,
                       L_w_other_over_micro_gamma: float,
                       logit_w_micro_gamma: float,
                       succession_signal: float, succession_pressure: float,
                       suitability_relative_advantage: float,
                       GT_ref: float, micro_growth_ref: float,
                       other_growth_ref: float, allocation_residual: float,
                       micro_capacity_flux: float, other_capacity_flux: float,
                       other_capacity_basis: float, other_capacity_env_multiplier: float,
                       other_recruitable_capacity_basis: float,
                       micro_capacity_slack_flux: float, other_capacity_slack_flux: float,
                       projection_distance_flux: float, unallocatable_budget_flux: float,
                       target_feasible: bool, capacity_budget_inconsistent: bool,
                       micro_capacity_binding: bool, other_capacity_binding: bool,
                       micro_growth: float, other_growth: float, micro_loss: float,
                       other_loss: float, composition_turnover_flux: float,
                       micro_inactivation_flux: float, other_inactivation_flux: float,
                       micro_inactive_lysis_flux: float, other_inactive_lysis_flux: float,
                       micro_inactive_settling_flux: float, other_inactive_settling_flux: float,
                       inactive_module_residual_flux: float,
                       late_memory: float, late_damage: float, late_raw_pressure: float,
                       late_decline_signal: float, late_saturation_gate: float,
                       late_exhaustion_gate: float, micro_inactivation_rate_h: float,
                       surface_factor: float, neg_corr: float, residual: float) -> None:
    bl = Bm + Bo
    bi = Bm_inactive + Bo_inactive
    bt_intact = bl + bi
    bl_next = Bm_next + Bo_next
    bi_next = Bm_inactive_next + Bo_inactive_next
    bt_intact_next = bl_next + bi_next
    if not rec["_started"]:
        rec["_started"] = True
        rec["first_BM"] = float(Bm)
        rec["first_BO"] = float(Bo)
        rec["first_BM_inactive"] = float(Bm_inactive)
        rec["first_BO_inactive"] = float(Bo_inactive)
        rec["first_BL"] = float(bl)
        rec["first_BI"] = float(bi)
        rec["first_BT_intact"] = float(bt_intact)
    rec["n_grid_hours"] += 1.0
    rec["sum_BM"] += float(Bm)
    rec["sum_BO"] += float(Bo)
    rec["sum_BM_inactive"] += float(Bm_inactive)
    rec["sum_BO_inactive"] += float(Bo_inactive)
    rec["sum_BL"] += float(bl)
    rec["sum_BI"] += float(bi)
    rec["sum_BT_intact"] += float(bt_intact)
    rec["last_BM"] = float(Bm_next)
    rec["last_BO"] = float(Bo_next)
    rec["last_BM_inactive"] = float(Bm_inactive_next)
    rec["last_BO_inactive"] = float(Bo_inactive_next)
    rec["last_BL"] = float(bl_next)
    rec["last_BI"] = float(bi_next)
    rec["last_BT_intact"] = float(bt_intact_next)
    rec["sum_micro_mu_max"] += float(Bm * mu_max)
    rec["sum_micro_mu_committed"] += float(Bm * mu_comm)
    rec["sum_micro_density_gate"] += float(Bm * density_micro)
    rec["sum_micro_active_gate"] += float(Bm * active)
    rec["sum_micro_rho_gate"] += float(Bm * safe_div(mu_comm, mu_max))
    rec["sum_micro_realized"] += float(Bm * mu_micro)
    rec["sum_other_mu_max"] += float(Bo * other_mu_max)
    rec["sum_other_temp_gate"] += float(Bo * temp_gate)
    rec["sum_other_light_gate"] += float(Bo * light_gate)
    rec["sum_other_N_gate"] += float(Bo * n_gate)
    rec["sum_other_P_gate"] += float(Bo * p_gate)
    rec["sum_other_nutrient_gate"] += float(Bo * nutrient_gate)
    rec["sum_other_env_product"] += float(Bo * env_product)
    rec["sum_other_env_single_limiter"] += float(Bo * env_single)
    rec["sum_other_density_gate"] += float(Bo * density_other)
    rec["sum_other_realized"] += float(Bo * mu_other)
    rec["sum_BT_phi_micro"] += float(bl * phi_micro)
    rec["sum_BT_phi_other"] += float(bl * phi_other)
    rec["sum_BT_phi_ratio_micro_over_other"] += float(bl * phi_ratio_micro_over_other)
    rec["sum_BT_ln_phi_micro_over_other"] += float(bl * ln_phi_micro_over_other)
    rec["sum_BT_phi_micro_lt_other"] += float(bl * (1.0 if phi_micro < phi_other else 0.0))
    rec["sum_BT_phi_ratio_lt_0p5"] += float(bl * (1.0 if phi_ratio_micro_over_other < 0.5 else 0.0))
    rec["sum_BT_phi_ratio_gt_2"] += float(bl * (1.0 if phi_ratio_micro_over_other > 2.0 else 0.0))
    rec["sum_BT_w_micro_allocation"] += float(bl * w_micro_allocation)
    rec["sum_BT_w_micro_phi_pure"] += float(bl * w_micro_phi_pure)
    rec["sum_BT_w_minus_p_micro"] += float(bl * w_minus_p_micro)
    rec["sum_BT_logit_w_minus_logit_p"] += float(bl * logit_w_minus_logit_p)
    rec["sum_BT_L_B_other_over_micro"] += float(bl * L_B_other_over_micro)
    rec["sum_BT_L_phi_other_over_micro"] += float(bl * L_phi_other_over_micro)
    rec["sum_BT_L_w_other_over_micro_gamma"] += float(bl * L_w_other_over_micro_gamma)
    rec["sum_BT_logit_w_micro_gamma"] += float(bl * logit_w_micro_gamma)
    rec["sum_BT_succession_signal"] += float(bl * succession_signal)
    rec["sum_BT_succession_pressure"] += float(bl * succession_pressure)
    rec["sum_BT_suitability_relative_advantage"] += float(bl * suitability_relative_advantage)
    rec["sum_GT_ref_flux"] += float(GT_ref)
    rec["sum_micro_growth_ref_flux"] += float(micro_growth_ref)
    rec["sum_other_growth_ref_flux"] += float(other_growth_ref)
    rec["sum_micro_growth_allocation_delta_flux"] += float(micro_growth - micro_growth_ref)
    rec["sum_other_growth_allocation_delta_flux"] += float(other_growth - other_growth_ref)
    rec["sum_allocation_residual_flux"] += float(allocation_residual)
    rec["sum_micro_capacity_flux"] += float(micro_capacity_flux)
    rec["sum_other_capacity_flux"] += float(other_capacity_flux)
    rec["sum_other_capacity_basis"] += float(other_capacity_basis)
    rec["sum_other_capacity_basis_env_multiplier"] += float(other_capacity_basis * other_capacity_env_multiplier)
    rec["sum_other_recruitable_capacity_basis"] += float(other_recruitable_capacity_basis)
    rec["sum_micro_capacity_slack_flux"] += float(micro_capacity_slack_flux)
    rec["sum_other_capacity_slack_flux"] += float(other_capacity_slack_flux)
    rec["sum_projection_distance_flux"] += float(projection_distance_flux)
    rec["sum_unallocatable_budget_flux"] += float(unallocatable_budget_flux)
    rec["sum_target_infeasible_flux"] += float(0.0 if target_feasible else GT_ref)
    rec["sum_capacity_budget_inconsistent_flux"] += float(GT_ref if capacity_budget_inconsistent else 0.0)
    rec["sum_micro_capacity_binding_hours"] += float(1.0 if micro_capacity_binding else 0.0)
    rec["sum_other_capacity_binding_hours"] += float(1.0 if other_capacity_binding else 0.0)
    rec["sum_target_feasible_hours"] += float(1.0 if target_feasible else 0.0)
    rec["sum_micro_growth_flux"] += float(micro_growth)
    rec["sum_other_growth_flux"] += float(other_growth)
    rec["sum_micro_loss_flux"] += float(micro_loss)
    rec["sum_other_loss_flux"] += float(other_loss)
    rec["sum_composition_turnover_flux"] += float(composition_turnover_flux)
    rec["sum_micro_inactivation_flux"] += float(micro_inactivation_flux)
    rec["sum_other_inactivation_flux"] += float(other_inactivation_flux)
    rec["sum_micro_inactive_lysis_flux"] += float(micro_inactive_lysis_flux)
    rec["sum_other_inactive_lysis_flux"] += float(other_inactive_lysis_flux)
    rec["sum_micro_inactive_settling_flux"] += float(micro_inactive_settling_flux)
    rec["sum_other_inactive_settling_flux"] += float(other_inactive_settling_flux)
    rec["sum_inactive_module_residual_flux"] += float(inactive_module_residual_flux)
    rec["sum_BM_late_memory"] += float(Bm * late_memory)
    rec["sum_BM_late_damage"] += float(Bm * late_damage)
    rec["sum_BM_late_raw_pressure"] += float(Bm * late_raw_pressure)
    rec["sum_BM_late_decline_signal"] += float(Bm * late_decline_signal)
    rec["sum_BM_late_saturation_gate"] += float(Bm * late_saturation_gate)
    rec["sum_BM_late_exhaustion_gate"] += float(Bm * late_exhaustion_gate)
    rec["sum_BM_micro_inactivation_rate_h"] += float(Bm * micro_inactivation_rate_h)
    rec["sum_BM_surface_factor"] += float(Bm * surface_factor)
    rec["sum_BM_surface_equiv"] += float(Bm * surface_factor)
    rec["sum_BT_surface_equiv"] += float(bl * surface_factor)
    rec["sum_negative_correction_mass"] += float(neg_corr)
    rec["sum_state_balance_residual"] += float(residual)

def rows_from_accums(accums: dict[tuple[str, str, str, str, str], dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for (scenario, scope, scope_id, date, phase), rec in accums.items():
        out = {c: float(rec.get(c, 0.0)) for c in ACCUM_COLS}
        out.update({
            "scenario": scenario,
            "scope": scope,
            "scope_id": scope_id,
            "date": date,
            "phase": phase,
        })
        rows.append(out)
    return rows

def simulate_grid_coupled_generator(
    gdf: pd.DataFrame,
    init_row: dict[str, Any],
    rs: Any,
    cfg: Any,
    station_names: list[str],
    source_state_scenario: str | None = None,
    collect_accums: bool = False,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    gdf = gdf.sort_values("time").reset_index(drop=True)
    times = pd.to_datetime(gdf["time"], utc=True)
    if len(times) > 1:
        dt = times.diff().dt.total_seconds().div(3600).fillna(1.0).clip(lower=1e-6).to_numpy()
    else:
        dt = np.array([1.0])

    temp_arr = pd.to_numeric(gdf["temp_C"], errors="coerce").to_numpy(float)
    light_arr = pd.to_numeric(gdf["light_Wm2"], errors="coerce").to_numpy(float)
    DIN_arr = pd.to_numeric(gdf["DIN_mgN_L"], errors="coerce").to_numpy(float)
    DIP_arr = pd.to_numeric(gdf["DIP_mgP_L"], errors="coerce").to_numpy(float)
    wind_arr = pd.to_numeric(gdf["wind_speed_mps"], errors="coerce").fillna(0.0).to_numpy(float) if "wind_speed_mps" in gdf.columns else np.zeros(len(gdf))
    lon = float(pd.to_numeric(gdf["lon"], errors="coerce").dropna().iloc[0]) if "lon" in gdf and gdf["lon"].notna().any() else np.nan
    lat = float(pd.to_numeric(gdf["lat"], errors="coerce").dropna().iloc[0]) if "lat" in gdf and gdf["lat"].notna().any() else np.nan
    depth = float(pd.to_numeric(gdf["depth_m"], errors="coerce").dropna().iloc[0]) if "depth_m" in gdf and gdf["depth_m"].notna().any() else np.nan
    basin = basin_from_lon(lon)

    states = {name: state_from_initial(init_row, cfg) for name in SCENARIOS}
    accums: dict[tuple[str, str, str, str, str], dict[str, Any]] = {}
    source_rows: list[dict[str, Any]] = []

    for k in range(len(gdf)):
        date = str(times.iloc[k].date())
        phase = phase_from_date(date)
        dth = float(dt[k])
        temp = float(temp_arr[k])
        light = max(float(light_arr[k]), 0.0)
        DIN = max(float(DIN_arr[k]), 0.0)
        DIP = max(float(DIP_arr[k]), 0.0)
        wind = float(wind_arr[k]) if np.isfinite(wind_arr[k]) else 0.0

        for scenario, meta in SCENARIOS.items():
            st = states[scenario]
            Bm = max(float(st["Bm"]), 0.0)
            Bo = max(float(st["Bo"]), 0.0)
            Bm_inactive = max(float(st.get("Bm_inactive", 0.0)), 0.0)
            Bo_inactive = max(float(st.get("Bo_inactive", 0.0)), 0.0)
            ftox = float(np.clip(st.get("ftox", cfg.ftox_default), cfg.ftox_min, cfg.ftox_max))
            ftox_gene_proxy_raw = float(st.get("ftox_gene_proxy_raw", np.nan))
            active_state = float(np.clip(st["active_state"], cfg.active_min_fraction, 1.0))
            allocation_state = float(np.clip(st["allocation_state"], 0.0, 1.0))
            release_pressure_state = float(np.clip(st["release_pressure_state"], 0.0, 1.0))
            low_light_memory = float(np.clip(st["low_light_memory"], 0.0, 1.0))
            aggregation_state = float(np.clip(st.get("aggregation_state", 0.0), 0.0, 1.0))
            high_biomass_memory_state = float(np.clip(st.get("high_biomass_memory_state", 0.0), 0.0, 1.0))
            dominance_memory_state = float(np.clip(st.get("dominance_memory_state", 0.0), 0.0, 1.0))
            bloom_patch_memory_state = float(np.clip(st.get("bloom_patch_memory_state", 0.0), 0.0, 1.0))
            scum_accumulation_memory_state = float(np.clip(st.get("scum_accumulation_memory_state", 0.0), 0.0, 1.0))
            succession_pressure_state = float(np.clip(st.get("succession_pressure_state", 0.0), 0.0, 1.0))
            late_inactivation_memory_state = float(np.clip(st.get("late_inactivation_memory_state", 0.0), 0.0, 1.0))
            late_viability_burden_h = max(float(st.get("late_viability_burden_h", 0.0)), 0.0)
            late_damage_burden_state = max(float(st.get("late_damage_burden_state", 0.0)), 0.0)
            mu_capacity_memory_h = float(st.get("mu_capacity_memory_h", np.nan))
            dynamic_capacity_memory_fraction = float(np.clip(
                st.get("dynamic_capacity_memory_fraction", 1.0), 0.0, 1.0
            ))
            growth_competence_state = float(np.clip(st.get("growth_competence_state", meta.get("rgc_initial_state", 1.0)), 0.0, 1.0))

            n_gate = sat(DIN, cfg.K_N_mgN_L)
            p_gate = sat(DIP, cfg.K_P_mgP_L)
            nutrient_gate = min(n_gate, p_gate)
            light_suit = sat(light, cfg.light_stress_half_Wm2)
            light_stress = 1.0 - light_suit
            other_light = sat(light, cfg.other_light_half_Wm2)
            temp_gate = temp_suitability(temp, cfg.temp_opt_C, cfg.temp_width_C)
            low_light_memory += dth * (light_stress - low_light_memory) / max(cfg.low_light_memory_tau_h, EPS)
            low_light_memory = float(np.clip(low_light_memory, 0.0, 1.0))

            inactive_chi = float(getattr(cfg, "inactive_density_contribution_chi", 0.0))
            Bm_density_equiv = Bm + inactive_chi * Bm_inactive
            Bo_density_equiv = Bo + inactive_chi * Bo_inactive
            density_feedback_enabled = bool(meta.get("density_feedback_enabled", True))
            if density_feedback_enabled:
                density_pressure = max(0.0, (Bm_density_equiv + cfg.alpha_other_on_micro * Bo_density_equiv) / max(cfg.carrying_capacity_state, EPS))
                density_micro = max(0.0, 1.0 - density_pressure)
                density_other = max(0.0, 1.0 - (Bo_density_equiv + cfg.alpha_micro_on_other * Bm_density_equiv) / max(cfg.carrying_capacity_state, EPS))
            else:
                density_pressure = 0.0
                density_micro = 1.0
                density_other = 1.0
            dominance = Bm / max(Bm + Bo, EPS)
            competition_pressure = float(np.clip(Bo / max(Bm + Bo, EPS), 0.0, 1.0))
            maturity_state = float(np.clip(
                0.60 * sigmoid(density_pressure, 0.45, 0.18)
                + 0.25 * sigmoid(dominance, 0.55, 0.15)
                + 0.15 * sigmoid(Bm / max(cfg.carrying_capacity_state, EPS), 0.55, 0.20),
                0.0, 1.0,
            ))
            env_nutrient_stress = 1.0 - nutrient_gate
            if bool(getattr(cfg, "enable_patch_memory", True)):
                high_biomass_signal = sigmoid(Bm, cfg.patch_B_micro_ref_state, cfg.patch_B_micro_width_state)
                dominance_signal = sigmoid(dominance, cfg.patch_dominance_ref, cfg.patch_dominance_width)
                patch_instant_signal = float(np.clip(0.55 * high_biomass_signal + 0.30 * dominance_signal + 0.15 * maturity_state, 0.0, 1.0))
                stability_for_memory = 1.0 / (1.0 + max(wind, 0.0) / max(cfg.aggregation_wind_half_mps, EPS))
                scum_instant_signal = float(np.clip(
                    bloom_patch_memory_state * stability_for_memory * (0.50 + 0.50 * release_pressure_state) * (0.50 + 0.50 * maturity_state),
                    0.0, 1.0,
                ))
                high_biomass_memory_state += dth * (high_biomass_signal - high_biomass_memory_state) / max(cfg.high_biomass_memory_tau_h, EPS)
                dominance_memory_state += dth * (dominance_signal - dominance_memory_state) / max(cfg.dominance_memory_tau_h, EPS)
                bloom_patch_memory_state += dth * (patch_instant_signal - bloom_patch_memory_state) / max(cfg.bloom_patch_memory_tau_h, EPS)
                scum_accumulation_memory_state += dth * (scum_instant_signal - scum_accumulation_memory_state) / max(cfg.scum_memory_tau_h, EPS)
                high_biomass_memory_state = float(np.clip(high_biomass_memory_state, 0.0, 1.0))
                dominance_memory_state = float(np.clip(dominance_memory_state, 0.0, 1.0))
                bloom_patch_memory_state = float(np.clip(bloom_patch_memory_state, 0.0, 1.0))
                scum_accumulation_memory_state = float(np.clip(scum_accumulation_memory_state, 0.0, 1.0))
            else:
                high_biomass_memory_state = dominance_memory_state = bloom_patch_memory_state = scum_accumulation_memory_state = 0.0
            allocation_target = float(np.clip(
                cfg.alloc_stress_weight * env_nutrient_stress
                + cfg.alloc_maturity_weight * maturity_state
                + cfg.alloc_competition_weight * competition_pressure
                + cfg.alloc_lowlight_weight * low_light_memory,
                0.0, 1.0,
            ))
            allocation_state += dth * (allocation_target - allocation_state) / max(cfg.alloc_tau_h, EPS)
            allocation_state = float(np.clip(allocation_state, 0.0, 1.0))
            rho = float(np.clip(
                cfg.rho_growth_high - allocation_state * (cfg.rho_growth_high - cfg.rho_growth_low),
                min(cfg.rho_growth_low, cfg.rho_growth_high),
                max(cfg.rho_growth_low, cfg.rho_growth_high),
            ))
            resp = rs.lookup(temp_C=temp, light_Wm2=light, DIN_mgN_L=DIN, DIP_mgP_L=DIP, rho_growth_commitment=rho, k_env=cfg.response_k_env)
            mu_max = max(float(resp["mu_max_h"]), 0.0)
            mu_comm = max(float(resp["mu_committed_h"]), 0.0)
            cap_upper_unit = max(float(resp["toxin_capacity_upper_per_gDW_h"]), 0.0)

            env_product = temp_gate * other_light * nutrient_gate
            env_single = temp_gate * min(other_light, n_gate, p_gate)

            # R17 is enabled only through scenario metadata. Capacity responds
            # to daylight ecModel decline relative to the existing 72 h memory;
            # darkness alone must not collapse K_eff.
            dynamic_capacity_enabled = bool(meta.get("dynamic_capacity_enabled", False))
            dynamic_k_feedback_enabled = bool(meta.get(
                "dynamic_k_feedback_enabled", dynamic_capacity_enabled
            ))
            replacement_growth_enabled = bool(meta.get(
                "replacement_growth_enabled", dynamic_capacity_enabled
            ))
            dynamic_excess_exit_enabled = bool(meta.get(
                "dynamic_excess_exit_enabled", dynamic_capacity_enabled
            ))
            dynamic_capacity_ratio = 1.0
            dynamic_capacity_target_fraction = 1.0
            dynamic_capacity_fraction = 1.0
            dynamic_capacity_state = float(cfg.carrying_capacity_state)
            dynamic_excess_rate_h = 0.0
            if dynamic_capacity_enabled:
                if light >= float(cfg.late_decline_light_min_Wm2):
                    dynamic_capacity_target_fraction = float(np.clip(
                        safe_div(mu_max, cfg.micro_suitability_ref_h, default=1.0), 0.0, 1.0
                    ))
                    dynamic_capacity_memory_fraction += dth * (
                        dynamic_capacity_target_fraction - dynamic_capacity_memory_fraction
                    ) / max(float(meta.get("dynamic_capacity_tau_h", cfg.late_mu_memory_tau_h)), EPS)
                    dynamic_capacity_memory_fraction = float(np.clip(
                        dynamic_capacity_memory_fraction, 0.0, 1.0
                    ))
                dynamic_capacity_ratio = dynamic_capacity_memory_fraction
                floor_fraction = float(np.clip(meta.get("dynamic_capacity_floor_fraction", 0.5), 0.0, 1.0))
                dynamic_capacity_fraction = floor_fraction + (1.0 - floor_fraction) * dynamic_capacity_ratio
                dynamic_capacity_state = max(float(cfg.carrying_capacity_state) * dynamic_capacity_fraction, EPS)
                if density_feedback_enabled:
                    density_denominator = (
                        dynamic_capacity_state
                        if dynamic_k_feedback_enabled
                        else max(float(cfg.carrying_capacity_state), EPS)
                    )
                    density_pressure = max(
                        0.0,
                        (Bm_density_equiv + cfg.alpha_other_on_micro * Bo_density_equiv) / density_denominator,
                    )
                    density_micro = max(0.0, 1.0 - density_pressure)
                    density_other = max(
                        0.0,
                        1.0 - (Bo_density_equiv + cfg.alpha_micro_on_other * Bm_density_equiv) / density_denominator,
                    )

            micro_potential_flux = max(0.0, mu_max * Bm)
            other_potential_flux = max(0.0, cfg.other_mu_max_h * env_product * Bo)
            if replacement_growth_enabled:
                micro_replacement_flux = min(micro_potential_flux, cfg.micro_base_loss_rate_h * Bm)
                other_replacement_flux = min(other_potential_flux, cfg.other_loss_rate_h * Bo)
                micro_growth_ref = micro_replacement_flux + density_micro * max(
                    0.0, micro_potential_flux - micro_replacement_flux
                )
                other_growth_ref = other_replacement_flux + density_other * max(
                    0.0, other_potential_flux - other_replacement_flux
                )
            else:
                micro_replacement_flux = 0.0
                other_replacement_flux = 0.0
                micro_growth_ref = max(0.0, mu_max * density_micro * Bm)
                other_growth_ref = max(0.0, cfg.other_mu_max_h * env_product * density_other * Bo)
            GT_ref_native = micro_growth_ref + other_growth_ref
            rgc_enabled = bool(meta.get("rgc_enabled", False))
            if rgc_enabled:
                rgc_live_basis = max(Bm + Bo, EPS)
                rgc_loss_ref = max(
                    0.0,
                    cfg.micro_base_loss_rate_h * Bm + cfg.other_loss_rate_h * Bo,
                )
                rgc_mu_potential_h = safe_div(GT_ref_native, rgc_live_basis, default=0.0)
                rgc_mu_loss_h = safe_div(rgc_loss_ref, rgc_live_basis, default=0.0)
                growth_competence_eq_state = float(np.clip(
                    safe_div(rgc_mu_potential_h, rgc_mu_potential_h + rgc_mu_loss_h + EPS, default=1.0),
                    0.0,
                    1.0,
                ))
                if growth_competence_eq_state > growth_competence_state:
                    rgc_tau_h = max(float(meta.get("rgc_tau_recovery_h", 72.0)), EPS)
                else:
                    rgc_tau_h = max(float(meta.get("rgc_tau_deactivation_h", 48.0)), EPS)
                growth_competence_state += dth * (growth_competence_eq_state - growth_competence_state) / rgc_tau_h
                growth_competence_state = float(np.clip(growth_competence_state, 0.0, 1.0))
                GT_ref = growth_competence_state * GT_ref_native
            else:
                growth_competence_state = 1.0
                growth_competence_eq_state = 1.0
                GT_ref = GT_ref_native

            phi_micro_current = float(np.clip(safe_div(mu_comm, mu_max, default=0.0), 0.0, 1.0))
            phi_other_product = float(np.clip(env_product, 0.0, 1.0))
            phi_other_bottleneck = float(np.clip(min(temp_gate, other_light, n_gate, p_gate), 0.0, 1.0))
            phi_other_bridge = float(np.clip(np.sqrt(max(phi_other_product, 0.0) * max(phi_other_bottleneck, 0.0)), 0.0, 1.0))

            phi_m_mode = str(meta.get("phi_M_mode", "current"))
            if phi_m_mode == "current":
                phi_micro = phi_micro_current
            elif phi_m_mode == "current_times_active":
                phi_micro = float(np.clip(phi_micro_current * active_state, 0.0, 1.0))
            elif phi_m_mode == "current_times_sqrt_active":
                phi_micro = float(np.clip(phi_micro_current * np.sqrt(max(active_state, 0.0)), 0.0, 1.0))
            else:
                raise ValueError(f"Unknown phi_M_mode={phi_m_mode!r}")

            phi_o_mode = str(meta.get("phi_O_mode", "product"))
            if phi_o_mode == "product":
                phi_other = phi_other_product
            elif phi_o_mode == "bottleneck":
                phi_other = phi_other_bottleneck
            elif phi_o_mode == "product_bottleneck_geomean":
                phi_other = phi_other_bridge
            else:
                raise ValueError(f"Unknown phi_O_mode={phi_o_mode!r}")

            denom_phi = Bm * phi_micro + Bo * phi_other
            w_phi = safe_div(Bm * phi_micro, denom_phi)
            if not np.isfinite(w_phi):
                w_phi = safe_div(Bm, Bm + Bo, default=0.0)
            w_phi_pure = safe_div(phi_micro, phi_micro + phi_other, default=0.0)
            p_live = float(np.clip(safe_div(Bm, Bm + Bo, default=0.0), 0.0, 1.0))
            phi_ratio_micro_over_other = safe_div(phi_micro, phi_other, default=np.nan)
            if not np.isfinite(phi_ratio_micro_over_other):
                phi_ratio_micro_over_other = 0.0 if phi_micro <= 0.0 and phi_other > 0.0 else 1e12
            ln_phi_micro_over_other = float(np.log((phi_micro + EPS) / (phi_other + EPS)))
            L_B_other_over_micro = float(np.log((Bo + EPS) / (Bm + EPS)))
            L_phi_other_over_micro = float(np.log((phi_other + EPS) / (phi_micro + EPS)))

            micro_potential_h = mu_max * density_micro
            other_realized_potential_h = cfg.other_mu_max_h * env_product * density_other
            micro_capacity_flux = max(0.0, micro_growth_ref)
            other_capacity_mode = str(meta.get("other_capacity_mode", "env_single_existing_biomass"))
            other_capacity_recruitable_fraction = max(float(meta.get("other_capacity_recruitable_fraction", 0.0)), 0.0)
            other_capacity_env_power = max(float(meta.get("other_capacity_env_power", 1.0)), 0.0)
            other_recruitable_capacity_basis = 0.0
            if other_capacity_mode == "env_single_existing_biomass":
                other_capacity_basis = Bo
                other_capacity_env_multiplier = env_single
            elif other_capacity_mode == "env_power_existing_biomass":
                other_capacity_basis = Bo
                other_capacity_env_multiplier = float(np.clip(env_single, 0.0, 1.0)) ** other_capacity_env_power
            elif other_capacity_mode == "intrinsic_existing_biomass":
                other_capacity_basis = Bo
                other_capacity_env_multiplier = 1.0
            elif other_capacity_mode == "intrinsic_recruitable_biomass":
                other_recruitable_capacity_basis = other_capacity_recruitable_fraction * (Bm + Bo)
                other_capacity_basis = Bo + other_recruitable_capacity_basis
                other_capacity_env_multiplier = 1.0
            else:
                raise ValueError(f"Unknown other_capacity_mode={other_capacity_mode!r}")
            other_capacity_flux = max(
                0.0,
                (
                    min(
                        cfg.other_mu_max_h * other_capacity_env_multiplier * other_capacity_basis,
                        cfg.other_loss_rate_h * other_capacity_basis,
                    )
                    + density_other * max(
                        0.0,
                        cfg.other_mu_max_h * other_capacity_env_multiplier * other_capacity_basis
                        - cfg.other_loss_rate_h * other_capacity_basis,
                    )
                )
                if replacement_growth_enabled
                else cfg.other_mu_max_h * other_capacity_env_multiplier * density_other * other_capacity_basis,
            )
            S_micro = float(np.clip(safe_div(micro_potential_h, cfg.micro_suitability_ref_h), 0.0, 1.0))
            S_other = float(np.clip(safe_div(other_realized_potential_h, cfg.other_mu_max_h), 0.0, 1.0))

            allocation_selectivity_enabled = bool(meta.get("allocation_selectivity_enabled", True))
            mode = str(meta["allocation_mode"])
            if mode == "matched_phi_allocation":
                gamma = float(meta.get("gamma_inventory", 1.0))
                if not allocation_selectivity_enabled:
                    w_micro_allocation = p_live
                elif abs(gamma - 1.0) <= 1e-12:
                    w_micro_allocation = float(np.clip(w_phi, 0.0, 1.0))
                else:
                    p_gamma = float(np.clip(p_live, 0.0, 1.0)) ** gamma
                    o_gamma = float(np.clip(1.0 - p_live, 0.0, 1.0)) ** gamma
                    denom_gamma = p_gamma * phi_micro + o_gamma * phi_other
                    w_micro_allocation = safe_div(p_gamma * phi_micro, denom_gamma, default=p_live)
                    w_micro_allocation = float(np.clip(w_micro_allocation, 0.0, 1.0))
                micro_growth = w_micro_allocation * GT_ref
                other_growth = (1.0 - w_micro_allocation) * GT_ref
                micro_capacity_slack_flux = max(0.0, micro_capacity_flux - micro_growth)
                other_capacity_slack_flux = max(0.0, other_capacity_flux - other_growth)
                projection_distance_flux = 0.0
                unallocatable_budget_flux = 0.0
                target_feasible = True
                capacity_budget_inconsistent = False
                micro_capacity_binding = bool(micro_capacity_flux - micro_growth <= 1e-12 and micro_growth > 1e-12)
                other_capacity_binding = bool(other_capacity_flux - other_growth <= 1e-12 and other_growth > 1e-12)
            elif mode == "conservative_projected_phi_allocation":
                gamma = float(meta.get("gamma_inventory", 1.0))
                if not allocation_selectivity_enabled:
                    w_micro_allocation = p_live
                else:
                    p_gamma = float(np.clip(p_live, 0.0, 1.0)) ** gamma
                    o_gamma = float(np.clip(1.0 - p_live, 0.0, 1.0)) ** gamma
                    denom_gamma = p_gamma * phi_micro + o_gamma * phi_other
                    w_micro_allocation = safe_div(p_gamma * phi_micro, denom_gamma, default=p_live)
                    w_micro_allocation = float(np.clip(w_micro_allocation, 0.0, 1.0))
                projection = project_two_species_conservative_production(
                    B_M=Bm,
                    B_O=Bo,
                    P_total_native=GT_ref,
                    w_M_target=w_micro_allocation,
                    C_M=micro_capacity_flux,
                    C_O=other_capacity_flux,
                    tolerance=1e-12,
                )
                micro_growth = float(projection.P_M_projected)
                other_growth = float(projection.P_O_projected)
                micro_capacity_slack_flux = float(projection.M_capacity_slack)
                other_capacity_slack_flux = float(projection.O_capacity_slack)
                projection_distance_flux = float(projection.projection_distance)
                unallocatable_budget_flux = float(projection.unallocatable_budget)
                target_feasible = bool(projection.target_feasible)
                capacity_budget_inconsistent = bool(projection.capacity_budget_inconsistent)
                micro_capacity_binding = bool(projection.M_capacity_binding)
                other_capacity_binding = bool(projection.O_capacity_binding)
            else:
                raise ValueError(f"Unknown allocation_mode={mode!r}")
            L_w_other_over_micro_gamma = gamma * L_B_other_over_micro + L_phi_other_over_micro
            logit_w_micro_gamma = -L_w_other_over_micro_gamma
            logit_p = float(np.log((p_live + EPS) / (1.0 - p_live + EPS)))
            w_minus_p_micro = w_micro_allocation - p_live
            logit_w_minus_logit_p = logit_w_micro_gamma - logit_p

            mu_micro = safe_div(micro_growth, Bm, default=0.0)
            mu_other = safe_div(other_growth, Bo, default=0.0)
            allocation_residual = (micro_growth + other_growth) - GT_ref

            if not np.isfinite(mu_capacity_memory_h):
                mu_capacity_memory_h = max(mu_max, 0.0)
            mu_capacity_memory_before_h = float(mu_capacity_memory_h)
            late_decline_daylight_gate = 1.0 if light >= float(cfg.late_decline_light_min_Wm2) else 0.0
            if late_decline_daylight_gate > 0.0:
                mu_capacity_memory_h += dth * (max(mu_max, 0.0) - mu_capacity_memory_h) / max(cfg.late_mu_memory_tau_h, EPS)
                mu_capacity_memory_h = float(max(mu_capacity_memory_h, 0.0))
                mu_capacity_decline_per_day = 24.0 * max(
                    0.0,
                    (mu_capacity_memory_before_h - mu_capacity_memory_h) / max(mu_capacity_memory_before_h, EPS) / max(dth, EPS),
                )
            else:
                mu_capacity_decline_per_day = 0.0
            late_saturation_pressure = float(np.clip(1.0 - density_micro, 0.0, 1.0)) ** max(float(cfg.late_saturation_power), 0.0)
            late_saturation_excess_state = max(0.0, late_saturation_pressure - float(cfg.late_saturation_ref))
            late_saturation_gate_state = sat(late_saturation_excess_state, cfg.late_saturation_width)
            late_potential_net_per_day = 24.0 * (mu_micro - float(cfg.micro_base_loss_rate_h))
            late_exhaustion_excess_state = max(0.0, float(cfg.late_exhaustion_rcrit_per_day) - late_potential_net_per_day)
            late_exhaustion_gate_state = sat(late_exhaustion_excess_state, cfg.late_exhaustion_width_per_day)
            late_decline_signal = sat(mu_capacity_decline_per_day, cfg.late_decline_half_per_day)
            late_inactivation_signal_mode = str(getattr(cfg, "late_inactivation_signal_mode", "metabolic_decline_saturation") or "metabolic_decline_saturation").strip().lower()
            late_inactivation_enabled = bool(meta.get(
                "late_inactivation_enabled",
                bool(getattr(cfg, "enable_late_inactivation_pool", False)),
            ))
            if late_inactivation_enabled and late_inactivation_signal_mode in {"metabolic_decline_saturation_exhaustion_damage", "damage_memory", "d3_3a"}:
                late_raw_viability_pressure_state = float(np.clip(
                    late_saturation_gate_state * late_decline_signal * late_exhaustion_gate_state,
                    0.0, 1.0,
                ))
                late_damage_burden_state += dth * (
                    float(cfg.late_damage_accum_rate_h) * late_raw_viability_pressure_state
                    - late_damage_burden_state / max(float(cfg.late_damage_repair_tau_h), EPS)
                )
                late_damage_burden_state = float(max(0.0, late_damage_burden_state))
                late_inactivation_signal_state = sat(late_damage_burden_state, cfg.late_damage_half_saturation)
            elif late_inactivation_enabled and late_inactivation_signal_mode in {"metabolic_decline_saturation_exhaustion", "and_gate", "d3_2c"}:
                late_raw_viability_pressure_state = float(np.clip(
                    late_saturation_gate_state * late_decline_signal * late_exhaustion_gate_state,
                    0.0, 1.0,
                ))
                late_inactivation_signal_state = late_raw_viability_pressure_state
                late_damage_burden_state = 0.0
            elif late_inactivation_enabled and late_inactivation_signal_mode in {"metabolic_decline_saturation", "decline_saturation", "mu_decline"}:
                late_raw_viability_pressure_state = float(np.clip(late_saturation_pressure * late_decline_signal, 0.0, 1.0))
                late_inactivation_signal_state = late_raw_viability_pressure_state
                late_damage_burden_state = 0.0
            elif (not late_inactivation_enabled) or late_inactivation_signal_mode in {"off", "none", "disabled"}:
                late_raw_viability_pressure_state = 0.0
                late_inactivation_signal_state = 0.0
                late_damage_burden_state = 0.0
            else:
                raise ValueError(f"Unknown late_inactivation_signal_mode={late_inactivation_signal_mode!r}")
            late_viability_burden_h = 0.0
            late_inactivation_memory_state += dth * (late_inactivation_signal_state - late_inactivation_memory_state) / max(cfg.late_inactivation_memory_tau_h, EPS)
            late_inactivation_memory_state = float(np.clip(late_inactivation_memory_state, 0.0, 1.0))

            release_pressure_target = float(np.clip(
                0.45 * sigmoid(density_pressure, 0.85, 0.20)
                + 0.35 * env_nutrient_stress
                + 0.20 * low_light_memory,
                0.0, 1.0,
            ))
            release_pressure_state += dth * (release_pressure_target - release_pressure_state) / max(cfg.release_pressure_tau_h, EPS)
            release_pressure_state = float(np.clip(release_pressure_state, 0.0, 1.0))
            active_capacity_index = float(np.clip(mu_comm / max(cfg.active_mu_ref_h, EPS), cfg.active_min_fraction, 1.0))
            active_release_factor = 1.0 - cfg.active_release_pressure_weight * release_pressure_state
            active_lowlight_factor = 1.0 - cfg.active_low_light_weight * low_light_memory
            active_target = float(np.clip(
                active_capacity_index * active_release_factor * active_lowlight_factor,
                cfg.active_min_fraction,
                1.0,
            ))
            active_state += dth * (active_target - active_state) / max(cfg.active_tau_h, EPS)
            active_state = float(np.clip(active_state, cfg.active_min_fraction, 1.0))

            ftox_factor = (
                1.0
                - cfg.ftox_growth_competition_sensitivity
                * competition_pressure
                * dth
                / max(cfg.ftox_tau_h, EPS)
            )
            ftox_target = float(np.clip(ftox * ftox_factor, cfg.ftox_min, cfg.ftox_max))
            ftox += dth * (ftox_target - ftox) / max(cfg.ftox_tau_h, EPS)
            ftox = float(np.clip(ftox, cfg.ftox_min, cfg.ftox_max))
            if np.isfinite(ftox_gene_proxy_raw):
                raw_target = max(1e-12, ftox_gene_proxy_raw * ftox_factor)
                ftox_gene_proxy_raw += dth * (raw_target - ftox_gene_proxy_raw) / max(cfg.ftox_tau_h, EPS)
                ftox_gene_proxy_raw = float(max(1e-12, ftox_gene_proxy_raw))

            X_micro = cfg.biomass_scale_gDW_m3_per_state * Bm
            raw_mcye_effective_fraction = float(np.clip(ftox, cfg.ftox_min, cfg.ftox_max))
            patch_memory_toxic_prior_increment = cfg.patch_memory_toxic_prior_boost_fraction * bloom_patch_memory_state
            toxic_producer_prior_state = float(np.clip(
                cfg.toxic_producer_prior_fraction
                + cfg.toxic_prior_maturity_boost_fraction * maturity_state
                + cfg.toxic_prior_release_boost_fraction * release_pressure_state * maturity_state
                + patch_memory_toxic_prior_increment,
                cfg.toxic_producer_prior_min,
                cfg.toxic_producer_prior_max,
            ))
            bloom_stage_gate = float(np.clip(
                cfg.bloom_stage_gate_min
                + (cfg.bloom_stage_gate_max - cfg.bloom_stage_gate_min) * maturity_state,
                cfg.bloom_stage_gate_min,
                cfg.bloom_stage_gate_max,
            ))
            release_maturity_gate = float(np.clip(
                cfg.release_maturity_gate_min
                + (cfg.release_maturity_gate_max - cfg.release_maturity_gate_min) * release_pressure_state * maturity_state,
                cfg.release_maturity_gate_min,
                cfg.release_maturity_gate_max,
            ))
            mcye_soft_factor = 1.0
            mcye_mode = str(getattr(cfg, "mcye_source_mode", "diagnostic")).strip().lower()
            if mcye_mode == "soft_constraint" and np.isfinite(raw_mcye_effective_fraction):
                ratio = raw_mcye_effective_fraction / max(cfg.mcye_soft_ref_effective_fraction, EPS)
                mcye_soft_factor = float(np.clip(ratio, cfg.mcye_soft_min_factor, cfg.mcye_soft_max_factor))
            elif mcye_mode not in {"off", "diagnostic", "soft_constraint"}:
                raise ValueError(
                    f"Unknown mcye_source_mode={cfg.mcye_source_mode!r}; expected off/diagnostic/soft_constraint"
                )
            effective_toxic_biomass_fraction = float(np.clip(
                toxic_producer_prior_state * bloom_stage_gate * release_maturity_gate * mcye_soft_factor,
                cfg.toxic_producer_prior_min,
                cfg.toxic_pool_fraction_max,
            ))
            ftox_for_synthesis = effective_toxic_biomass_fraction
            ftox_for_pool = effective_toxic_biomass_fraction
            X_tox_pool_basis = max(0.0, X_micro * ftox_for_pool)
            B_tox_live_gDW_m3 = max(0.0, X_tox_pool_basis)
            B_tox_inactive_gDW_m3 = max(0.0, cfg.biomass_scale_gDW_m3_per_state * Bm_inactive * ftox_for_pool)
            G_tox_live_gDW_m3_h = max(0.0, cfg.biomass_scale_gDW_m3_per_state * micro_growth * ftox_for_pool)
            cap_realized_unit = max(0.0, cap_upper_unit * float(np.clip(cfg.capacity_use_fraction, 0.0, 1.0)))
            MCLR_capacity_ceiling_flux_mmol_m3_h = max(0.0, cap_realized_unit * B_tox_live_gDW_m3)

            other_potential_h = cfg.other_mu_max_h * other_light * nutrient_gate * temp_gate
            suitability_difference_state = S_other - S_micro
            suitability_relative_advantage_state = safe_div(
                suitability_difference_state,
                S_other + S_micro + EPS,
                default=0.0,
            )
            suitability_log_advantage_state = float(np.log((S_other + EPS) / (S_micro + EPS)))
            succession_pressure_mode = str(getattr(cfg, "succession_pressure_mode", "potential_h") or "potential_h").strip().lower()
            if succession_pressure_mode in {"relative_suitability", "suitability", "delta_s"}:
                suitability_metric = str(getattr(cfg, "succession_suitability_metric", "raw_delta") or "raw_delta").strip().lower()
                if suitability_metric in {"raw_delta", "delta", "delta_s"}:
                    succession_metric_state = suitability_difference_state
                elif suitability_metric in {"relative_advantage", "arel", "a_rel", "relative"}:
                    succession_metric_state = suitability_relative_advantage_state
                elif suitability_metric in {"log_advantage", "alog", "a_log", "log"}:
                    succession_metric_state = suitability_log_advantage_state
                else:
                    raise ValueError(f"Unknown succession_suitability_metric={suitability_metric!r}")
                trigger_function = str(getattr(cfg, "succession_trigger_function", "linear_clip") or "linear_clip").strip().lower()
                succession_excess_state = max(0.0, succession_metric_state - cfg.succession_suitability_delta_ref)
                if trigger_function in {"linear_clip", "linear", "legacy"}:
                    succession_signal = float(np.clip(
                        (succession_metric_state - cfg.succession_suitability_delta_ref) / max(cfg.succession_suitability_width, EPS),
                        0.0, 1.0,
                    ))
                elif trigger_function in {"positive_saturation", "positive_sat", "saturation", "michaelis"}:
                    succession_signal = sat(succession_excess_state, cfg.succession_suitability_width)
                else:
                    raise ValueError(f"Unknown succession_trigger_function={trigger_function!r}")
            elif succession_pressure_mode in {"potential_h", "current", "legacy"}:
                succession_metric_state = other_potential_h - micro_potential_h
                succession_excess_state = max(0.0, succession_metric_state - cfg.succession_advantage_ref_h)
                succession_signal = float(np.clip(
                    (succession_metric_state - cfg.succession_advantage_ref_h) / max(cfg.succession_advantage_width_h, EPS),
                    0.0, 1.0,
                ))
            else:
                raise ValueError(f"Unknown succession_pressure_mode={succession_pressure_mode!r}")
            succession_pressure_state += dth * (succession_signal - succession_pressure_state) / max(cfg.succession_tau_h, EPS)
            succession_pressure_state = float(np.clip(succession_pressure_state, 0.0, 1.0))

            micro_loss_rate = float(cfg.micro_base_loss_rate_h)
            other_loss_rate = float(cfg.other_loss_rate_h)
            micro_loss = max(0.0, micro_loss_rate * Bm)
            other_loss = max(0.0, other_loss_rate * Bo)
            dynamic_excess_total_flux = 0.0
            dynamic_excess_micro_flux = 0.0
            dynamic_excess_other_flux = 0.0
            if dynamic_capacity_enabled and dynamic_excess_exit_enabled:
                excess_tau_h = max(float(meta.get("dynamic_excess_tau_h", cfg.late_mu_memory_tau_h)), EPS)
                dynamic_excess_rate_h = (1.0 - dynamic_capacity_fraction) / excess_tau_h
                live_total = max(Bm + Bo, 0.0)
                excess_inventory = max(0.0, live_total - dynamic_capacity_state)
                dynamic_excess_total_flux = min(
                    dynamic_excess_rate_h * excess_inventory, live_total / max(dth, EPS)
                )
                dynamic_excess_micro_flux = dynamic_excess_total_flux * safe_div(Bm, live_total, default=0.0)
                dynamic_excess_other_flux = dynamic_excess_total_flux - dynamic_excess_micro_flux
                micro_loss += dynamic_excess_micro_flux
                other_loss += dynamic_excess_other_flux
            composition_turnover_rate_h = 0.0
            composition_turnover_flux = 0.0
            succession_effect_mode = str(getattr(cfg, "succession_effect_mode", "loss") or "loss").strip().lower()
            if succession_effect_mode in {"composition_turnover", "turnover"}:
                p_turn = float(np.clip(safe_div(Bm, Bm + Bo, default=0.0), 0.0, 1.0))
                composition_turnover_rate_h = cfg.composition_turnover_rate_h * succession_pressure_state
                composition_turnover_flux = max(0.0, composition_turnover_rate_h * p_turn * (1.0 - p_turn) * (Bm + Bo))
                composition_turnover_flux = min(composition_turnover_flux, Bm / max(dth, EPS))

            micro_inactivation_rate_h = 0.0
            other_inactivation_rate_h = 0.0
            micro_inactivation_flux = 0.0
            other_inactivation_flux = 0.0
            micro_inactive_lysis_flux = 0.0
            other_inactive_lysis_flux = 0.0
            micro_inactive_settling_flux = 0.0
            other_inactive_settling_flux = 0.0
            micro_inactive_exit = 0.0
            other_inactive_exit = 0.0
            if late_inactivation_enabled:
                micro_inactivation_rate_h = max(0.0, cfg.micro_inactivation_max_rate_h * late_inactivation_memory_state)
                other_inactivation_rate_h = max(0.0, cfg.other_inactivation_max_rate_h * late_inactivation_memory_state)
                micro_inactivation_flux = min(micro_inactivation_rate_h * Bm, Bm / max(dth, EPS))
                other_inactivation_flux = min(other_inactivation_rate_h * Bo, Bo / max(dth, EPS))
                micro_inactive_lysis_flux = max(0.0, cfg.micro_inactive_lysis_rate_h * Bm_inactive)
                micro_inactive_settling_flux = max(0.0, cfg.micro_inactive_settling_rate_h * Bm_inactive)
                micro_inactive_exit = micro_inactive_lysis_flux + micro_inactive_settling_flux
                if dth * micro_inactive_exit > Bm_inactive and micro_inactive_exit > 0.0:
                    scale = Bm_inactive / max(dth * micro_inactive_exit, EPS)
                    micro_inactive_lysis_flux *= scale
                    micro_inactive_settling_flux *= scale
                other_inactive_lysis_flux = max(0.0, cfg.other_inactive_lysis_rate_h * Bo_inactive)
                other_inactive_settling_flux = max(0.0, cfg.other_inactive_settling_rate_h * Bo_inactive)
                other_inactive_exit = other_inactive_lysis_flux + other_inactive_settling_flux
                if dth * other_inactive_exit > Bo_inactive and other_inactive_exit > 0.0:
                    scale = Bo_inactive / max(dth * other_inactive_exit, EPS)
                    other_inactive_lysis_flux *= scale
                    other_inactive_settling_flux *= scale
            inactive_module_residual_flux = 0.0

            Bm_raw = Bm + dth * (micro_growth - micro_loss - composition_turnover_flux - micro_inactivation_flux)
            Bo_raw = Bo + dth * (other_growth - other_loss + composition_turnover_flux - other_inactivation_flux)
            Bm_inactive_raw = Bm_inactive + dth * (micro_inactivation_flux - micro_inactive_lysis_flux - micro_inactive_settling_flux)
            Bo_inactive_raw = Bo_inactive + dth * (other_inactivation_flux - other_inactive_lysis_flux - other_inactive_settling_flux)
            neg_corr = max(0.0, -Bm_raw) + max(0.0, -Bo_raw) + max(0.0, -Bm_inactive_raw) + max(0.0, -Bo_inactive_raw)
            Bm_next = max(0.0, Bm_raw)
            Bo_next = max(0.0, Bo_raw)
            Bm_inactive_next = max(0.0, Bm_inactive_raw)
            Bo_inactive_next = max(0.0, Bo_inactive_raw)
            expected_next = (
                Bm + Bo + Bm_inactive + Bo_inactive
                + dth * (
                    micro_growth + other_growth
                    - micro_loss - other_loss
                    - micro_inactive_lysis_flux - micro_inactive_settling_flux
                    - other_inactive_lysis_flux - other_inactive_settling_flux
                )
            )
            residual = (Bm_next + Bo_next + Bm_inactive_next + Bo_inactive_next) - expected_next - neg_corr

            if bool(cfg.enable_aggregation_state):
                stability = 1.0 / (1.0 + max(wind, 0.0) / max(cfg.aggregation_wind_half_mps, EPS))
                biomass_patch = sigmoid(Bm / max(cfg.carrying_capacity_state, EPS), 0.35, 0.16) * dominance * stability
                release_patch = release_pressure_state * (0.5 + 0.5 * maturity_state) * stability
                memory_patch = scum_accumulation_memory_state * (0.5 + 0.5 * bloom_patch_memory_state)
                agg_target = float(np.clip(0.50 * biomass_patch + 0.25 * release_patch + 0.25 * memory_patch, 0.0, 1.0))
                aggregation_state += dth * (agg_target - aggregation_state) / max(cfg.aggregation_tau_h, EPS)
                aggregation_state = float(np.clip(aggregation_state, 0.0, 1.0))
            else:
                aggregation_state = 0.0
            surface_release_gate = float(np.clip((1.0 - cfg.surface_sample_release_weight) + cfg.surface_sample_release_weight * release_pressure_state, 0.0, 1.0))
            surface_maturity_gate = float(np.clip((1.0 - cfg.surface_sample_maturity_weight) + cfg.surface_sample_maturity_weight * maturity_state, 0.0, 1.0))
            patch_memory_surface_contribution = cfg.patch_memory_surface_extra * scum_accumulation_memory_state * surface_release_gate * surface_maturity_gate
            surface_expression_factor = float(np.clip(
                1.0 + cfg.aggregation_surface_amplification * aggregation_state
                + cfg.surface_sample_enrichment_max_extra * aggregation_state * surface_release_gate * surface_maturity_gate
                + patch_memory_surface_contribution,
                1.0, cfg.surface_sample_factor_max,
            ))

            if source_state_scenario is not None and scenario == source_state_scenario:
                source_rows.append({
                    "time": gdf["time"].iloc[k],
                    "grid_id": str(gdf["grid_id"].iloc[k]),
                    "lon": lon,
                    "lat": lat,
                    "depth_m": depth,
                    "wind_speed_mps": wind,
                    "B_micro_state": Bm,
                    "B_other_state": Bo,
                    "B_total_state": Bm + Bo,
                    "B_micro_inactive_state": Bm_inactive,
                    "B_other_inactive_state": Bo_inactive,
                    "B_live_total_state": Bm + Bo,
                    "B_watercolumn_total_with_inactive_state": Bm + Bo + Bm_inactive + Bo_inactive,
                    "p_micro_live_state": safe_div(Bm, Bm + Bo, default=0.0),
                    "active_state": active_state,
                    "allocation_state": allocation_state,
                    "release_pressure_state": release_pressure_state,
                    "growth_competence_state": growth_competence_state,
                    "growth_competence_eq_state": growth_competence_eq_state,
                    "growth_competence_enabled": int(rgc_enabled),
                    "surface_expression_factor": surface_expression_factor,
                    "aggregation_state": aggregation_state,
                    "f_toxigenic_state": ftox,
                    "f_toxigenic_for_synthesis_state": ftox_for_synthesis,
                    "f_toxigenic_effective_pool_state": ftox_for_pool,
                    "raw_mcye_effective_fraction_state": raw_mcye_effective_fraction,
                    "toxic_producer_prior_state": toxic_producer_prior_state,
                    "bloom_stage_gate_state": bloom_stage_gate,
                    "release_maturity_gate_state": release_maturity_gate,
                    "mcye_soft_factor_state": mcye_soft_factor,
                    "effective_toxic_biomass_fraction_state": effective_toxic_biomass_fraction,
                    "B_tox_live_gDW_m3": B_tox_live_gDW_m3,
                    "B_tox_inactive_gDW_m3": B_tox_inactive_gDW_m3,
                    "G_tox_live_gDW_m3_h": G_tox_live_gDW_m3_h,
                    "F_tox_biomass_L_to_I_gDW_m3_h": max(0.0, cfg.biomass_scale_gDW_m3_per_state * micro_inactivation_flux * ftox_for_pool),
                    "F_tox_biomass_turnover_gDW_m3_h": max(0.0, cfg.biomass_scale_gDW_m3_per_state * composition_turnover_flux * ftox_for_pool),
                    "F_tox_biomass_clear_gDW_m3_h": max(0.0, cfg.biomass_scale_gDW_m3_per_state * micro_inactive_exit * ftox_for_pool),
                    "MCLR_capacity_ceiling_flux_mmol_m3_h": MCLR_capacity_ceiling_flux_mmol_m3_h,
                    "toxin_capacity_upper_per_gDW_h": cap_upper_unit,
                    "toxin_capacity_realized_per_gDW_h": cap_realized_unit,
                    "micro_growth_state_h": micro_growth,
                    "other_growth_state_h": other_growth,
                    "micro_loss_state_h": micro_loss,
                    "other_loss_state_h": other_loss,
                    "dynamic_capacity_enabled": int(dynamic_capacity_enabled),
                    "dynamic_k_feedback_enabled": int(dynamic_k_feedback_enabled),
                    "replacement_growth_enabled": int(replacement_growth_enabled),
                    "dynamic_excess_exit_enabled": int(dynamic_excess_exit_enabled),
                    "density_feedback_enabled": int(density_feedback_enabled),
                    "allocation_selectivity_enabled": int(allocation_selectivity_enabled),
                    "late_inactivation_enabled": int(late_inactivation_enabled),
                    "dynamic_capacity_target_fraction": dynamic_capacity_target_fraction,
                    "dynamic_capacity_ratio": dynamic_capacity_ratio,
                    "dynamic_capacity_fraction": dynamic_capacity_fraction,
                    "dynamic_capacity_state": dynamic_capacity_state,
                    "dynamic_micro_replacement_flux_state_h": micro_replacement_flux,
                    "dynamic_other_replacement_flux_state_h": other_replacement_flux,
                    "dynamic_excess_rate_h": dynamic_excess_rate_h,
                    "dynamic_excess_total_flux_state_h": dynamic_excess_total_flux,
                    "dynamic_excess_micro_flux_state_h": dynamic_excess_micro_flux,
                    "dynamic_excess_other_flux_state_h": dynamic_excess_other_flux,
                    "composition_turnover_flux_state_h": composition_turnover_flux,
                    "micro_inactivation_flux_state_h": micro_inactivation_flux,
                    "micro_inactive_lysis_flux_state_h": micro_inactive_lysis_flux,
                    "micro_inactive_settling_flux_state_h": micro_inactive_settling_flux,
                    "late_H_state": late_inactivation_memory_state,
                    "late_damage_burden_state": late_damage_burden_state,
                    "other_capacity_mode": other_capacity_mode,
                    "other_capacity_env_power": other_capacity_env_power,
                    "other_capacity_env_multiplier": other_capacity_env_multiplier,
                    "other_capacity_flux_state_h": other_capacity_flux,
                    "allocation_mode": mode,
                    "allocation_residual_state_h": allocation_residual,
                    "state_balance_residual": residual,
                    "negative_correction_mass_state": neg_corr,
                })

            scopes = ([("lake", "all"), ("basin", basin)] + [("station_nearest_grid", s) for s in station_names]) if collect_accums else []
            for scope, scope_id in scopes:
                key = (scenario, scope, scope_id, date, phase)
                if key not in accums:
                    accums[key] = new_daily_accum()
                    accums[key]["n_grid_days"] = 1.0
                update_daily_accum(
                    accums[key],
                    Bm=Bm, Bo=Bo, Bm_inactive=Bm_inactive, Bo_inactive=Bo_inactive,
                    Bm_next=Bm_next, Bo_next=Bo_next,
                    Bm_inactive_next=Bm_inactive_next, Bo_inactive_next=Bo_inactive_next,
                    mu_max=mu_max, mu_comm=mu_comm, density_micro=density_micro,
                    active=active_state, rho_gate=safe_div(mu_comm, mu_max),
                    mu_micro=mu_micro, other_mu_max=cfg.other_mu_max_h,
                    temp_gate=temp_gate, light_gate=other_light, n_gate=n_gate,
                    p_gate=p_gate, nutrient_gate=nutrient_gate,
                    env_product=env_product, env_single=env_single,
                    density_other=density_other, mu_other=mu_other,
                    phi_micro=phi_micro, phi_other=phi_other,
                    phi_ratio_micro_over_other=phi_ratio_micro_over_other,
                    ln_phi_micro_over_other=ln_phi_micro_over_other,
                    w_micro_allocation=w_micro_allocation, w_micro_phi_pure=w_phi_pure,
                    w_minus_p_micro=w_minus_p_micro,
                    logit_w_minus_logit_p=logit_w_minus_logit_p,
                    L_B_other_over_micro=L_B_other_over_micro,
                    L_phi_other_over_micro=L_phi_other_over_micro,
                    L_w_other_over_micro_gamma=L_w_other_over_micro_gamma,
                    logit_w_micro_gamma=logit_w_micro_gamma,
                    succession_signal=succession_signal, succession_pressure=succession_pressure_state,
                    suitability_relative_advantage=suitability_relative_advantage_state,
                    GT_ref=GT_ref,
                    micro_growth_ref=micro_growth_ref, other_growth_ref=other_growth_ref,
                    allocation_residual=allocation_residual,
                    micro_capacity_flux=micro_capacity_flux,
                    other_capacity_flux=other_capacity_flux,
                    other_capacity_basis=other_capacity_basis,
                    other_capacity_env_multiplier=other_capacity_env_multiplier,
                    other_recruitable_capacity_basis=other_recruitable_capacity_basis,
                    micro_capacity_slack_flux=micro_capacity_slack_flux,
                    other_capacity_slack_flux=other_capacity_slack_flux,
                    projection_distance_flux=projection_distance_flux,
                    unallocatable_budget_flux=unallocatable_budget_flux,
                    target_feasible=target_feasible,
                    capacity_budget_inconsistent=capacity_budget_inconsistent,
                    micro_capacity_binding=micro_capacity_binding,
                    other_capacity_binding=other_capacity_binding,
                    micro_growth=micro_growth, other_growth=other_growth,
                    micro_loss=micro_loss, other_loss=other_loss,
                    composition_turnover_flux=composition_turnover_flux,
                    micro_inactivation_flux=micro_inactivation_flux,
                    other_inactivation_flux=other_inactivation_flux,
                    micro_inactive_lysis_flux=micro_inactive_lysis_flux,
                    other_inactive_lysis_flux=other_inactive_lysis_flux,
                    micro_inactive_settling_flux=micro_inactive_settling_flux,
                    other_inactive_settling_flux=other_inactive_settling_flux,
                    inactive_module_residual_flux=inactive_module_residual_flux,
                    late_memory=late_inactivation_memory_state,
                    late_damage=late_damage_burden_state,
                    late_raw_pressure=late_raw_viability_pressure_state,
                    late_decline_signal=late_decline_signal,
                    late_saturation_gate=late_saturation_gate_state,
                    late_exhaustion_gate=late_exhaustion_gate_state,
                    micro_inactivation_rate_h=micro_inactivation_rate_h,
                    surface_factor=surface_expression_factor,
                    neg_corr=neg_corr, residual=residual,
                )

            st.update({
                "Bm": Bm_next,
                "Bo": Bo_next,
                "Bm_inactive": Bm_inactive_next,
                "Bo_inactive": Bo_inactive_next,
                "ftox": ftox,
                "ftox_gene_proxy_raw": ftox_gene_proxy_raw if np.isfinite(ftox_gene_proxy_raw) else np.nan,
                "active_state": active_state,
                "allocation_state": allocation_state,
                "release_pressure_state": release_pressure_state,
                "low_light_memory": low_light_memory,
                "aggregation_state": aggregation_state,
                "high_biomass_memory_state": high_biomass_memory_state,
                "dominance_memory_state": dominance_memory_state,
                "bloom_patch_memory_state": bloom_patch_memory_state,
                "scum_accumulation_memory_state": scum_accumulation_memory_state,
                "succession_pressure_state": succession_pressure_state,
                "late_inactivation_memory_state": late_inactivation_memory_state,
                "late_viability_burden_h": late_viability_burden_h,
                "late_damage_burden_state": late_damage_burden_state,
                "mu_capacity_memory_h": mu_capacity_memory_h,
                "dynamic_capacity_memory_fraction": dynamic_capacity_memory_fraction,
                "growth_competence_state": growth_competence_state,
            })

        if source_state_scenario is None:
            raise ValueError("coupled generator requires source_state_scenario")
        coupled_state = states[source_state_scenario]
        payload = {
            "time": gdf["time"].iloc[k],
            "grid_id": str(gdf["grid_id"].iloc[k]),
            "Bm_next": float(coupled_state["Bm"]),
            "Bo_next": float(coupled_state["Bo"]),
            "Bm_inactive_next": float(coupled_state.get("Bm_inactive", 0.0)),
            "Bo_inactive_next": float(coupled_state.get("Bo_inactive", 0.0)),
            "source_row": source_rows[-1],
        }
        incoming = yield payload
        if incoming is not None:
            for key in ["Bm", "Bo"]:
                if key in incoming:
                    coupled_state[key] = max(float(incoming[key]), 0.0)
        source_rows.clear()

    return rows_from_accums(accums), source_rows
