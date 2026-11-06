"""Query the archived metabolic response table without an upstream solver.

This module retains the adopted response-query arithmetic and cache behavior.
The published entry point is the response table, not a metabolic-network build.
See data/response/README.txt for support, units and cache/replay requirements.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

import numpy as np
import pandas as pd

try:
    from scipy.spatial import cKDTree
except Exception:  # fallback retained from the archived implementation
    cKDTree = None


@dataclass
class ResponseConfig:
    """Only configuration fields consumed by the archived response class.

    response_k_env is retained as configuration metadata; lookup's k_env
    argument, whose archived default is 4, controls each individual query.
    """

    exact_night_light_threshold_Wm2: float = 1e-9
    response_k_env: int = 4
    require_native_v16_response: bool = False


def first_col(columns: Iterable[str], candidates: list[str]) -> Optional[str]:
    cols = list(columns)
    exact = {c: c for c in cols}
    lower = {str(c).strip().lower(): c for c in cols}
    for cand in candidates:
        if cand in exact:
            return exact[cand]
        key = str(cand).strip().lower()
        if key in lower:
            return lower[key]
    return None


def canonicalize(df: pd.DataFrame, mapping: dict[str, list[str]], *, required: list[str], source: str) -> pd.DataFrame:
    out = df.copy()
    ren: dict[str, str] = {}
    for canonical, candidates in mapping.items():
        if canonical in out.columns:
            continue
        c = first_col(out.columns, candidates)
        if c is not None:
            ren[c] = canonical
    if ren:
        out = out.rename(columns=ren)
    missing = [c for c in required if c not in out.columns]
    if missing:
        raise ValueError(f"{source} missing required columns {missing}. Available columns: {list(df.columns)[:140]}")
    return out


RESPONSE_COLUMNS = {
    "env_id": ["env_id", "base_env_id", "matched_env_id"],
    "temp_C": ["temp_C", "temperature_C", "T_C", "temperature"],
    "light_Wm2": ["light_Wm2", "shortwave_Wm2", "PAR_Wm2", "I_Wm2", "light"],
    "DIN_mgN_L": ["DIN_mgN_L", "DIN", "NO3_mgN_L", "nitrate_mgN_L"],
    "DIP_mgP_L": ["DIP_mgP_L", "DIP", "SRP_mgP_L", "PO4_mgP_L"],
    "rho_growth_commitment": ["rho_growth_commitment", "rho", "growth_commitment", "growth_retention", "epsilon"],
    "mu_max_h": ["mu_max_h", "mu_capacity_h", "mu_max", "growth_capacity_h"],
    "mu_committed_h": ["mu_committed_h", "mu_min_committed_h", "growth_committed_h"],
    "toxin_capacity_upper_per_gDW_h": ["toxin_capacity_upper_per_gDW_h", "toxin_capacity_mmol_gDW_h", "toxin_capacity_per_gDW_h", "toxin_flux_max_mmol_gDW_h"],
    "mu_toxin_solution_h": ["mu_toxin_solution_h", "mu_joint_toxin_solution_h", "growth_at_toxin_solution_h"],
    "growth_constraint_target_raw": ["growth_constraint_target_raw", "biomass_flux_target_raw"],
    "biomass_flux_toxin_solution_raw": ["biomass_flux_toxin_solution_raw", "biomass_flux_joint_solution_raw"],
    "growth_constraint_residual_raw": ["growth_constraint_residual_raw", "biomass_flux_constraint_residual_raw"],
    "growth_constraint_mode": ["growth_constraint_mode", "toxin_growth_constraint_mode"],
    "toxin_solution_status": ["toxin_solution_status", "solver_status_joint_toxin"],
    "response_mode": ["response_mode", "v16_response_mode"],
}


class ResponseSurfaceV16:
    """State-conditioned response-surface lookup with env + rho axis."""

    def __init__(self, path: Path, cfg: ResponseConfig):
        self.path = Path(path)
        self.cfg = cfg
        self.lookup_cache: dict[tuple, dict[str, float | str]] = {}
        self.max_cache_size = 300000

        raw = pd.read_csv(path, low_memory=False)
        df = canonicalize(raw, RESPONSE_COLUMNS, required=["mu_max_h", "toxin_capacity_upper_per_gDW_h"], source=str(path))
        if "env_id" not in df.columns:
            df["env_id"] = np.arange(len(df)).astype(str)
        if "rho_growth_commitment" not in df.columns:
            raise ValueError(
                "V16 response surface must contain rho_growth_commitment. "
                "For engineering-only smoke tests, build a compatibility response with 03_build_V16_response_surface_compat.py."
            )
        if "mu_committed_h" not in df.columns:
            df["mu_committed_h"] = pd.to_numeric(df["mu_max_h"], errors="coerce") * pd.to_numeric(df["rho_growth_commitment"], errors="coerce")
        if "response_mode" not in df.columns:
            df["response_mode"] = "unknown"
        for c in ["temp_C", "light_Wm2", "DIN_mgN_L", "DIP_mgP_L", "rho_growth_commitment"]:
            if c not in df.columns:
                df[c] = 0.0
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0.0)
        optional_nonnegative_numeric = ["mu_toxin_solution_h", "growth_constraint_target_raw", "biomass_flux_toxin_solution_raw"]
        optional_signed_numeric = ["growth_constraint_residual_raw"]
        optional_numeric = optional_nonnegative_numeric + optional_signed_numeric
        optional_text = ["growth_constraint_mode", "toxin_solution_status"]
        for c in optional_numeric:
            if c not in df.columns:
                df[c] = np.nan
        for c in optional_text:
            if c not in df.columns:
                df[c] = ""
            df[c] = df[c].astype(str)
        for c in ["mu_max_h", "mu_committed_h", "toxin_capacity_upper_per_gDW_h"] + optional_nonnegative_numeric:
            df[c] = pd.to_numeric(df[c], errors="coerce").clip(lower=0.0)
        for c in optional_signed_numeric:
            df[c] = pd.to_numeric(df[c], errors="coerce")
        # Collapse exact duplicates, keeping conservative max capacity for the same state-conditioned point.
        key_cols = ["temp_C", "light_Wm2", "DIN_mgN_L", "DIP_mgP_L", "rho_growth_commitment"]
        mode_text = ";".join(sorted(set(df["response_mode"].astype(str).unique())))[:200]
        df = df.groupby(key_cols, as_index=False).agg({
            "mu_max_h": "max",
            "mu_committed_h": "max",
            "toxin_capacity_upper_per_gDW_h": "max",
            "mu_toxin_solution_h": "max",
            "growth_constraint_target_raw": "max",
            "biomass_flux_toxin_solution_raw": "max",
            "growth_constraint_residual_raw": "max",
            "growth_constraint_mode": "first",
            "toxin_solution_status": "first",
        })
        df["env_id"] = [f"V16ENV_{i:05d}" for i in range(len(df))]
        df["response_mode"] = mode_text
        self.df = df.reset_index(drop=True)
        self.env_cols = key_cols
        scales = self.df[self.env_cols].std(ddof=0).replace(0, 1.0).fillna(1.0)
        # Make rho distance reasonably important even if numeric variance is small.
        if "rho_growth_commitment" in scales.index:
            scales.loc["rho_growth_commitment"] = max(float(scales.loc["rho_growth_commitment"]), 0.10)
        self.env_values = self.df[self.env_cols].to_numpy(float)
        self.env_scales = scales.to_numpy(float)
        self.env_ids = self.df["env_id"].astype(str).to_numpy()
        self.mu_max_values = self.df["mu_max_h"].to_numpy(float)
        self.mu_committed_values = self.df["mu_committed_h"].to_numpy(float)
        self.toxin_upper_values = self.df["toxin_capacity_upper_per_gDW_h"].to_numpy(float)
        self.mu_toxin_solution_values = self.df["mu_toxin_solution_h"].to_numpy(float)
        self.growth_constraint_target_raw_values = self.df["growth_constraint_target_raw"].to_numpy(float)
        self.biomass_flux_toxin_solution_raw_values = self.df["biomass_flux_toxin_solution_raw"].to_numpy(float)
        self.growth_constraint_residual_raw_values = self.df["growth_constraint_residual_raw"].to_numpy(float)
        self.growth_constraint_modes = self.df["growth_constraint_mode"].astype(str).to_numpy()
        self.toxin_solution_status_values = self.df["toxin_solution_status"].astype(str).to_numpy()
        self.response_mode = mode_text
        self.is_native = "native" in mode_text.lower() and "compat" not in mode_text.lower()
        if cfg.require_native_v16_response and not self.is_native:
            raise ValueError(f"Config requires native V16 response, but response_mode={mode_text!r} in {path}")
        # KDTree on standardized coordinates. Fallback to numpy argpartition.
        self.tree = None
        if cKDTree is not None and len(self.env_values):
            self.tree = cKDTree(self.env_values / self.env_scales)

    @staticmethod
    def _round_key(temp_C: float, light_Wm2: float, DIN_mgN_L: float, DIP_mgP_L: float, rho: float, k_env: int) -> tuple:
        return (
            round(float(temp_C), 3), round(float(light_Wm2), 3), round(float(DIN_mgN_L), 5),
            round(float(DIP_mgP_L), 6), round(float(rho), 3), int(k_env)
        )

    def lookup(self, *, temp_C: float, light_Wm2: float, DIN_mgN_L: float, DIP_mgP_L: float, rho_growth_commitment: float, k_env: int = 4) -> dict[str, float | str]:
        light = max(float(light_Wm2), 0.0)
        rho = float(np.clip(rho_growth_commitment, 0.0, 1.0))
        if light <= self.cfg.exact_night_light_threshold_Wm2:
            return {
                "matched_env_ids": "EXACT_NIGHT_GATE",
                "matched_env_distance": 0.0,
                "mu_max_h": 0.0,
                "mu_committed_h": 0.0,
                "toxin_capacity_upper_per_gDW_h": 0.0,
                "mu_toxin_solution_h": 0.0,
                "growth_constraint_target_raw": 0.0,
                "biomass_flux_toxin_solution_raw": 0.0,
                "growth_constraint_residual_raw": 0.0,
                "growth_constraint_mode": "night_zero",
                "toxin_solution_status": "night_zero",
                "response_mode": self.response_mode,
            }
        key = self._round_key(temp_C, light, DIN_mgN_L, DIP_mgP_L, rho, k_env)
        cached = self.lookup_cache.get(key)
        if cached is not None:
            return cached
        query = np.array([temp_C, light, DIN_mgN_L, DIP_mgP_L, rho], dtype=float)
        q_std = query / self.env_scales
        k = max(1, min(int(k_env), len(self.env_values)))
        if self.tree is not None:
            dist, idx = self.tree.query(q_std, k=k)
            idx = np.atleast_1d(idx).astype(int)
            dist = np.atleast_1d(dist).astype(float)
        else:
            dist_all = np.sqrt(np.nansum(((self.env_values - query) / self.env_scales) ** 2, axis=1))
            if k < len(dist_all):
                idx = np.argpartition(dist_all, k - 1)[:k]
                idx = idx[np.argsort(dist_all[idx])]
            else:
                idx = np.arange(len(dist_all))
            dist = dist_all[idx]
        # Avoid borrowing exact night rows for nonzero-light queries if possible.
        daytime = self.env_values[idx, 1] > self.cfg.exact_night_light_threshold_Wm2
        if daytime.any():
            idx = idx[daytime]
            dist = dist[daytime]
        w = 1.0 / (dist + 1e-9)
        w = w / w.sum()
        mu_max = float(np.average(self.mu_max_values[idx], weights=w))
        mu_comm = float(np.average(self.mu_committed_values[idx], weights=w))
        tox = float(np.average(self.toxin_upper_values[idx], weights=w))
        mu_toxin = float(np.average(self.mu_toxin_solution_values[idx], weights=w))
        growth_target_raw = float(np.average(self.growth_constraint_target_raw_values[idx], weights=w))
        biomass_toxin_raw = float(np.average(self.biomass_flux_toxin_solution_raw_values[idx], weights=w))
        growth_resid_raw = float(np.average(self.growth_constraint_residual_raw_values[idx], weights=w))
        mean_dist = float(np.average(dist, weights=w))
        growth_modes = ";".join(sorted(set(self.growth_constraint_modes[idx].astype(str))))[:200]
        toxin_status = ";".join(sorted(set(self.toxin_solution_status_values[idx].astype(str))))[:200]
        out = {
            "matched_env_ids": ";".join(self.env_ids[idx].astype(str)),
            "matched_env_distance": mean_dist,
            "mu_max_h": max(0.0, mu_max),
            "mu_committed_h": max(0.0, mu_comm),
            "toxin_capacity_upper_per_gDW_h": max(0.0, tox),
            "mu_toxin_solution_h": max(0.0, mu_toxin) if np.isfinite(mu_toxin) else float("nan"),
            "growth_constraint_target_raw": max(0.0, growth_target_raw) if np.isfinite(growth_target_raw) else float("nan"),
            "biomass_flux_toxin_solution_raw": max(0.0, biomass_toxin_raw) if np.isfinite(biomass_toxin_raw) else float("nan"),
            "growth_constraint_residual_raw": growth_resid_raw if np.isfinite(growth_resid_raw) else float("nan"),
            "growth_constraint_mode": growth_modes,
            "toxin_solution_status": toxin_status,
            "response_mode": self.response_mode,
        }
        if len(self.lookup_cache) < self.max_cache_size:
            self.lookup_cache[key] = out
        return out
