"""Published response-table interface and ecological target equations.

This module reads saved environmental response values. It contains no metabolic
network, optimization, network reconstruction, or complete lake simulation.
See the accompanying mathematical specification for update order and scope.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

COORDINATES = ["temp_C", "light_Wm2", "DIN_mgN_L", "DIP_mgP_L", "rho_growth_commitment"]
OUTPUTS = ["mu_max_h", "mu_committed_h", "toxin_capacity_upper_per_gDW_h"]
EPS = 1e-30


def clip(value, lower, upper):
    return float(np.clip(value, lower, upper))


def lag_step(state, target, time_constant_h, dt_h=1.0, lower=0.0, upper=1.0):
    """Bounded forward Euler update; state and target have the same units."""
    return clip(state + dt_h * (target - state) / max(time_constant_h, EPS), lower, upper)


def saturation(value, half_saturation):
    value = max(float(value), 0.0)
    return value / max(value + half_saturation, EPS)


def sigmoid(value, center, width):
    return float(1.0 / (1.0 + np.exp(-(value - center) / width)))


def temperature_factor(temp_C):
    """Empirical Gaussian, optimum 26 C, width 7.5 C."""
    return float(np.exp(-0.5 * ((float(temp_C) - 26.0) / 7.5) ** 2))


class SavedResponseQuery:
    """Exact-input-key lookup of saved capacity summaries, not new optimization."""

    def __init__(self, response_csv):
        raw = pd.read_csv(response_csv, low_memory=False)
        self.frame = raw.groupby(COORDINATES, as_index=False)[OUTPUTS].max()
        self.scales = self.frame[COORDINATES].std(ddof=0).replace(0, 1.0).fillna(1.0).to_numpy(float)
        self.scales[-1] = max(self.scales[-1], 0.10)
        self.coordinates = self.frame[COORDINATES].to_numpy(float)
        self.values = self.frame[OUTPUTS].to_numpy(float)
        self.tree = cKDTree(self.coordinates / self.scales)
        self.cache = {}

    def query(self, temp_C, light_Wm2, DIN_mgN_L, DIP_mgP_L, rho):
        light = max(float(light_Wm2), 0.0)
        if light <= 1e-9:
            return dict(zip(OUTPUTS, [0.0, 0.0, 0.0]))
        x = np.array([float(temp_C), light, max(float(DIN_mgN_L), 0.0),
                      max(float(DIP_mgP_L), 0.0), clip(rho, 0.0, 1.0)])
        key = tuple(x)
        if key in self.cache:
            return self.cache[key].copy()
        distance, index = self.tree.query(x / self.scales, k=min(4, len(self.frame)))
        distance = np.atleast_1d(distance).astype(float)
        index = np.atleast_1d(index).astype(int)
        daylight = self.coordinates[index, 1] > 1e-9
        if daylight.any():
            index, distance = index[daylight], distance[daylight]
        weight = 1.0 / (distance + 1e-9)
        weight /= weight.sum()
        values = np.maximum(np.average(self.values[index], axis=0, weights=weight), 0.0)
        out = dict(zip(OUTPUTS, map(float, values)))
        self.cache[key] = out
        return out.copy()


def allocation_target(B_micro, B_other, B_micro_inactive, B_other_inactive,
                      DIN_mgN_L, DIP_mgP_L, low_light_history):
    """Reference-scale maturity and allocation target; biomass in state units."""
    total = max(B_micro + B_other, EPS)
    fraction = B_micro / total
    density = (B_micro + B_micro_inactive + 0.8 * (B_other + B_other_inactive)) / 0.30
    maturity = clip(0.60 * sigmoid(density, 0.45, 0.18)
                    + 0.25 * sigmoid(fraction, 0.55, 0.15)
                    + 0.15 * sigmoid(B_micro / 0.30, 0.55, 0.20), 0.0, 1.0)
    nutrient_stress = 1.0 - min(saturation(DIN_mgN_L, 0.35), saturation(DIP_mgP_L, 0.030))
    target = clip(0.35 * nutrient_stress + 0.35 * maturity
                  + 0.15 * (1.0 - fraction) + 0.15 * low_light_history, 0.0, 1.0)
    return {"maturity": maturity, "nutrient_stress": nutrient_stress,
            "allocation_target": target}


def growth_commitment(allocation_history):
    return 0.80 - 0.55 * clip(allocation_history, 0.0, 1.0)


def density_capacity_step(memory_fraction, mu_max_h, light_Wm2, dt_h=1.0):
    """Daylight 72 h capacity memory; hold unchanged below 20 W m-2."""
    state = clip(memory_fraction, 0.0, 1.0)
    if light_Wm2 >= 20.0:
        target = clip(mu_max_h / 0.029147432795074, 0.0, 1.0)
        state = lag_step(state, target, 72.0, dt_h)
    return {"memory_fraction": state, "effective_density_scale": 0.30 * (0.5 + 0.5 * state)}


def empirical_growth_capacity(temp_C, light_Wm2, DIN_mgN_L, DIP_mgP_L, rho):
    """Frozen capacity comparator: inherited scale, no fitting to observations."""
    if light_Wm2 <= 1e-9:
        return {"mu_max_h": 0.0, "mu_committed_h": 0.0}
    capacity = (12.233206914782599 * 0.029147432795074 * temperature_factor(temp_C)
                * saturation(light_Wm2, 80.0)
                * min(saturation(DIN_mgN_L, 0.35), saturation(DIP_mgP_L, 0.030)))
    return {"mu_max_h": capacity, "mu_committed_h": clip(rho, 0.0, 1.0) * capacity}
