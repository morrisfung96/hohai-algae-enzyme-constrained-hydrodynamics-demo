"""Frozen conservative offline graph transport and mobile-pool losses.

These functions operate on caller-supplied geometry and retained velocities.
They do not reconstruct solver-native historical FVCOM volume fluxes.
"""
from __future__ import annotations
import math
from typing import Any
import numpy as np

def _mass_flux_many_fields(
    mass: np.ndarray,
    concentration: np.ndarray,
    u: np.ndarray,
    v: np.ndarray,
    edges: dict[str, np.ndarray],
    depth: np.ndarray,
    face_width: np.ndarray,
    dt_s: float,
    transport_scale: float,
    max_total_out_fraction: float,
    max_edge_out_fraction: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mass = np.maximum(np.where(np.isfinite(mass), mass, 0.0), 0.0)
    concentration = np.maximum(np.where(np.isfinite(concentration), concentration, 0.0), 0.0)
    n = mass.shape[1]
    src = edges["src"]
    dst = edges["dst"]
    ex = edges["ex"]
    ey = edges["ey"]
    u = np.where(np.isfinite(u), u, 0.0)
    v = np.where(np.isfinite(v), v, 0.0)
    proj = u[src] * ex + v[src] * ey
    active = proj > 0
    mean_depth = 0.5 * (depth[src] + depth[dst])
    volume_flux = np.zeros_like(proj, dtype=float)
    volume_flux[active] = proj[active] * face_width[active] * mean_depth[active] * dt_s * float(transport_scale)
    raw = concentration[:, src] * volume_flux[None, :]
    raw = np.maximum(np.where(np.isfinite(raw), raw, 0.0), 0.0)
    raw = np.minimum(raw, float(max_edge_out_fraction) * mass[:, src])
    out_raw = np.vstack([np.bincount(src, weights=raw[i], minlength=n) for i in range(raw.shape[0])])
    cap = float(max_total_out_fraction) * mass
    scale = np.ones_like(mass)
    over = out_raw > cap
    scale[over] = cap[over] / np.maximum(out_raw[over], 1e-30)
    amount = raw * scale[:, src]
    out = np.vstack([np.bincount(src, weights=amount[i], minlength=n) for i in range(amount.shape[0])])
    inc = np.vstack([np.bincount(dst, weights=amount[i], minlength=n) for i in range(amount.shape[0])])
    return np.maximum(mass - out + inc, 0.0), out, inc

def _c4_transport_many(
    mass: np.ndarray,
    ti_flow: int,
    dt_h: float,
    geom: dict[str, Any],
    cfg: Any,
    pool_kind: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    volume = geom["volume"]
    dt_s = max(float(dt_h), 0.0) * 3600.0
    conc = mass / np.maximum(volume[None, :], 1e-30)
    if pool_kind == "dead":
        return _mass_flux_many_fields(
            mass, conc, geom["u_particle"][ti_flow, :], geom["v_particle"][ti_flow, :],
            geom["edges"], geom["depth"], geom["face_width"], dt_s,
            cfg.transport_scale_particle, cfg.max_total_out_fraction_particle, cfg.max_edge_out_fraction
        )
    return _mass_flux_many_fields(
        mass, conc, geom["u_dissolved"][ti_flow, :], geom["v_dissolved"][ti_flow, :],
        geom["edges"], geom["depth"], geom["face_width"], dt_s,
        cfg.transport_scale_dissolved, cfg.max_total_out_fraction_dissolved, cfg.max_edge_out_fraction
    )

def _c5_transport_loss_step(dead: np.ndarray, diss: np.ndarray, ti_flow: int, dt_h: float, geom: dict[str, Any], cfg: Any) -> tuple[np.ndarray, np.ndarray, dict[str, Any]]:
    dead_arr, dead_out, dead_in = _c4_transport_many(dead[None, :], ti_flow, dt_h, geom, cfg, "dead")
    diss_arr, diss_out, diss_in = _c4_transport_many(diss[None, :], ti_flow, dt_h, geom, cfg, "dissolved")
    dead_t = dead_arr[0]
    diss_t = diss_arr[0]
    dead_rate = min(float(cfg.particulate_settling_base_h), float(cfg.max_particulate_loss_fraction_per_h))
    diss_rate = min(float(cfg.dissolved_decay_base_h), float(cfg.max_dissolved_loss_fraction_per_h))
    dead_loss = dead_t * (1.0 - math.exp(-dead_rate * max(dt_h, 0.0)))
    diss_loss = diss_t * (1.0 - math.exp(-diss_rate * max(dt_h, 0.0)))
    dead_f = np.maximum(dead_t - dead_loss, 0.0)
    diss_f = np.maximum(diss_t - diss_loss, 0.0)
    diag = {
        "dead_transport_in": dead_in[0],
        "dead_transport_out": dead_out[0],
        "dissolved_transport_in": diss_in[0],
        "dissolved_transport_out": diss_out[0],
        "particulate_settling_loss": dead_loss,
        "dissolved_degradation_loss": diss_loss,
        "transport_closure": float(np.nansum(dead_in - dead_out + diss_in - diss_out)),
    }
    return dead_f, diss_f, diag

transport_loss_step = _c5_transport_loss_step
