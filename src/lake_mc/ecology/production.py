"""Exact source extraction from scripts/v16_3_conservative_production_allocation.py; algorithm retained.
No FACHB905/SBML/LP imports. See ecology extraction provenance.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np



@dataclass(frozen=True)
class ProjectionResult:
    P_M_projected: float
    P_O_projected: float
    production_residual: float
    M_capacity_binding: bool
    O_capacity_binding: bool
    target_feasible: bool
    capacity_budget_inconsistent: bool
    projection_distance: float
    unallocatable_budget: float
    diagnostic_code: str
    lower_bound_M: float
    upper_bound_M: float
    P_M_target: float
    P_O_target: float
    M_capacity_slack: float
    O_capacity_slack: float

def _finite_nonnegative(name: str, value: float) -> float:
    try:
        x = float(value)
    except Exception as exc:  # pragma: no cover - defensive
        raise ValueError(f"{name} must be numeric") from exc
    if not np.isfinite(x):
        raise ValueError(f"{name} must be finite, got {value!r}")
    if x < 0.0:
        raise ValueError(f"{name} must be non-negative, got {value!r}")
    return x

def project_two_species_conservative_production(
    *,
    B_M: float,
    B_O: float,
    P_total_native: float,
    w_M_target: float,
    C_M: float,
    C_O: float,
    tolerance: float = 1e-12,
) -> ProjectionResult:
    """Project target production onto capacity bounds while conserving total.

    ``B_M`` and ``B_O`` are accepted for lineage/diagnostic completeness; the
    projection itself depends only on the total budget, target weight and
    capacity bounds.
    """
    _finite_nonnegative("B_M", B_M)
    _finite_nonnegative("B_O", B_O)
    P_total = _finite_nonnegative("P_total_native", P_total_native)
    C_M = _finite_nonnegative("C_M", C_M)
    C_O = _finite_nonnegative("C_O", C_O)
    tol = _finite_nonnegative("tolerance", tolerance)
    try:
        w = float(w_M_target)
    except Exception as exc:  # pragma: no cover - defensive
        raise ValueError("w_M_target must be numeric") from exc
    if not np.isfinite(w):
        raise ValueError(f"w_M_target must be finite, got {w_M_target!r}")
    w = float(np.clip(w, 0.0, 1.0))

    P_M_target = w * P_total
    P_O_target = P_total - P_M_target

    if P_total <= tol:
        return ProjectionResult(
            P_M_projected=0.0,
            P_O_projected=0.0,
            production_residual=0.0,
            M_capacity_binding=False,
            O_capacity_binding=False,
            target_feasible=True,
            capacity_budget_inconsistent=False,
            projection_distance=0.0,
            unallocatable_budget=0.0,
            diagnostic_code="zero_total_production",
            lower_bound_M=0.0,
            upper_bound_M=0.0,
            P_M_target=0.0,
            P_O_target=0.0,
            M_capacity_slack=C_M,
            O_capacity_slack=C_O,
        )

    lower = max(0.0, P_total - C_O)
    upper = min(C_M, P_total)
    capacity_budget_inconsistent = P_total > C_M + C_O + tol
    if capacity_budget_inconsistent or lower > upper + tol:
        # No feasible allocation can conserve the native total budget.  Return
        # the closest capacity-saturating allocation and expose the residual.
        P_M = min(C_M, P_total)
        P_O = min(C_O, max(0.0, P_total - P_M))
        residual = (P_M + P_O) - P_total
        return ProjectionResult(
            P_M_projected=float(P_M),
            P_O_projected=float(P_O),
            production_residual=float(residual),
            M_capacity_binding=bool(C_M - P_M <= tol),
            O_capacity_binding=bool(C_O - P_O <= tol),
            target_feasible=False,
            capacity_budget_inconsistent=True,
            projection_distance=float(np.hypot(P_M - P_M_target, P_O - P_O_target)),
            unallocatable_budget=float(max(0.0, P_total - C_M - C_O)),
            diagnostic_code="capacity_budget_inconsistent",
            lower_bound_M=float(lower),
            upper_bound_M=float(upper),
            P_M_target=float(P_M_target),
            P_O_target=float(P_O_target),
            M_capacity_slack=float(max(0.0, C_M - P_M)),
            O_capacity_slack=float(max(0.0, C_O - P_O)),
        )

    P_M = float(np.clip(P_M_target, lower, upper))
    P_O = float(P_total - P_M)
    residual = (P_M + P_O) - P_total
    target_feasible = (
        P_M_target >= lower - tol
        and P_M_target <= upper + tol
        and P_M_target <= C_M + tol
        and P_O_target <= C_O + tol
    )
    if abs(residual) > tol:
        diagnostic = "projection_residual_exceeds_tolerance"
    elif target_feasible:
        diagnostic = "target_feasible"
    elif abs(P_M - lower) <= tol:
        diagnostic = "O_capacity_binding_target_projected"
    elif abs(P_M - upper) <= tol:
        diagnostic = "M_capacity_binding_target_projected"
    else:
        diagnostic = "target_projected_inside_bounds"

    return ProjectionResult(
        P_M_projected=P_M,
        P_O_projected=P_O,
        production_residual=float(residual),
        M_capacity_binding=bool(C_M - P_M <= tol and P_M > tol),
        O_capacity_binding=bool(C_O - P_O <= tol and P_O > tol),
        target_feasible=bool(target_feasible),
        capacity_budget_inconsistent=False,
        projection_distance=float(np.hypot(P_M - P_M_target, P_O - P_O_target)),
        unallocatable_budget=0.0,
        diagnostic_code=diagnostic,
        lower_bound_M=float(lower),
        upper_bound_M=float(upper),
        P_M_target=float(P_M_target),
        P_O_target=float(P_O_target),
        M_capacity_slack=float(max(0.0, C_M - P_M)),
        O_capacity_slack=float(max(0.0, C_O - P_O)),
    )
