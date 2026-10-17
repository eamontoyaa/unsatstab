"""Water-table rise approximations based on mass balance."""

from __future__ import annotations

import numpy as np


def cumulative_simpson_midpoint(
    y: np.ndarray,
    y_midpoint: np.ndarray,
    x: np.ndarray,
) -> np.ndarray:
    """Integrate values using Simpson's rule and explicit interval midpoints."""

    y = np.asarray(y, dtype=float)
    y_midpoint = np.asarray(y_midpoint, dtype=float)
    x = np.asarray(x, dtype=float)
    if y.shape != x.shape or y_midpoint.shape != (max(len(x) - 1, 0),):
        raise ValueError("Expected len(y) == len(x) and len(y_midpoint) == len(x) - 1.")

    out = np.zeros_like(y)
    for i in range(1, len(y)):
        dt = x[i] - x[i - 1]
        out[i] = out[i - 1] + dt * (y[i - 1] + 4.0 * y_midpoint[i - 1] + y[i]) / 6.0
    return out


def fillable_volume_curve(
    z: np.ndarray,
    theta: np.ndarray,
    *,
    theta_s: float,
) -> np.ndarray:
    """Compute fillable pore volume from each depth to the column base.

    V_f(z) = integral_z^z_w0 [theta_s - theta(xi)] dxi
    """

    z = np.asarray(z, dtype=float)
    theta = np.asarray(theta, dtype=float)

    if z.shape != theta.shape:
        raise ValueError("z and theta must have the same shape.")

    fillable = np.maximum(theta_s - theta, 0.0)
    volume = np.zeros_like(z)

    for i in range(len(z) - 2, -1, -1):
        dz = z[i + 1] - z[i]
        volume[i] = volume[i + 1] + 0.5 * (fillable[i] + fillable[i + 1]) * dz

    return volume


def water_table_rise_from_volume(
    available_volume: float,
    z: np.ndarray,
    fillable_volume: np.ndarray,
    *,
    initial_water_table_depth: float,
) -> float:
    """Estimate water-table rise by inverting fillable pore volume."""

    if available_volume <= 0.0:
        return 0.0

    if available_volume >= fillable_volume[0]:
        return float(initial_water_table_depth)

    for i in range(len(z) - 1):
        high = fillable_volume[i]
        low = fillable_volume[i + 1]
        if high >= available_volume >= low:
            denom = high - low
            if abs(denom) < 1.0e-15:
                z_w = z[i]
            else:
                fraction = (high - available_volume) / denom
                z_w = z[i] + fraction * (z[i + 1] - z[i])
            return float(max(initial_water_table_depth - z_w, 0.0))

    return 0.0
