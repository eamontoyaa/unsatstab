"""Infinite-slope stability calculations."""

from __future__ import annotations

import numpy as np

from .parameters import HydraulicProperties, StrengthProperties


def bishop_chi(theta: np.ndarray, *, theta_r: float, theta_s: float) -> np.ndarray:
    """Bishop/Vanapalli-style effective-stress parameter."""
    return np.clip((theta - theta_r) / (theta_s - theta_r), 0.0, 1.0)


def factor_of_safety_profile(
    z: np.ndarray,
    pressure_head: np.ndarray,
    theta: np.ndarray,
    *,
    hydro: HydraulicProperties,
    strength: StrengthProperties,
) -> np.ndarray:
    """Compute infinite-slope factor of safety along a depth profile.

    FS = tan(phi)/tan(delta)
         + [c' - chi * gamma_w * psi * tan(phi)]
           / [gamma * z * sin(delta) * cos(delta)]

    Adapted from Baum et al. (2010), Eq. (18), with the unsaturated
    effective-stress correction in section 6, paragraph 42; also Baum et al.
    (2008), Eq. (18). Constant soil unit weight and numerical clipping are
    package choices.

    The sign convention is:
    - positive pressure head reduces stability;
    - negative pressure head represents suction and increases stability through
      the chi-scaled term.
    """

    z = np.asarray(z, dtype=float)
    psi = np.asarray(pressure_head, dtype=float)
    theta = np.asarray(theta, dtype=float)

    chi = np.where(
        psi < 0.0,
        bishop_chi(theta, theta_r=hydro.theta_r, theta_s=hydro.theta_s),
        1.0,
    )

    if abs(np.tan(hydro.slope_angle_rad)) <= 1.0e-12:
        return np.full_like(z, strength.max_factor_of_safety, dtype=float)

    denominator = (
        strength.unit_weight_soil
        * np.maximum(z, 1.0e-6)
        * np.sin(hydro.slope_angle_rad)
        * np.cos(hydro.slope_angle_rad)
    )

    friction = np.tan(strength.friction_angle_rad) / np.tan(hydro.slope_angle_rad)
    pore_pressure = (
        -chi * strength.unit_weight_water * psi * np.tan(strength.friction_angle_rad)
        / denominator
    )
    cohesion = strength.cohesion / denominator

    # The frictional strength contribution cannot be less than zero in this formulation.
    fs = np.where(friction + pore_pressure < 0.0, cohesion, friction + pore_pressure + cohesion)
    fs = np.where(z <= 1.0e-2, strength.max_factor_of_safety, fs)

    return np.clip(fs, 0.0, strength.max_factor_of_safety)
