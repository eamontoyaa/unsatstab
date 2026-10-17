"""Analytical unsaturated infiltration response.

This module provides the analytical unsaturated-zone response associated with
Gardner and Srivastava–Yeh formulations.
"""

from __future__ import annotations

from dataclasses import replace
from functools import lru_cache
import numpy as np
from scipy.integrate import quad
from scipy.special import erfc, erfcx

from .parameters import HydraulicProperties
from .rainfall import RainfallSeries
from .roots import compute_roots


def gardner_water_content(
    pressure_head: np.ndarray | float,
    *,
    alpha: float,
    theta_s: float,
    theta_r: float,
) -> np.ndarray:
    """Return the intrinsic Gardner retention curve.

    This is the soil-water characteristic relation used for Figure 4 in Baum
    et al. (2010), Eq. (4e), with psi_0 = 0; also Baum et al. (2008),
    Eq. (4C): ``theta = theta_r + (theta_s - theta_r) exp(alpha psi)``
    for non-positive pressure head.  It deliberately uses the intrinsic
    ``alpha`` (not the slope-transformed value used by the 1D flow solution).
    """

    values = (alpha, theta_s, theta_r)
    if not all(np.isfinite(value) for value in values):
        raise ValueError("Gardner parameters must be finite.")
    if alpha <= 0.0 or not 0.0 <= theta_r < theta_s <= 1.0:
        raise ValueError("Expected alpha > 0 and 0 <= theta_r < theta_s <= 1.")

    psi = np.asarray(pressure_head, dtype=float)
    theta = theta_r + (theta_s - theta_r) * np.exp(alpha * np.minimum(psi, 0.0))
    return np.clip(theta, theta_r, theta_s)


def initial_conductivity(
    z: np.ndarray | float,
    *,
    hydro: HydraulicProperties,
    depth: float | None = None,
) -> np.ndarray | float:
    """Initial hydraulic-conductivity profile.

    K_0(z) = I_ZLT - (I_ZLT - K_s) * exp[-a * (z_w0 - z)]

    Baum et al. (2010), Eq. (5c), corrected by Baum and Godt (2013),
    with psi_0 = 0 and renamed depth/slope parameters.
    """

    z_w0 = hydro.unsaturated_domain_depth if depth is None else depth
    return hydro.long_term_infiltration - (
        hydro.long_term_infiltration - hydro.saturated_hydraulic_conductivity
    ) * np.exp(-hydro.transformed_alpha * (z_w0 - z))


def pressure_from_conductivity(
    conductivity: np.ndarray,
    *,
    hydro: HydraulicProperties,
) -> np.ndarray:
    """Convert conductivity to pressure head with psi_0 = 0.

    The TRIGRS set uses Baum et al. (2010), Eq. (9). The Srivastava-Yeh
    set inverts the intrinsic Gardner law, psi=log(K/Ks)/alpha.
    """

    k = np.clip(conductivity, 1.0e-30, hydro.saturated_hydraulic_conductivity)
    coefficient = (
        1.0 / hydro.alpha if hydro.equation_set == "srivastava_yeh"
        else np.cos(hydro.slope_angle_rad) / hydro.transformed_alpha
    )
    psi = coefficient * np.log(k / hydro.saturated_hydraulic_conductivity)
    return np.minimum(psi, 0.0)


def water_content_from_pressure(
    pressure_head: np.ndarray,
    *,
    hydro: HydraulicProperties,
) -> np.ndarray:
    """Compute flow-profile water content from pressure head.

    The TRIGRS set uses slope-transformed alpha in its reconstructed
    retention curve. The Srivastava-Yeh set uses intrinsic alpha, preserving
    theta=theta_r+(theta_s-theta_r)*K/Ks in the unsaturated flow domain.
    """

    return gardner_water_content(
        pressure_head,
        alpha=hydro.alpha if hydro.equation_set == "srivastava_yeh" else hydro.transformed_alpha,
        theta_s=hydro.theta_s,
        theta_r=hydro.theta_r,
    )


def initial_profile(hydro: HydraulicProperties) -> dict[str, np.ndarray]:
    """Return the initial unsaturated profile down to the initial water table."""

    z = np.linspace(0.0, hydro.unsaturated_domain_depth, hydro.n_depths + 1)
    k = initial_conductivity(z, hydro=hydro)
    psi = pressure_from_conductivity(k, hydro=hydro)
    theta = water_content_from_pressure(psi, hydro=hydro)

    return {
        "z": z,
        "pressure_head": psi,
        "volumetric_water_content": theta,
        "hydraulic_conductivity": k,
    }


@lru_cache(maxsize=32)
def _refined_roots(alpha: float, depth: float, count: int) -> np.ndarray:
    roots = compute_roots(alpha, depth, count)
    roots.setflags(write=False)
    return roots


def _resolved_roots(t_star: float, alpha: float, depth: float, roots: np.ndarray):
    """Resolve transient series damping, or indicate a half-space limit.

    Beyond 2048 required terms, t*/(a*d)^2 is so small that the lower
    boundary's contribution is exponentially negligible. Avoid allocating an
    unbounded number of roots as an elapsed forcing time approaches zero.
    """
    length = alpha * depth
    count = int(np.ceil(length * np.sqrt((36.0 + length / 2.0) / t_star) / np.pi)) + 2
    if count <= len(roots):
        return roots
    if count > 2048:
        return None
    return _refined_roots(alpha, depth, count)


def _half_space_step(t_star: float, x: float) -> float:
    """Unit-flux step for K_t=K_xx-K_x on a half-space, initially zero.

    Integrate the inverse-Laplace impulse response for the Robin condition
    K-K_x=1 at x=0. Substitution s=t*v² removes the integrable square-root
    singularity. erfcx prevents overflow in exponential/erfc products.
    """
    def integrand(v):
        if v == 0.0:
            return 2.0 * np.sqrt(t_star / np.pi) if x == 0.0 else 0.0
        root_s = np.sqrt(t_star) * v
        exponent = x / 2.0 - x*x / (4.0 * root_s*root_s) - root_s*root_s / 4.0
        if exponent < -745.0:
            return 0.0
        argument = x / (2.0 * root_s) + root_s / 2.0
        return np.exp(exponent) * (
            2.0 * np.sqrt(t_star / np.pi) - t_star * v * erfcx(argument)
        )
    return float(quad(integrand, 0.0, 1.0, epsabs=1e-12, epsrel=1e-10)[0])


def _small_time_basal_flux(
    t_star: float,
    *,
    step_flux: float,
    hydro: HydraulicProperties,
    depth: float,
) -> float:
    """Early-time basal-flux approximation using erfc terms."""

    z = depth
    a = hydro.transformed_alpha
    qa = hydro.long_term_infiltration

    f0 = 2.0 * (step_flux - qa) * np.exp(a * z / 2.0)

    arg1 = a * z / (2.0 * np.sqrt(t_star))
    f1b = np.exp(arg1 * np.sqrt(t_star) + t_star / 4.0)
    f1c = erfc(arg1 + np.sqrt(t_star / 4.0))

    arg2 = a * (2.0 * depth - z) / (2.0 * np.sqrt(t_star))
    f2b = np.exp(arg2 * np.sqrt(t_star) + t_star / 4.0)
    f2c = erfc(arg2 + np.sqrt(t_star / 4.0))

    return float(f0 * (f1b * f1c - f2b * f2c) / 2.0)


def _basal_step_response(
    tau: float,
    *,
    step_flux: float,
    hydro: HydraulicProperties,
    depth: float,
    roots: np.ndarray,
) -> float:
    """Basal-flux response to one flux step."""

    if tau <= 0.0:
        return 0.0

    t_star = hydro.time_factor * tau
    # At the basal boundary the early-time approximation is zero. Its
    # validity also depends on travel across the unsaturated column: a small
    # t_star alone can represent a substantial response in a shallow column.
    # Retain the configured cutoff, but shorten it when a*depth is below 1.
    alpha_depth = hydro.transformed_alpha * depth
    early_cutoff = hydro.small_time_threshold * min(1.0, alpha_depth**2)
    if hydro.equation_set == "srivastava_yeh":
        roots = _resolved_roots(t_star, hydro.transformed_alpha, depth, roots)
        if roots is None:
            return float(hydro.long_term_infiltration)
    elif t_star < early_cutoff:
        return _small_time_basal_flux(
            t_star,
            step_flux=step_flux,
            hydro=hydro,
            depth=depth,
        )

    a = hydro.transformed_alpha
    qa = hydro.long_term_infiltration

    denominator = 1.0 + a * depth / 2.0 + 2.0 * a * depth * roots**2
    numerator = roots * np.sin(roots * a * depth) * np.exp(-(roots**2) * t_star)
    series_sum = np.sum(numerator / denominator)

    q = step_flux - 4.0 * (step_flux - qa) * np.exp(a * depth / 2.0) * series_sum * np.exp(
        -t_star / 4.0
    )

    return float(q)


def basal_flux_at_time(
    time: float,
    hydro: HydraulicProperties,
    rainfall: RainfallSeries,
    *,
    depth: float | None = None,
    roots: np.ndarray | None = None,
) -> float:
    """Compute basal flux at a time in hours by step superposition."""

    z_w0 = hydro.unsaturated_domain_depth if depth is None else depth
    roots = compute_roots(hydro.transformed_alpha, z_w0, hydro.n_roots) if roots is None else roots

    # Superpose changes, not absolute states: the equilibrium is included once.
    response_hydro = replace(hydro, long_term_infiltration=0.0)
    q = hydro.long_term_infiltration
    for change, t0 in zip(rainfall.flux_changes(hydro.long_term_infiltration), rainfall.breakpoints):
        q += change * _basal_step_response(
            time - t0,
            step_flux=1.0,
            hydro=response_hydro,
            depth=z_w0,
            roots=roots,
        )
    return float(np.clip(q, 0.0, hydro.saturated_hydraulic_conductivity))


def _conductivity_step_response(
    tau: float,
    z: float,
    *,
    step_flux: float,
    hydro: HydraulicProperties,
    depth: float,
    roots: np.ndarray,
) -> float:
    """Hydraulic-conductivity response to one flux step."""

    if tau <= 1.0e-10 / 3600.0:
        return 0.0

    a = hydro.transformed_alpha
    qa = hydro.long_term_infiltration
    ks = hydro.saturated_hydraulic_conductivity
    t_star = hydro.time_factor * tau

    ck0 = initial_conductivity(z, hydro=hydro, depth=depth)

    if hydro.equation_set == "srivastava_yeh":
        roots = _resolved_roots(t_star, a, depth, roots)
        if roots is None:
            return float(ck0 + (step_flux - qa) * _half_space_step(t_star, a*z))
    elif t_star < hydro.small_time_threshold:
        f0 = 2.0 * (step_flux - qa) * np.exp(a * (depth - z) / 2.0)

        arg1 = a * z / (2.0 * np.sqrt(t_star))
        f1 = erfc(arg1) - np.exp(arg1 * np.sqrt(t_star) + t_star / 4.0) * erfc(
            arg1 + np.sqrt(t_star / 4.0)
        )

        arg2 = a * (2.0 * depth - z) / (2.0 * np.sqrt(t_star))
        f2 = -erfc(arg2) + np.exp(arg2 * np.sqrt(t_star) + t_star / 4.0) * erfc(
            arg2 + np.sqrt(t_star / 4.0)
        )

        return float(ck0 + f0 * (f1 + f2) / 4.0)

    top = (
        np.sin(roots * a * (depth - z))
        * np.sin(roots * a * depth)
        * np.exp(-(roots**2) * t_star)
    )
    denominator = 1.0 + a * depth / 2.0 + 2.0 * a * depth * roots**2
    series_sum = np.sum(top / denominator)

    ck1 = step_flux - (step_flux - ks) * np.exp(-a * (depth - z))
    f1 = 2.0 * (step_flux - qa) * np.exp(a * z / 2.0) * np.exp(-t_star / 4.0)

    return float(ck1 - 2.0 * f1 * series_sum)


def profile_at_time(
    time: float,
    hydro: HydraulicProperties,
    rainfall: RainfallSeries,
    *,
    depth: float | None = None,
    roots: np.ndarray | None = None,
    evaluation_depths: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Compute the fixed-domain unsaturated response at time in hours.

    Optional evaluation_depths are vertical depths within [0, initial NF].
    The column solver constructs saturated depths separately up to maximum_depth.
    """

    z_w0 = hydro.unsaturated_domain_depth if depth is None else depth
    roots = compute_roots(hydro.transformed_alpha, z_w0, hydro.n_roots) if roots is None else roots

    z = (np.linspace(0.0, z_w0, hydro.n_depths + 1) if evaluation_depths is None
         else np.asarray(evaluation_depths, dtype=float))
    if (z.ndim != 1 or z.size == 0 or not np.all(np.isfinite(z))
            or np.any(z < 0) or np.any(z > z_w0) or np.any(np.diff(z) <= 0)):
        raise ValueError("evaluation_depths must be increasing depths within the analytical domain.")
    k = np.empty_like(z)
    response_hydro = replace(hydro, long_term_infiltration=0.0)
    changes = rainfall.flux_changes(hydro.long_term_infiltration)

    for i, zi in enumerate(z):
        if zi >= z_w0:
            k[i] = hydro.saturated_hydraulic_conductivity
            continue

        if time <= rainfall.breakpoints[0]:
            k[i] = initial_conductivity(zi, hydro=hydro, depth=z_w0)
            continue

        ki = initial_conductivity(zi, hydro=hydro, depth=z_w0)
        response_initial = initial_conductivity(zi, hydro=response_hydro, depth=z_w0)
        for change, t0 in zip(changes, rainfall.breakpoints):
            if time - t0 <= 1.0e-10 / 3600.0 or change == 0.0:
                continue
            ki += change * (_conductivity_step_response(
                time - t0,
                zi,
                step_flux=1.0,
                hydro=response_hydro,
                depth=z_w0,
                roots=roots,
            ) - response_initial)

        k[i] = ki

    k = np.clip(k, 1.0e-30, hydro.saturated_hydraulic_conductivity)
    psi = pressure_from_conductivity(k, hydro=hydro)
    theta = water_content_from_pressure(psi, hydro=hydro)

    return {
        "z": z,
        "pressure_head": psi,
        "volumetric_water_content": theta,
        "hydraulic_conductivity": k,
    }
