"""High-level one-column solver."""

from __future__ import annotations

from dataclasses import replace, fields
import numpy as np

from ._version import __version__
from .parameters import HydraulicProperties, StrengthProperties
from .rainfall import AntecedentRainfall, RainfallSeries
from .boundaries import BasalBoundary
from .roots import compute_roots
from .analytical import basal_flux_at_time, profile_at_time
from .water_table import (
    cumulative_simpson_midpoint,
    fillable_volume_curve,
    water_table_rise_from_volume,
)
from .stability import factor_of_safety_profile
from .results import ColumnResult


def _apply_boundary_pressure(
    profile: dict[str, np.ndarray],
    *,
    hydro: HydraulicProperties,
    boundary: BasalBoundary,
    water_table_depth: float,
) -> dict[str, np.ndarray]:
    """Build total pressure, water-content, and conductivity profiles."""

    z = profile["z"]
    psi_unsat = profile["pressure_head"].copy()
    theta = profile["volumetric_water_content"].copy()
    conductivity = profile["hydraulic_conductivity"].copy()

    positive_pressure = boundary.positive_pressure(
        z,
        hydro=hydro,
        water_table_depth=water_table_depth,
    )

    pressure_total = psi_unsat.copy()
    saturated = z >= water_table_depth
    pressure_total[saturated] = positive_pressure[saturated]
    theta[saturated] = hydro.theta_s
    conductivity[saturated] = hydro.saturated_hydraulic_conductivity

    return {
        "z": z,
        "pressure_head_unsaturated": psi_unsat,
        "pressure_head_total": pressure_total,
        "volumetric_water_content": theta,
        "hydraulic_conductivity": conductivity,
        "water_table_depth": np.full_like(z, water_table_depth),
    }


def _run_history(
    *,
    hydro: HydraulicProperties,
    rainfall: RainfallSeries,
    strength: StrengthProperties,
    boundary: BasalBoundary,
    output_times: np.ndarray,
    depth_to_failure_plane: float | None = None,
    storage_time_step: float | None = None,
    antecedent_drainage_efficiency: float | None = None,
    antecedent_end_time_h: float | None = None,
) -> ColumnResult:
    """Solve a combined history, with an optional internal drainage transition."""
    output_times = np.asarray(output_times, dtype=float)
    # Validate that the requested outputs form a non-empty time vector.
    if output_times.ndim != 1 or len(output_times) == 0:
        raise ValueError("output_times must be a non-empty one-dimensional array.")
    # Validate finite, non-negative output times.
    if not np.all(np.isfinite(output_times)) or output_times[0] < 0.0:
        raise ValueError("output_times must be finite and non-negative.")
    # Validate the zero-time state required by the cumulative storage balance.
    if not np.isclose(output_times[0], 0.0, rtol=0.0, atol=1.0e-12):
        raise ValueError("output_times must begin at zero for a complete storage balance.")
    # Validate chronological ordering of all requested outputs.
    if np.any(np.diff(output_times) <= 0.0):
        raise ValueError("output_times must be strictly increasing.")
    if storage_time_step is not None and (
        not np.isfinite(storage_time_step) or storage_time_step <= 0
    ):
        raise ValueError("storage_time_step must be finite and positive.")
    if (antecedent_drainage_efficiency is None) != (antecedent_end_time_h is None):
        raise ValueError("Supply antecedent_drainage_efficiency and antecedent_end_time_h together.")
    antecedent_boundary = boundary
    if antecedent_end_time_h is not None:
        if (not np.isfinite(antecedent_end_time_h) or antecedent_end_time_h <= 0
                or antecedent_end_time_h > output_times[-1]):
            raise ValueError("antecedent_end_time_h must be finite, positive and within the simulated output interval.")
        # Reuse the boundary's validated [0, 1] efficiency and pressure settings.
        antecedent_boundary = replace(boundary, drainage_efficiency=antecedent_drainage_efficiency)

    depth = hydro.unsaturated_domain_depth
    maximum_depth = hydro.maximum_depth
    failure_depth = maximum_depth if depth_to_failure_plane is None else depth_to_failure_plane
    if not np.isfinite(failure_depth) or not 0.0 < failure_depth <= maximum_depth:
        raise ValueError("depth_to_failure_plane must be within (0, maximum_depth].")
    precipitation = rainfall
    # TRIGRS-style local infiltrability: cap each unit's rainfall input at Ks.
    # The difference is discarded here; there is no surface store or routing.
    rainfall = precipitation.capped(hydro.saturated_hydraulic_conductivity)
    roots = compute_roots(hydro.transformed_alpha, depth, hydro.n_roots)

    basal_flux = np.array(
        [
            basal_flux_at_time(t, hydro, rainfall, depth=depth, roots=roots)
            for t in output_times
        ]
    )

    capacity = boundary.drainage_capacity(hydro)
    antecedent_capacity = antecedent_boundary.drainage_capacity(hydro)

    def capacity_at(t):
        t = np.asarray(t)
        if antecedent_end_time_h is None:
            return np.full(t.shape, capacity)
        return np.where(t < antecedent_end_time_h, antecedent_capacity, capacity)

    drainage_capacity = capacity_at(output_times)
    excess_flux = np.maximum(basal_flux - drainage_capacity, 0.0)
    integration_times = np.unique(np.r_[
        output_times, rainfall.breakpoints[
            (rainfall.breakpoints > 0.0) & (rainfall.breakpoints < output_times[-1])
        ],
        [antecedent_end_time_h] if antecedent_end_time_h is not None else [],
    ])
    if storage_time_step is not None:
        segments = [np.linspace(left, right, int(np.ceil((right-left) / storage_time_step)) + 1)
                    for left, right in zip(integration_times[:-1], integration_times[1:])]
        if segments:
            integration_times = np.unique(np.concatenate(segments))
    integration_basal_flux = np.array([
        basal_flux_at_time(t, hydro, rainfall, depth=depth, roots=roots)
        for t in integration_times
    ])
    midpoint_times = 0.5 * (integration_times[:-1] + integration_times[1:])
    midpoint_basal_flux = np.array(
        [
            basal_flux_at_time(t, hydro, rainfall, depth=depth, roots=roots)
            for t in midpoint_times
        ]
    )
    # Capacity can jump at the transition. All three Simpson samples of an
    # interval use that interval's capacity, selected by its midpoint. Using
    # a right-sided endpoint capacity for the preceding interval leaks event
    # accumulation into the antecedent state.
    interval_capacity = capacity_at(midpoint_times)
    left_rate = np.maximum(integration_basal_flux[:-1] - interval_capacity, 0.0)
    mid_rate = np.maximum(midpoint_basal_flux - interval_capacity, 0.0)
    right_rate = np.maximum(integration_basal_flux[1:] - interval_capacity, 0.0)
    interval_seconds = np.diff(integration_times) * 3600.0  # convert hours to seconds
    increments = interval_seconds * (
        left_rate + 4.0 * mid_rate + right_rate
    ) / 6.0 - interval_seconds * hydro.long_term_infiltration
    # Keep the storage floor throughout the history. A dry antecedent period
    # must not create a negative storage debt that consumes the next storm.
    corrected = np.zeros_like(integration_times)
    for i, increment in enumerate(increments, start=1):
        corrected[i] = max(corrected[i-1] + increment, 0.0)
    cumulative_excess = corrected[np.searchsorted(integration_times, output_times)]

    unsaturated_profiles = [
        profile_at_time(time, hydro, rainfall, depth=depth, roots=roots)
        for time in output_times
    ]
    # Reporting geometry is independent of the initial analytical flow domain.
    # Keep that domain's profile for storage inversion, including an unreported
    # deeper water table in the infinite-domain case.
    report_z = np.linspace(0.0, maximum_depth, hydro.n_depths + 1)
    if depth < maximum_depth:
        report_z = np.unique(np.r_[report_z, depth])
    report_unsaturated = []
    for time, storage_profile in zip(output_times, unsaturated_profiles):
        if np.array_equal(report_z, storage_profile["z"]):
            report_unsaturated.append(storage_profile)
            continue
        inside = report_z <= depth
        partial = profile_at_time(time, hydro, rainfall, depth=depth, roots=roots,
                                  evaluation_depths=report_z[inside])
        profile = {"z": report_z.copy()}
        for key, saturated_value in (("pressure_head", 0.0),
                                      ("volumetric_water_content", hydro.theta_s),
                                      ("hydraulic_conductivity", hydro.saturated_hydraulic_conductivity)):
            profile[key] = np.full_like(report_z, saturated_value)
            profile[key][inside] = partial[key]
        report_unsaturated.append(profile)

    water_table_rise = np.array(
        [
            water_table_rise_from_volume(
                volume,
                profile["z"],
                fillable_volume_curve(
                    profile["z"],
                    profile["volumetric_water_content"],
                    theta_s=hydro.theta_s,
                ),
                initial_water_table_depth=depth,
            )
            for volume, profile in zip(cumulative_excess, unsaturated_profiles)
        ]
    )
    water_table_depth = depth - water_table_rise

    pressure_head_at_failure = np.zeros_like(output_times)

    profiles = []
    fs_at_failure = np.zeros_like(output_times)
    rainfall_flux = np.array([precipitation.flux_at(t) for t in output_times])
    infiltration_flux = np.array([rainfall.flux_at(t) for t in output_times])
    rainfall_excess_flux = np.maximum(rainfall_flux - infiltration_flux, 0.0)

    for i, unsat in enumerate(report_unsaturated):
        total = _apply_boundary_pressure(
            unsat,
            hydro=hydro,
            boundary=boundary,
            water_table_depth=water_table_depth[i],
        )
        fs = factor_of_safety_profile(
            total["z"],
            total["pressure_head_total"],
            total["volumetric_water_content"],
            hydro=hydro,
            strength=strength,
        )
        total["factor_of_safety"] = fs
        profiles.append(total)
        fs_at_failure[i] = np.interp(failure_depth, total["z"], fs)
        pressure_head_at_failure[i] = np.interp(
            failure_depth, total["z"], total["pressure_head_total"]
        )

    return ColumnResult(
        time=output_times,
        rainfall_flux=rainfall_flux,
        infiltration_flux=infiltration_flux,
        rainfall_excess_flux=rainfall_excess_flux,
        basal_flux=basal_flux,
        drainage_capacity=drainage_capacity,
        excess_flux=excess_flux,
        cumulative_excess=cumulative_excess,
        water_table_rise=water_table_rise,
        water_table_depth=water_table_depth,
        pressure_head_at_failure_plane=pressure_head_at_failure,
        factor_of_safety_at_failure_plane=fs_at_failure,
        profiles=profiles,
        metadata={
            "package": "unsatstab",
            "version": __version__,
            "equation_set": hydro.equation_set,
            "model_scope": "single_column",
            "boundary": boundary.name,
            "drainage_efficiency": boundary.drainage_efficiency,
            "long_term_infiltration_m_s": float(hydro.long_term_infiltration),
            "surface_infiltration_cap_m_s": float(hydro.saturated_hydraulic_conductivity),
            "surface_excess_handling": "discarded_locally_without_storage_or_routing",
            "storage_time_step": storage_time_step,
            "antecedent_drainage_efficiency": antecedent_drainage_efficiency,
            "antecedent_end_time_h": antecedent_end_time_h,
            "initial_water_table_depth_m": float(depth),
            "maximum_depth_m": float(maximum_depth),
            "n_roots": int(hydro.n_roots),
            "n_depths": int(hydro.n_depths),
            "analytical_domain_depth_m": float(depth),
            "failure_plane_depth_m": float(failure_depth),
            "note": (
                "One-dimensional analytical infiltration response with a storage-based "
                "moving-water-table approximation."
            ),
        },
    )


def run_single_column(
    *,
    hydro: HydraulicProperties,
    rainfall: RainfallSeries,
    strength: StrengthProperties,
    boundary: BasalBoundary,
    output_times: np.ndarray,
    depth_to_failure_plane: float | None = None,
    antecedent: AntecedentRainfall | None = None,
    storage_time_step: float | None = None,
) -> ColumnResult:
    """Simulate a triggering event, optionally preceded by antecedent rain.

    Parameters
    ----------
    hydro : HydraulicProperties
        Hydraulic properties and geometry in SI units. maximum_depth defines
        the reporting column; initial water-table depth is independent. ``long_term_infiltration`` is
        the equilibrium background before any antecedent period.
    rainfall : RainfallSeries
        Gross triggering precipitation, with rainfall breakpoints relative to the trigger.
        Infiltration is capped at local Ks; rainfall excess is discarded.
    strength : StrengthProperties
        Effective cohesion, friction angle and unit weights.
    boundary : BasalBoundary
        Event drainage capacity and newly saturated pressure assumptions.
    output_times : numpy.ndarray
        Instants at which time series and depth profiles are returned, in hours
        relative to the triggering event. Must be finite, strictly increasing,
        and start at zero; the last value sets the simulated event horizon.
        These nodes also seed the storage integration grid. Internal storage
        refinement does not add returned profiles.
    depth_to_failure_plane : float, optional
        Vertical depth in (0, maximum_depth] for extracting pressure and FS
        time series. Defaults to maximum_depth. It does not change flow,
        storage, or the full-column hydraulic and stability profiles.
    antecedent : AntecedentRainfall, optional
        Previous rainfall history and its drainage efficiency. Use
        ``AntecedentRainfall.uniform(total_rainfall_mm, duration, ...)``
        for a uniform history or wrap a ``RainfallSeries`` for a hyetograph.
        Moisture, conductivity and basal storage continue into the event.
        Output times stay event-relative; water-table rise and storage reference
        the beginning of the complete history.
    storage_time_step : float, optional
        Maximum interval in hours for integrating the basal storage balance.
        Intervals longer than this are subdivided; each uses Simpson endpoint
        and midpoint flux evaluations. Smaller values can improve retained
        storage, water-table rise and FS accuracy at greater computational cost.
        This does not set rainfall intervals or saved output times. Default
        None adds no subdivision without antecedent; with antecedent, None
        selects 1 h. Output nodes, rainfall breakpoints, and any drainage
        transition are always included in the storage integration grid.

    Returns
    -------
    ColumnResult
        Event time series, profiles and model/budget diagnostics. Depth
        increases downward. FS at the declared plane is distinct from a
        minimum over all depths. Antecedent conditions can change FS at zero.

    Notes
    -----
    No evapotranspiration, interception or numerical Richards solver is used.
    The analytical time factor does not change with initial wetness. Basal
    capacity is eta * Ks * sin(slope)**2, not a drained rainfall fraction.
    Even eta=1 does not enforce fixed NF. A different antecedent efficiency
    is an explicit boundary assumption. The approximate moving-NF saturation
    overlay, surface cap and lack of separate aquifer recession still apply.
    """
    if not isinstance(rainfall, RainfallSeries):
        raise TypeError("rainfall must be a RainfallSeries.")
    if not isinstance(hydro, HydraulicProperties):
        raise TypeError("hydro must be HydraulicProperties.")
    if hydro.maximum_depth is None or hydro.slope_angle_deg is None:
        raise ValueError("Resolved slope_angle_deg and maximum_depth are required for a column run.")
    if not isinstance(boundary, BasalBoundary):
        raise TypeError("boundary must be a BasalBoundary.")
    requested_water_table_depth = hydro.unsaturated_domain_depth
    effective_depth = (min(requested_water_table_depth, hydro.maximum_depth)
                       if boundary.has_finite_base else requested_water_table_depth)
    hydro = replace(hydro, initial_water_table_depth=effective_depth)
    if antecedent is None:
        result = _run_history(
            hydro=hydro, rainfall=rainfall, strength=strength, boundary=boundary,
            output_times=output_times, depth_to_failure_plane=depth_to_failure_plane,
            storage_time_step=storage_time_step,
        )
    else:
        if not isinstance(antecedent, AntecedentRainfall):
            raise TypeError("antecedent must be an AntecedentRainfall or None.")
        result = _run_with_antecedent(
            hydro=hydro, rainfall=rainfall, strength=strength, boundary=boundary,
            output_times=output_times, depth_to_failure_plane=depth_to_failure_plane,
            antecedent=antecedent,
            storage_time_step=1.0 if storage_time_step is None else storage_time_step,
        )
    result.metadata.update({
        "requested_initial_water_table_depth_m": float(requested_water_table_depth),
        "initial_water_table_depth_was_clamped": bool(effective_depth != requested_water_table_depth),
        "finite_base": bool(boundary.has_finite_base),
    })
    return result


def _run_with_antecedent(
    *,
    hydro: HydraulicProperties,
    rainfall: RainfallSeries,
    strength: StrengthProperties,
    boundary: BasalBoundary,
    output_times: np.ndarray,
    depth_to_failure_plane: float | None = None,
    antecedent: AntecedentRainfall,
    storage_time_step: float,
) -> ColumnResult:
    """Return event-relative outputs from the full analytical forcing history."""
    event_times = np.asarray(output_times, dtype=float)
    if (event_times.ndim != 1 or len(event_times) == 0
            or not np.all(np.isfinite(event_times)) or event_times[0] != 0.0
            or np.any(np.diff(event_times) <= 0)):
        raise ValueError("output_times must be finite, strictly increasing event-relative hours beginning at zero.")
    antecedent_efficiency = (
        boundary.drainage_efficiency if antecedent.drainage_efficiency is None
        else antecedent.drainage_efficiency
    )
    antecedent_boundary = replace(boundary, drainage_efficiency=antecedent_efficiency)
    antecedent_total_mm = antecedent.total_rainfall_mm
    history = antecedent.prepend_to(rainfall)
    duration_h = antecedent.duration
    duration_s = duration_h * 3600.0  # convert hours to seconds for SI flux budgets
    result = _run_history(
        hydro=hydro, rainfall=history, strength=strength, boundary=boundary,
        output_times=np.r_[0.0, duration_h + event_times],
        depth_to_failure_plane=depth_to_failure_plane,
        storage_time_step=storage_time_step,
        antecedent_drainage_efficiency=antecedent_efficiency,
        antecedent_end_time_h=duration_h,
    )
    # Drop the prehistory start output, keeping the event-start state intact.
    event_fields = {
        field.name: getattr(result, field.name)[1:]
        for field in fields(ColumnResult)
        if field.name not in {"metadata", "time"}
    }
    event_result = replace(result, time=event_times.copy(), **event_fields)
    # Quantify the analytical unsaturated response before any moving-NF overlay.
    antecedent_forcing = history.capped(hydro.saturated_hydraulic_conductivity)
    roots = compute_roots(hydro.transformed_alpha, hydro.unsaturated_domain_depth, hydro.n_roots)
    raw_start = profile_at_time(0.0, hydro, antecedent_forcing, roots=roots)
    raw_end = profile_at_time(duration_h, hydro, antecedent_forcing, roots=roots)
    raw_water_change = np.trapezoid(
        raw_end["volumetric_water_content"] - raw_start["volumetric_water_content"],
        raw_end["z"],
    )
    budget_step_h = storage_time_step
    budget_times = np.unique(np.r_[
        np.linspace(0.0, duration_h, int(np.ceil(duration_h / budget_step_h)) + 1),
        antecedent.rainfall.breakpoints,
    ])
    budget_midpoints = 0.5 * (budget_times[:-1] + budget_times[1:])
    budget_flux = np.array([
        basal_flux_at_time(t, hydro, antecedent_forcing, roots=roots) for t in budget_times
    ])
    budget_mid_flux = np.array([
        basal_flux_at_time(t, hydro, antecedent_forcing, roots=roots) for t in budget_midpoints
    ])
    basal_transfer = cumulative_simpson_midpoint(budget_flux, budget_mid_flux, budget_times * 3600.0)[-1]  # convert hours to seconds for m/s integration
    antecedent_capacity = antecedent_boundary.drainage_capacity(hydro)
    capacity_drainage = cumulative_simpson_midpoint(
        np.minimum(budget_flux, antecedent_capacity),
        np.minimum(budget_mid_flux, antecedent_capacity), budget_times * 3600.0,  # convert hours to seconds
    )[-1]
    basal_storage = event_result.cumulative_excess[0]
    # Storage follows the inherited background correction and zero floor.
    # Report outflow allocated by that balance; it is independent of theta.
    basal_drainage = basal_transfer - basal_storage
    # The budget covers the full analytical domain, independently of the
    # output-depth window or selected observation plane.
    start_profile = _apply_boundary_pressure(raw_start, hydro=hydro, boundary=boundary,
                                             water_table_depth=hydro.unsaturated_domain_depth)
    end_profile = _apply_boundary_pressure(raw_end, hydro=hydro, boundary=boundary,
                                           water_table_depth=event_result.water_table_depth[0])
    profile_water_change = np.trapezoid(
        end_profile["volumetric_water_content"] - start_profile["volumetric_water_content"],
        raw_end["z"],
    )
    antecedent_infiltration = antecedent.rainfall.capped(hydro.saturated_hydraulic_conductivity)
    effective_volume = float(np.dot(
        np.diff(antecedent_infiltration.breakpoints), antecedent_infiltration.fluxes,
    )) * 3600.0  # convert hours to seconds for the m/s input budget
    # Report closure error, rather than redefining drainage to force an exact
    # balance. Moving-water-table reconstruction and finite sampling are approximate.
    balance_residual = effective_volume - basal_drainage - profile_water_change
    event_result.metadata = {
        **result.metadata,
        "antecedent_rainfall_mm": float(antecedent_total_mm),
        "antecedent_infiltration_mm": float(effective_volume * 1000.0),
        "antecedent_duration_h": antecedent.duration,
        "antecedent_flux_m_s": float(antecedent_total_mm / 1000.0 / duration_s),
        "antecedent_rainfall_flux_m_s": float(antecedent_total_mm / 1000.0 / duration_s),
        "antecedent_infiltration_flux_m_s": float(effective_volume / duration_s),
        "event_start_time_h_in_history": float(duration_h),
        "time_reference": "event_start",
        "antecedent_drainage_efficiency": float(antecedent_efficiency),
        "antecedent_drainage_capacity_m_s": float(antecedent_capacity),
        "antecedent_lower_boundary": "capacity_limited_drainage",
        "antecedent_effective_infiltration_mm": float(effective_volume * 1000.0),
        "antecedent_profile_water_change_mm": float(profile_water_change * 1000.0),
        "antecedent_unsaturated_response_water_change_mm": float(raw_water_change * 1000.0),
        "antecedent_rainfall_excess_mm": float(antecedent_total_mm - effective_volume * 1000.0),
        "antecedent_basal_transfer_mm": float(basal_transfer * 1000.0),
        "antecedent_basal_storage_mm": float(basal_storage * 1000.0),
        "antecedent_storage_above_pore_capacity_mm": float(max(
            basal_storage - fillable_volume_curve(
                raw_end["z"], raw_end["volumetric_water_content"], theta_s=hydro.theta_s,
            )[0], 0.0,
        ) * 1000.0),
        "antecedent_capacity_drainage_mm": float(capacity_drainage * 1000.0),
        "antecedent_drainage_mm": float(basal_drainage * 1000.0),
        "antecedent_background_correction_mm": float((basal_drainage - capacity_drainage) * 1000.0),
        "antecedent_water_balance_residual_mm": float(balance_residual * 1000.0),
        "storage_reference": "antecedent_start",
        "event_initial_water_table_depth_m": float(event_result.water_table_depth[0]),
        "event_initial_factor_of_safety": float(event_result.factor_of_safety_at_failure_plane[0]),
    }
    return event_result
