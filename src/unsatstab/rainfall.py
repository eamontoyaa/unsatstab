"""Rainfall and infiltration forcing."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class RainfallSeries:
    """Piecewise-constant rainfall-intensity series.

    The object stores time breakpoints in hours and one rainfall rate per
    interval in m/s. The column solver caps this rate at the local ``Ks`` to
    obtain effective infiltration; the excess is discarded locally.
    """

    breakpoints: np.ndarray
    fluxes: np.ndarray

    def __post_init__(self) -> None:
        breakpoints = np.asarray(self.breakpoints, dtype=float)
        fluxes = np.asarray(self.fluxes, dtype=float)

        if breakpoints.ndim != 1 or fluxes.ndim != 1:
            raise ValueError("breakpoints and fluxes must be one-dimensional.")
        if len(breakpoints) != len(fluxes) + 1:
            raise ValueError("breakpoints must have length len(fluxes) + 1.")
        if len(fluxes) == 0:
            raise ValueError("At least one rainfall interval is required.")
        if not np.all(np.isfinite(breakpoints)) or not np.all(np.isfinite(fluxes)):
            raise ValueError("breakpoints and fluxes must contain only finite values.")
        if breakpoints[0] < 0.0:
            raise ValueError("breakpoints cannot start before time zero.")
        if np.any(np.diff(breakpoints) <= 0.0):
            raise ValueError("breakpoints must be strictly increasing.")
        if np.any(fluxes < 0.0):
            raise ValueError("fluxes must be non-negative.")

        object.__setattr__(self, "breakpoints", breakpoints)
        object.__setattr__(self, "fluxes", fluxes)

    @classmethod
    def constant(
        cls,
        intensity: float,
        duration: float,
        *,
        intensity_units: str = "mm/h",
        duration_units: str = "h",
    ) -> "RainfallSeries":
        """Create a constant event; duration defaults to hours.

        Alternative duration_units (min or s) are explicitly converted to hours.
        """

        if not np.isfinite(intensity) or intensity < 0.0:
            raise ValueError("intensity must be finite and non-negative.")
        if not np.isfinite(duration) or duration <= 0.0:
            raise ValueError("duration must be finite and positive.")

        if intensity_units == "mm/h":
            flux = intensity / 1000.0 / 3600.0
        elif intensity_units == "m/s":
            flux = intensity
        else:
            raise ValueError("intensity_units must be 'mm/h' or 'm/s'.")

        if duration_units == "h":
            duration_h = duration
        elif duration_units == "min":
            duration_h = duration / 60.0  # convert minutes to hours
        elif duration_units == "s":
            duration_h = duration / 3600.0  # convert seconds to hours
        else:
            raise ValueError("duration_units must be 'h', 'min' or 's'.")

        return cls(
            breakpoints=np.array([0.0, duration_h], dtype=float),
            fluxes=np.array([flux], dtype=float),
        )

    @classmethod
    def from_intervals(
        cls,
        intensities: np.ndarray,
        durations: np.ndarray,
        *,
        intensity_units: str = "mm/h",
        duration_units: str = "h",
    ) -> "RainfallSeries":
        """Create rainfall from successive intensity/duration pairs.

        Both arrays must be one-dimensional, non-empty, and equally sized.
        Intensities must be finite and non-negative; durations must be finite
        and positive. The history starts at zero and its breakpoints are the
        cumulative durations. Include zero-intensity intervals for dry gaps.

        Accepted intensity units are ``mm/h`` and ``m/s``; duration units are
        ``h``, ``min`` and ``s``. Stored values always use hours and m/s.
        """
        rates = np.asarray(intensities, dtype=float)
        lengths = np.asarray(durations, dtype=float)
        if (rates.ndim != 1 or lengths.ndim != 1 or rates.size == 0
                or rates.shape != lengths.shape):
            raise ValueError("intensities and durations must be non-empty one-dimensional arrays of equal length.")
        if not np.all(np.isfinite(rates)) or np.any(rates < 0):
            raise ValueError("intensities must be finite and non-negative.")
        if not np.all(np.isfinite(lengths)) or np.any(lengths <= 0):
            raise ValueError("durations must be finite and positive.")
        if intensity_units == "mm/h":
            rates = rates / 3.6e6  # convert mm/h to m/s
        elif intensity_units != "m/s":
            raise ValueError("intensity_units must be 'mm/h' or 'm/s'.")
        duration_factors = {"h": 1.0, "min": 1.0 / 60.0, "s": 1.0 / 3600.0}
        if duration_units not in duration_factors:
            raise ValueError("duration_units must be 'h', 'min' or 's'.")
        lengths = lengths * duration_factors[duration_units]  # convert durations to hours
        return cls(np.r_[0.0, np.cumsum(lengths)], rates)

    @classmethod
    def from_intensities(
        cls,
        intensities: np.ndarray,
        time_step: float,
        *,
        intensity_units: str = "mm/h",
        time_units: str = "h",
        start_time: float = 0.0,
    ) -> "RainfallSeries":
        """Create regular intervals; time_step and start_time default to hours.

        time_units converts the interval spacing only; start_time is always hours.
        """

        values = np.asarray(intensities, dtype=float)
        if values.ndim != 1 or len(values) == 0:
            raise ValueError("intensities must be a non-empty one-dimensional array.")
        if not np.isfinite(time_step) or time_step <= 0.0:
            raise ValueError("time_step must be finite and positive.")
        if not np.isfinite(start_time) or start_time < 0.0:
            raise ValueError("start_time must be finite and non-negative.")

        if intensity_units == "mm/h":
            fluxes = values / 1000.0 / 3600.0
        elif intensity_units == "m/s":
            fluxes = values
        else:
            raise ValueError("intensity_units must be 'mm/h' or 'm/s'.")

        factors = {"min": 1.0 / 60.0, "h": 1.0, "s": 1.0 / 3600.0}
        if time_units not in factors:
            raise ValueError("time_units must be 'min', 'h', or 's'.")
        step_h = time_step * factors[time_units]
        breakpoints = start_time + np.arange(len(values) + 1, dtype=float) * step_h
        return cls(breakpoints=breakpoints, fluxes=fluxes)

    def capped(self, maximum_flux: float) -> "RainfallSeries":
        """Return a copy whose rainfall rate cannot exceed ``maximum_flux``.

        The solver uses this operation with local ``Ks`` to create effective
        infiltration; removed rainfall is reported as ``rainfall_excess_flux``.
        """

        if not np.isfinite(maximum_flux) or maximum_flux < 0.0:
            raise ValueError("maximum_flux must be finite and non-negative.")
        return RainfallSeries(self.breakpoints.copy(), np.minimum(self.fluxes, maximum_flux))

    def with_antecedent_rainfall(
        self, total_rainfall_mm: float, duration: float, *, duration_units: str = "h",
    ) -> "RainfallSeries":
        """Return combined forcing for uniform antecedent rain and this event.

        For event-relative simulation with stage-specific drainage, pass an
        ``AntecedentRainfall`` object to ``run_single_column`` instead.
        """
        return AntecedentRainfall.uniform(total_rainfall_mm, duration, duration_units=duration_units).prepend_to(self)

    def flux_changes(self, background_flux: float) -> np.ndarray:
        """Flux increments relative to the equilibrium before the first interval.

        After the final breakpoint the forcing is zero. Each increment is
        applied once by the analytical step response, preserving continuity
        of the internal state at forcing transitions.
        """
        return np.diff(np.r_[background_flux, self.fluxes, 0.0])

    def flux_at(self, time: float) -> float:
        """Return instantaneous rainfall rate at time in hours using intervals [start, end).

        At an interior breakpoint the next interval is active. At the final
        breakpoint the event has ended and flux is zero. This is not an
        interval average; plot samples with steps rather than linear segments.
        """
        if time < self.breakpoints[0] or time >= self.breakpoints[-1]:
            return 0.0

        idx = np.searchsorted(self.breakpoints, time, side="right") - 1
        idx = min(max(idx, 0), len(self.fluxes) - 1)
        return float(self.fluxes[idx])

    @property
    def duration(self) -> float:
        """Total duration in hours."""
        return float(self.breakpoints[-1] - self.breakpoints[0])

    @property
    def max_flux(self) -> float:
        """Maximum rainfall rate in m/s."""
        return float(np.max(self.fluxes))

    @property
    def max_flux_mm_h(self) -> float:
        """Maximum rainfall rate in mm/h."""
        return float(self.max_flux * 1000.0 * 3600.0)


@dataclass(frozen=True)
class AntecedentRainfall:
    """Rainfall history preceding a trigger, with optional basal drainage.

    Parameters
    ----------
    rainfall : RainfallSeries
        Gross antecedent rainfall beginning at zero. Its final breakpoint is
        the trigger's time origin, even if the triggering series starts later.
        Include explicit zero-rain intervals for dry periods within this window.
    drainage_efficiency : float, optional
        Capacity multiplier in [0, 1] during the antecedent. ``None`` uses the
        event boundary's efficiency. This is not a fraction of rainfall drained;
        even 1 cannot guarantee a fixed water table under intense input.

    Notes
    -----
    Use ``uniform`` when only total rainfall and duration are known. The solver
    caps all rates at local Ks, without evapotranspiration or interception.
    This object does not change the pre-history ``long_term_infiltration``.
    """

    rainfall: RainfallSeries
    drainage_efficiency: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.rainfall, RainfallSeries):
            raise TypeError("antecedent rainfall must be a RainfallSeries.")
        if self.rainfall.breakpoints[0] != 0.0:
            raise ValueError("antecedent rainfall must begin at zero; include dry intervals explicitly.")
        if self.drainage_efficiency is not None and (
            not np.isfinite(self.drainage_efficiency)
            or not 0.0 <= self.drainage_efficiency <= 1.0
        ):
            raise ValueError("drainage_efficiency must be finite and between zero and one.")

    @classmethod
    def uniform(
        cls, total_rainfall_mm: float, duration: float,
        *, drainage_efficiency: float | None = None, duration_units: str = "h",
    ) -> "AntecedentRainfall":
        """Uniform gross rainfall from depth in mm and duration in hours.

        Other duration units require an explicit duration_units: days, min or s.
        The canonical stored history and duration always use hours.
        """
        if not np.isfinite(total_rainfall_mm) or total_rainfall_mm < 0.0:
            raise ValueError("total_rainfall_mm must be finite and non-negative.")
        if not np.isfinite(duration) or duration <= 0.0:
            raise ValueError("duration must be finite and positive.")
        factors = {"h": 1.0, "days": 24.0, "min": 1.0 / 60.0, "s": 1.0 / 3600.0}
        if duration_units not in factors:
            raise ValueError("duration_units must be 'h', 'days', 'min' or 's'.")
        duration_h = duration * factors[duration_units]  # convert duration to hours
        return cls(RainfallSeries.constant(total_rainfall_mm / duration_h, duration_h),
                   drainage_efficiency=drainage_efficiency)

    @property
    def duration(self) -> float:
        """Antecedent period in hours."""
        return self.rainfall.duration

    @property
    def duration_days(self) -> float:
        """Explicit calendar conversion of the canonical duration to days."""
        return self.duration / 24.0

    @property
    def total_rainfall_mm(self) -> float:
        """Integrated gross precipitation in mm, including all intervals."""
        return float(np.dot(np.diff(self.rainfall.breakpoints), self.rainfall.fluxes) * 3600.0 * 1000.0)

    def prepend_to(self, event: RainfallSeries) -> RainfallSeries:
        """Combine histories while preserving any dry gap before the event."""
        if not isinstance(event, RainfallSeries):
            raise TypeError("event must be a RainfallSeries.")
        shifted = event.breakpoints + self.duration
        if event.breakpoints[0] == 0.0:
            return RainfallSeries(
                np.r_[self.rainfall.breakpoints, shifted[1:]],
                np.r_[self.rainfall.fluxes, event.fluxes],
            )
        return RainfallSeries(
            np.r_[self.rainfall.breakpoints, shifted],
            np.r_[self.rainfall.fluxes, 0.0, event.fluxes],
        )
