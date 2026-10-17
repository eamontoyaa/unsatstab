"""Result containers for notebooks and exploratory analysis."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd


@dataclass
class ColumnResult:
    """Simulation result for one representative column.

    Attributes
    ----------
    time : numpy.ndarray
        Event-relative output hours, shape (n_times,), including simulations
        with an optional antecedent history.
    rainfall_flux, infiltration_flux, rainfall_excess_flux : numpy.ndarray
        Gross rainfall, effective infiltration after the local Ks cap, and
        rainfall discarded above that cap, in m/s (shape (n_times,)).
    basal_flux, drainage_capacity, excess_flux : numpy.ndarray
        Transfer at the effective initial NF (base of the analytical unsaturated
        domain), drainage capacity through the lower domain, and undrained
        transfer in m/s, each with shape (n_times,). ``excess_flux`` is basal drainage
        excess, not rainfall excess.
    cumulative_excess : numpy.ndarray
        Undrained storage after background correction and a non-negative floor,
        in metres, accumulated from history start. A supplied antecedent history
        carries basal storage into the event without resetting it.
    water_table_rise, water_table_depth : numpy.ndarray
        Rise Delta_h above z_w0 and current water-table depth z_w,
        in metres. A supplied antecedent history may already have NF rise at event
        time zero, depending on antecedent rainfall and drainage capacity.
    pressure_head_at_failure_plane : numpy.ndarray
        Pressure head at the selected observation plane, in metres (base by default).
    factor_of_safety_at_failure_plane : numpy.ndarray
        Dimensionless FS at the selected plane (base by default), including the cap.
    profiles : list[dict[str, numpy.ndarray]]
        One complete profile per output time, from zero to maximum_depth.
        There are n_depths + 1 uniform nodes, plus the initial water-table
        node when it lies inside the column; z and pressure heads are metres, conductivity is m/s,
        volumetric water content is m³/m³ and factor of safety is dimensionless.
    metadata : dict
        Model version, scope, boundary controls and antecedent diagnostics.
    """

    time: np.ndarray
    rainfall_flux: np.ndarray
    infiltration_flux: np.ndarray
    rainfall_excess_flux: np.ndarray
    basal_flux: np.ndarray
    drainage_capacity: np.ndarray
    excess_flux: np.ndarray
    cumulative_excess: np.ndarray
    water_table_rise: np.ndarray
    water_table_depth: np.ndarray
    pressure_head_at_failure_plane: np.ndarray
    factor_of_safety_at_failure_plane: np.ndarray
    profiles: list[dict[str, np.ndarray]]
    metadata: dict

    def time_dataframe(self) -> pd.DataFrame:
        """Return time-series outputs as a pandas DataFrame."""

        return pd.DataFrame(
            {
                "time_h": self.time,
                "rainfall_flux_m_s": self.rainfall_flux,
                "rainfall_flux_mm_h": self.rainfall_flux * 1000.0 * 3600.0,
                "infiltration_flux_m_s": self.infiltration_flux,
                "infiltration_flux_mm_h": self.infiltration_flux * 1000.0 * 3600.0,
                "rainfall_excess_flux_m_s": self.rainfall_excess_flux,
                "rainfall_excess_flux_mm_h": self.rainfall_excess_flux * 1000.0 * 3600.0,
                "basal_flux_m_s": self.basal_flux,
                "basal_flux_mm_h": self.basal_flux * 1000.0 * 3600.0,
                "drainage_capacity_m_s": self.drainage_capacity,
                "drainage_capacity_mm_h": self.drainage_capacity * 1000.0 * 3600.0,
                "excess_flux_m_s": self.excess_flux,
                "excess_flux_mm_h": self.excess_flux * 1000.0 * 3600.0,
                "cumulative_excess_m": self.cumulative_excess,
                "water_table_rise_m": self.water_table_rise,
                "water_table_depth_m": self.water_table_depth,
                "pressure_head_at_failure_plane_m": self.pressure_head_at_failure_plane,
                "factor_of_safety_at_failure_plane": self.factor_of_safety_at_failure_plane,
            }
        )

    def profile_dataframe(self) -> pd.DataFrame:
        """Return profile outputs as a long-format DataFrame."""

        frames = []
        for time, profile in zip(self.time, self.profiles):
            z = profile["z"]
            frames.append(
                pd.DataFrame(
                    {
                        "time_h": time,
                        "z_m": z,
                        "pressure_head_unsaturated_m": profile["pressure_head_unsaturated"],
                        "pressure_head_total_m": profile["pressure_head_total"],
                        "theta": profile["volumetric_water_content"],
                        "hydraulic_conductivity_m_s": profile["hydraulic_conductivity"],
                        "factor_of_safety": profile["factor_of_safety"],
                        "water_table_depth_m": profile["water_table_depth"],
                    }
                )
            )

        return pd.concat(frames, ignore_index=True)

    def arrays(self) -> dict[str, np.ndarray]:
        """Return core outputs as NumPy arrays."""

        z = self.profiles[0]["z"]
        return {
            "time": self.time,
            "z": z,
            "rainfall_flux": self.rainfall_flux,
            "infiltration_flux": self.infiltration_flux,
            "rainfall_excess_flux": self.rainfall_excess_flux,
            "basal_flux": self.basal_flux,
            "drainage_capacity": self.drainage_capacity,
            "excess_flux": self.excess_flux,
            "cumulative_excess": self.cumulative_excess,
            "water_table_rise": self.water_table_rise,
            "water_table_depth": self.water_table_depth,
            "pressure_head_at_failure_plane": self.pressure_head_at_failure_plane,
            "factor_of_safety_at_failure_plane": self.factor_of_safety_at_failure_plane,
            "pressure_head_total": np.vstack(
                [p["pressure_head_total"] for p in self.profiles]
            ),
            "theta": np.vstack([p["volumetric_water_content"] for p in self.profiles]),
            "hydraulic_conductivity": np.vstack(
                [p["hydraulic_conductivity"] for p in self.profiles]
            ),
            "factor_of_safety": np.vstack([p["factor_of_safety"] for p in self.profiles]),
        }

    def plot_time_series(self, column: str, *, ax=None, **kwargs):
        """Plot a column from `time_dataframe()` against time in hours."""

        import matplotlib.pyplot as plt

        df = self.time_dataframe()
        if column not in df:
            raise KeyError(f"{column!r} is not a valid time-series column.")

        if ax is None:
            _, ax = plt.subplots(figsize=kwargs.pop("figsize", (7, 4.5)))

        ax.plot(df["time_h"], df[column], **kwargs)
        ax.set_xlabel("Time [h]")
        ax.set_ylabel(column)
        ax.grid(True, alpha=0.3)
        return ax

    def plot_profile(self, variable: str, *, times_h: list[float], ax=None, **kwargs):
        """Plot selected depth profiles.

        Supported variables are:
        - "pressure_head_total"
        - "theta"
        - "hydraulic_conductivity"
        - "factor_of_safety"
        """

        import matplotlib.pyplot as plt

        arrays = self.arrays()
        if variable not in arrays:
            raise KeyError(f"{variable!r} is not available in result arrays.")

        if ax is None:
            _, ax = plt.subplots(figsize=kwargs.pop("figsize", (6, 5)))

        time_h_all = self.time
        for target in times_h:
            idx = int(np.argmin(np.abs(time_h_all - target)))
            ax.plot(arrays[variable][idx], arrays["z"], label=f"{time_h_all[idx]:g} h", **kwargs)

        ax.invert_yaxis()
        ax.set_ylabel("Depth [m]")
        ax.set_xlabel(variable)
        ax.grid(True, alpha=0.3)
        ax.legend()
        return ax
