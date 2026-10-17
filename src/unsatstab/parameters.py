"""Parameter objects for one-column unsaturated infiltration and stability models."""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True, kw_only=True)
class HydraulicProperties:
    """Hydraulic and geometric properties for one representative column.

    Parameters
    ----------
    slope_angle_deg:
        Slope angle in degrees.
    maximum_depth:
        Maximum vertical depth of the reported column in metres. It is the
        finite basal boundary depth, or the evaluation limit for an infinite
        basal domain. Omit geometry when using this object as a material template.
    saturated_hydraulic_conductivity:
        Saturated hydraulic conductivity, Ks, in m/s.
    alpha:
        Gardner exponential hydraulic parameter in 1/m.
    theta_s:
        Saturated volumetric water content.
    theta_r:
        Residual volumetric water content.
    long_term_infiltration:
        Equilibrium background infiltration flux in m/s BEFORE the simulated
        history. A finite antecedent rainfall window belongs in RainfallSeries,
        not in this parameter; matching the two rates produces no initial
        wetting perturbation of the unsaturated profile.
    initial_water_table_depth:
        Initial water-table depth z_w0 in metres. If omitted, it defaults to the
        maximum depth. A finite boundary clamps a deeper requested water table
        to the maximum depth; the solver records requested and effective depths.
    n_roots:
        Number of roots used in the analytical unsaturated series.
    n_depths:
        Number of vertical intervals used in profile outputs.
    small_time_threshold:
        TRIGRS-set dimensionless threshold for early-time approximations. Basal flux
        additionally multiplies this cutoff by min(1, (a * z_w0)**2) to avoid
        delaying the response of shallow columns; profile conductivity uses
        the unscaled threshold. The Srivastava-Yeh set instead resolves the
        transient series, using the half-space limit at extremely short times.
    equation_set:
        ``"baum_etal"`` retains the TRIGRS-style pressure and retention
        conversions. ``"srivastava_yeh"`` uses the intrinsic Gardner closure
        and a background pressure gradient consistent with the transformed
        Darcy law. The latter extends the horizontal 1991 formulation to a
        slope; it is not a different saturated diffusion solver.
    """

    slope_angle_deg: float | None = None
    maximum_depth: float | None = None
    saturated_hydraulic_conductivity: float
    alpha: float
    theta_s: float
    theta_r: float
    long_term_infiltration: float = 0.0
    initial_water_table_depth: float | None = None
    n_roots: int = 250
    n_depths: int = 250
    small_time_threshold: float = 1.0e-4
    equation_set: str = "srivastava_yeh"

    def __post_init__(self) -> None:
        if self.equation_set not in ("baum_etal", "srivastava_yeh"):
            raise ValueError("equation_set must be 'baum_etal' or 'srivastava_yeh'.")
        values = {
            "saturated_hydraulic_conductivity": self.saturated_hydraulic_conductivity,
            "alpha": self.alpha,
            "theta_s": self.theta_s,
            "theta_r": self.theta_r,
            "long_term_infiltration": self.long_term_infiltration,
            "small_time_threshold": self.small_time_threshold,
        }
        if not all(np.isfinite(value) for value in values.values()):
            raise ValueError("All hydraulic parameters must be finite.")
        if self.slope_angle_deg is not None and (not np.isfinite(self.slope_angle_deg) or not 0.0 <= self.slope_angle_deg < 90.0):
            raise ValueError("slope_angle_deg must be in [0, 90).")
        if self.maximum_depth is not None and (not np.isfinite(self.maximum_depth) or self.maximum_depth <= 0.0):
            raise ValueError("maximum_depth must be positive.")
        if self.initial_water_table_depth is not None and (
            not np.isfinite(self.initial_water_table_depth) or self.initial_water_table_depth <= 0.0
        ):
            raise ValueError("initial_water_table_depth must be finite and positive.")
        if self.saturated_hydraulic_conductivity <= 0.0:
            raise ValueError("saturated_hydraulic_conductivity must be positive.")
        if self.alpha <= 0.0:
            raise ValueError("alpha must be positive.")
        if not 0.0 <= self.theta_r < self.theta_s <= 1.0:
            raise ValueError("Water contents must satisfy 0 <= theta_r < theta_s <= 1.")
        if not 0.0 <= self.long_term_infiltration <= self.saturated_hydraulic_conductivity:
            raise ValueError(
                "long_term_infiltration must be between zero and saturated hydraulic conductivity."
            )
        if not isinstance(self.n_roots, int) or self.n_roots <= 0:
            raise ValueError("n_roots must be a positive integer.")
        if not isinstance(self.n_depths, int) or self.n_depths <= 0:
            raise ValueError("n_depths must be a positive integer.")
        if self.small_time_threshold <= 0.0:
            raise ValueError("small_time_threshold must be positive.")

    @property
    def unsaturated_domain_depth(self) -> float:
        """Requested flow-reference depth; finite-base adjustment occurs in the solver.

        Material templates without geometry cannot resolve this depth.
        """
        if self.initial_water_table_depth is None:
            if self.maximum_depth is None:
                raise ValueError("maximum_depth is required to resolve column geometry.")
            return float(self.maximum_depth)
        return float(self.initial_water_table_depth)

    @property
    def slope_angle_rad(self) -> float:
        """Slope angle in radians."""
        if self.slope_angle_deg is None:
            raise ValueError("slope_angle_deg is required to resolve column geometry.")
        return float(np.deg2rad(self.slope_angle_deg))

    @property
    def transformed_alpha(self) -> float:
        """Slope-corrected Gardner parameter, alpha*cos(slope)^2."""
        return float(self.alpha * np.cos(self.slope_angle_rad) ** 2)

    @property
    def time_factor(self) -> float:
        """Time factor in h⁻¹: multiplying elapsed hours gives dimensionless t*."""
        return float(
            3600.0  # convert the SI inverse-second coefficient to inverse hours
            * self.transformed_alpha
            * self.saturated_hydraulic_conductivity
            / (self.theta_s - self.theta_r)
        )

    @property
    def beta(self) -> float:
        """Background pressure-head gradient parameter.

        The TRIGRS setting follows:

            beta = cos(slope)^2 - Izlt / Ks

        The Srivastava-Yeh setting uses cos(slope)^2 * (1 - Izlt / Ks),
        so saturated Darcy flux at this gradient equals the background input.
        """
        if self.equation_set == "srivastava_yeh":
            return float(np.cos(self.slope_angle_rad) ** 2 * (
                1.0 - self.long_term_infiltration / self.saturated_hydraulic_conductivity
            ))
        return float(
            np.cos(self.slope_angle_rad) ** 2
            - self.long_term_infiltration / self.saturated_hydraulic_conductivity
        )

    @property
    def ks_mm_h(self) -> float:
        """Saturated hydraulic conductivity in mm/h."""
        return float(self.saturated_hydraulic_conductivity * 1000.0 * 3600.0)

    @property
    def izlt_mm_h(self) -> float:
        """Long-term infiltration flux in mm/h."""
        return float(self.long_term_infiltration * 1000.0 * 3600.0)


@dataclass(frozen=True)
class StrengthProperties:
    """Strength and unit-weight parameters for infinite-slope stability.

    Parameters
    ----------
    cohesion:
        Effective cohesion c' in Pa.
    friction_angle_deg:
        Effective friction angle phi' in degrees.
    unit_weight_soil:
        Total soil unit weight in N/m³.
    unit_weight_water:
        Water unit weight in N/m³.
    max_factor_of_safety:
        Upper cap for FS values, useful for plotting.
    """

    cohesion: float
    friction_angle_deg: float
    unit_weight_soil: float
    unit_weight_water: float = 9_810.0
    max_factor_of_safety: float = 10.0

    def __post_init__(self) -> None:
        values = (
            self.cohesion,
            self.friction_angle_deg,
            self.unit_weight_soil,
            self.unit_weight_water,
            self.max_factor_of_safety,
        )
        if not all(np.isfinite(value) for value in values):
            raise ValueError("All strength parameters must be finite.")
        if self.cohesion < 0.0:
            raise ValueError("cohesion must be non-negative.")
        if not 0.0 <= self.friction_angle_deg < 90.0:
            raise ValueError("friction_angle_deg must be in [0, 90).")
        if self.unit_weight_soil <= 0.0 or self.unit_weight_water <= 0.0:
            raise ValueError("Unit weights must be positive.")
        if self.max_factor_of_safety <= 0.0:
            raise ValueError("max_factor_of_safety must be positive.")

    @property
    def friction_angle_rad(self) -> float:
        """Effective friction angle in radians."""
        return float(np.deg2rad(self.friction_angle_deg))
