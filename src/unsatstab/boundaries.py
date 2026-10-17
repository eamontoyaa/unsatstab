"""Practical basal-boundary models.

The classes in this module represent two useful end-member concepts:

- `InfiniteDrainageBoundary`: a pervious/semi-infinite basal domain with high
  drainage capacity.
- `FinitePerchedBoundary`: a finite or restrictive basal domain that promotes
  perched-water buildup and water-table rise.

Both boundaries expose the same methods, so notebooks can compare them easily.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .parameters import HydraulicProperties


@dataclass(frozen=True)
class BasalBoundary:
    """Base class for practical basal-boundary behavior.

    Parameters
    ----------
    drainage_efficiency:
        Multiplier in [0, 1] controlling effective basal drainage capacity,
        not the fraction of incoming water that drains. Higher values increase
        capacity. Even 1 does not prevent accumulation if input exceeds capacity.
    name:
        Descriptive boundary name used in result metadata.
    """

    drainage_efficiency: float
    name: str

    def __post_init__(self) -> None:
        if not np.isfinite(self.drainage_efficiency) or not 0.0 <= self.drainage_efficiency <= 1.0:
            raise ValueError("drainage_efficiency must be between zero and one.")

    @property
    def has_finite_base(self) -> bool:
        """Whether maximum_depth represents a finite restrictive base."""
        return False

    def drainage_capacity(self, hydro: HydraulicProperties) -> float:
        """Return effective drainage capacity in m/s.

        This capacity adapts Baum et al. (2010), Eqs. (10a)--(10c),
        using the sign correction to Eq. (10a) in Baum and Godt (2013):

            q_capacity = eta * Ks * sin(slope)**2

        where eta is `drainage_efficiency`. This external capacity assumption
        is retained in both equation sets; it is not re-derived from the
        selected reconstructed pressure gradient.
        """
        gravity_capacity = hydro.saturated_hydraulic_conductivity * np.sin(hydro.slope_angle_rad)**2
        capacity = self.drainage_efficiency * max(gravity_capacity, 0.0)
        return float(max(capacity, 0.0))

    def excess_flux(self, basal_flux: np.ndarray, hydro: HydraulicProperties) -> np.ndarray:
        """Return the part of basal flux that accumulates instead of draining."""
        capacity = self.drainage_capacity(hydro)
        return np.maximum(np.asarray(basal_flux, dtype=float) - capacity, 0.0)

    def positive_pressure(
        self,
        z: np.ndarray,
        *,
        hydro: HydraulicProperties,
        water_table_depth: float,
    ) -> np.ndarray:
        """Pressure head below a moving water table.

        Algebraic form of Baum et al. (2010), Eq. (3b), and Baum et al.
        (2008), Eq. (3B), extended beyond their initially upward-flow case
        using the selected background gradient. This is not saturated diffusion.

        z_w is the current water-table depth, measured downward. Above the
        water table the returned pressure is zero. Below the water
        table:

            psi = beta * (z - z_w)
        """
        z = np.asarray(z, dtype=float)
        psi = np.zeros_like(z)
        below = z >= water_table_depth
        psi[below] = (
            hydro.beta
            * (z[below] - water_table_depth)
        )
        return psi


@dataclass(frozen=True)
class InfiniteDrainageBoundary(BasalBoundary):
    """Pervious/semi-infinite basal boundary.

    This boundary represents efficient drainage into a deep or highly pervious
    domain. In most moderate rainfall cases, little or no water-table rise is
    produced because drainage capacity is high.
    """

    drainage_efficiency: float = 1.0
    name: str = "infinite_drainage"


@dataclass(frozen=True)
class FinitePerchedBoundary(BasalBoundary):
    """Finite/restrictive basal boundary.

    This boundary represents a lower-permeability contact or finite saturated
    layer that restricts drainage and promotes perched-water accumulation.
    """

    drainage_efficiency: float = 0.1
    name: str = "finite_perched"

    @property
    def has_finite_base(self) -> bool:
        return True
