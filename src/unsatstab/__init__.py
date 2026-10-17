"""Spatially distributed unsaturated infiltration and slope-stability modeling."""

from ._version import __version__
from .parameters import HydraulicProperties, StrengthProperties
from .rainfall import AntecedentRainfall, RainfallSeries
from .boundaries import BasalBoundary, InfiniteDrainageBoundary, FinitePerchedBoundary
from .solver import run_single_column
from .results import ColumnResult
from .analytical import gardner_water_content
from .slope_units import (
    SlopeUnitSchema,
    SlopeUnitsResult,
    read_slope_units,
    materials_from_dataframe,
    run_slope_units,
    run_slope_units_from_gpkg,
)

__all__ = [
    "HydraulicProperties",
    "StrengthProperties",
    "RainfallSeries",
    "AntecedentRainfall",
    "BasalBoundary",
    "InfiniteDrainageBoundary",
    "FinitePerchedBoundary",
    "run_single_column",
    "ColumnResult",
    "gardner_water_content",
    "SlopeUnitSchema",
    "SlopeUnitsResult",
    "read_slope_units",
    "materials_from_dataframe",
    "run_slope_units",
    "run_slope_units_from_gpkg",
]
