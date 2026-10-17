"""Spatial simulations over slope units with reusable material properties."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from collections.abc import Mapping
from typing import Hashable

import re

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

from .boundaries import BasalBoundary, FinitePerchedBoundary, InfiniteDrainageBoundary
from .parameters import HydraulicProperties, StrengthProperties
from .rainfall import AntecedentRainfall, RainfallSeries
from .results import ColumnResult
from .solver import run_single_column


@dataclass(frozen=True)
class SlopeUnitSchema:
    """Column names for required geometry and optional per-unit overrides.

    Material properties are inherited from the material_id mapping. Optional
    non-null row properties override those defaults. Rainfall pairs int_0,
    dur_0, ... use mm/h and hours unless run-level units specify otherwise.
    """

    unit_id: str = "su_id"
    material_id: str = "material_id"
    slope_angle_deg: str = "slope_angle_deg"
    maximum_depth_m: str = "maximum_depth_m"
    initial_water_table_depth_m: str = "initial_water_table_depth_m"
    depth_to_failure_plane_m: str | None = "failure_depth_m"
    long_term_infiltration_m_s: str = "long_term_infiltration_m_s"
    boundary: str | None = "boundary"
    drainage_efficiency: str | None = "drainage_efficiency"

    @property
    def required_columns(self) -> tuple[str, ...]:
        return (self.unit_id, self.material_id, self.slope_angle_deg,
                self.maximum_depth_m, self.initial_water_table_depth_m)


@dataclass
class SlopeUnitsResult:
    """Column results indexed by slope-unit identifier."""

    slope_units: pd.DataFrame
    schema: SlopeUnitSchema
    columns: dict[Hashable, ColumnResult]

    def summary(self) -> pd.DataFrame:
        """Return input attributes with deterministic stability diagnostics."""
        records = []
        for unit_id, result in self.columns.items():
            fs = result.arrays()["factor_of_safety"]
            time_index, depth_index = np.unravel_index(np.argmin(fs), fs.shape)
            plane_fs = result.factor_of_safety_at_failure_plane
            records.append(
                {
                    self.schema.unit_id: unit_id,
                    "minimum_factor_of_safety": float(fs[time_index, depth_index]),
                    "time_of_minimum_factor_of_safety_h": float(result.time[time_index]),
                    "depth_of_minimum_factor_of_safety_m": float(result.profiles[0]["z"][depth_index]),
                    "minimum_factor_of_safety_at_failure_plane": float(np.min(plane_fs)),
                    "failure_plane_failed": int(np.any(plane_fs <= 1.0)),
                    "maximum_water_table_rise_m": float(np.max(result.water_table_rise)),
                    "boundary_model": str(result.metadata["boundary"]),
                    "effective_initial_water_table_depth_m": float(result.metadata["initial_water_table_depth_m"]),
                    "initial_water_table_depth_was_clamped": bool(result.metadata["initial_water_table_depth_was_clamped"]),
                    "failure_plane_depth_m": float(result.metadata["failure_plane_depth_m"]),
                }
            )
        diagnostics = pd.DataFrame.from_records(records)
        return self.slope_units.merge(
            diagnostics, on=self.schema.unit_id, how="left", validate="one_to_one"
        )

    def result_for(self, unit_id: Hashable) -> ColumnResult:
        """Return the time/depth result for one slope unit."""
        return self.columns[unit_id]

    def to_geopackage(self, path: str | Path, *, layer: str = "unsatstab_results") -> None:
        """Write summary attributes and geometry to a GeoPackage."""
        _require_geopandas()
        summary = self.summary()
        if not hasattr(summary, "to_file"):
            raise TypeError("A GeoDataFrame input is required to write a GeoPackage.")
        summary.to_file(path, layer=layer, driver="GPKG")


def _require_geopandas():
    try:
        import geopandas as gpd
    except ImportError as error:  # pragma: no cover - dependency error is user-facing
        raise ImportError("Install geopandas to read or write slope-unit GeoPackages.") from error
    return gpd


def read_slope_units(path: str | Path, *, layer: str | None = None):
    """Read a GeoPackage layer with one record per slope unit."""
    gpd = _require_geopandas()
    return gpd.read_file(path, layer=layer)


def _finite_float(row: Mapping, column: str, *, unit_id: Hashable) -> float:
    try:
        value = float(row[column])
    except (TypeError, ValueError, KeyError) as error:
        raise ValueError(f"Slope unit {unit_id!r}: {column!r} must be numeric.") from error
    if not np.isfinite(value):
        raise ValueError(f"Slope unit {unit_id!r}: {column!r} must be finite.")
    return value


def _resolve_boundary(row: Mapping, schema: SlopeUnitSchema, default: BasalBoundary, unit_id):
    kind = default.name
    if schema.boundary and schema.boundary in row and pd.notna(row[schema.boundary]):
        kind = str(row[schema.boundary]).strip().lower()
    efficiency = default.drainage_efficiency
    if schema.drainage_efficiency and schema.drainage_efficiency in row and pd.notna(row[schema.drainage_efficiency]):
        efficiency = _finite_float(row, schema.drainage_efficiency, unit_id=unit_id)
    if kind in {"finite", "finite_perched"}:
        return FinitePerchedBoundary(drainage_efficiency=efficiency)
    if kind in {"infinite", "infinite_drainage"}:
        return InfiniteDrainageBoundary(drainage_efficiency=efficiency)
    raise ValueError(f"Slope unit {unit_id!r}: boundary must be 'finite' or 'infinite', not {kind!r}.")


# These optional row columns use the same units as the material table.
_HYDRO_OVERRIDES = {
    "saturated_hydraulic_conductivity_m_s": "saturated_hydraulic_conductivity",
    "alpha_m_inv": "alpha", "theta_s": "theta_s", "theta_r": "theta_r",
    "long_term_infiltration_m_s": "long_term_infiltration",
}
_STRENGTH_OVERRIDES = {
    "cohesion_pa": "cohesion", "friction_angle_deg": "friction_angle_deg",
    "unit_weight_soil_n_m3": "unit_weight_soil", "unit_weight_water_n_m3": "unit_weight_water",
}


def materials_from_dataframe(table: pd.DataFrame) -> dict:
    """Build material_id -> {hydro, strength, optional boundary} objects.

    Geometry is omitted from hydraulic templates; each SU supplies it. Required
    hydraulic columns are saturated_hydraulic_conductivity_m_s, alpha_m_inv,
    theta_s and theta_r; strength columns are cohesion_pa, friction_angle_deg
    and unit_weight_soil_n_m3. Background infiltration defaults to zero and
    water unit weight to 9810 N/m³. Optional boundary and drainage_efficiency
    columns configure a material-level boundary.
    """
    if not isinstance(table, pd.DataFrame) or table.empty:
        raise ValueError("material_properties must be a non-empty DataFrame.")
    required = {"material_id", "saturated_hydraulic_conductivity_m_s", "alpha_m_inv",
                "theta_s", "theta_r", "cohesion_pa", "friction_angle_deg", "unit_weight_soil_n_m3"}
    missing = required - set(table.columns)
    if missing:
        raise KeyError(f"Missing required material columns: {sorted(missing)}")
    if table.material_id.isna().any() or table.material_id.duplicated().any():
        raise ValueError("material_id must contain unique, non-null identifiers in the material table.")
    result = {}
    for row in table.to_dict(orient="records"):
        unit_id = row["material_id"]
        hydraulic = {field: _finite_float(row, column, unit_id=unit_id)
                     for column, field in _HYDRO_OVERRIDES.items()
                     if column in row and pd.notna(row[column])}
        strength = {field: _finite_float(row, column, unit_id=unit_id)
                    for column, field in _STRENGTH_OVERRIDES.items()
                    if column in row and pd.notna(row[column])}
        if "equation_set" in row and pd.notna(row["equation_set"]):
            hydraulic["equation_set"] = str(row["equation_set"])
        result[unit_id] = {"hydro": HydraulicProperties(**hydraulic),
                           "strength": StrengthProperties(**strength)}
        if any(key in row and pd.notna(row[key]) for key in ("boundary", "drainage_efficiency")):
            result[unit_id]["boundary"] = _resolve_boundary(row, SlopeUnitSchema(), FinitePerchedBoundary(), unit_id)
    return result


def _rainfall_pairs(columns) -> list[tuple[str, str]]:
    indices = {int(m.group(1)) for c in columns if (m := re.fullmatch(r"int_(\d+)", str(c)))}
    durations = {int(m.group(1)) for c in columns if (m := re.fullmatch(r"dur_(\d+)", str(c)))}
    if not indices and not durations:
        return []
    if indices != durations or indices != set(range(max(indices | durations) + 1)):
        raise ValueError("Rainfall columns must form consecutive int_0/dur_0, int_1/dur_1, ... pairs.")
    return [(f"int_{i}", f"dur_{i}") for i in sorted(indices)]


def _property_overrides(row, columns, aliases, unit_id):
    overrides = {field: _finite_float(row, column, unit_id=unit_id)
                 for column, field in columns.items()
                 if column in row and pd.notna(row[column])}
    for field in aliases:
        if field in row and pd.notna(row[field]):
            value = _finite_float(row, field, unit_id=unit_id)
            if field in overrides and value != overrides[field]:
                raise ValueError(f"Slope unit {unit_id!r}: conflicting aliases for {field!r}.")
            overrides[field] = value
    return overrides


def _simulate_one(row, material, *, rainfall, output_times, schema, default_boundary,
                  n_roots, n_depths, storage_time_step, antecedent, equation_set,
                  rainfall_pairs, rainfall_intensity_units, rainfall_duration_units):
    unit_id = row[schema.unit_id]
    overrides = _property_overrides(
        row, _HYDRO_OVERRIDES,
        ("saturated_hydraulic_conductivity", "alpha", "long_term_infiltration"), unit_id)
    # A custom schema may map background infiltration to another column.
    if schema.long_term_infiltration_m_s in row and pd.notna(row[schema.long_term_infiltration_m_s]):
        overrides["long_term_infiltration"] = _finite_float(row, schema.long_term_infiltration_m_s, unit_id=unit_id)
    overrides.update(
        slope_angle_deg=_finite_float(row, schema.slope_angle_deg, unit_id=unit_id),
        maximum_depth=_finite_float(row, schema.maximum_depth_m, unit_id=unit_id),
        initial_water_table_depth=_finite_float(row, schema.initial_water_table_depth_m, unit_id=unit_id),
    )
    if n_roots is not None:
        overrides["n_roots"] = n_roots
    if n_depths is not None:
        overrides["n_depths"] = n_depths
    if equation_set is not None:
        overrides["equation_set"] = equation_set
    hydro = replace(material["hydro"], **overrides)
    strength_overrides = _property_overrides(
        row, _STRENGTH_OVERRIDES, ("cohesion", "unit_weight_soil", "unit_weight_water"), unit_id)
    strength = replace(material["strength"], **strength_overrides)
    boundary = _resolve_boundary(row, schema, material.get("boundary", default_boundary), unit_id)
    forcing = material.get("rainfall", rainfall)
    values = [row[column] for pair in rainfall_pairs for column in pair]
    if values and any(pd.notna(v) for v in values):
        try:
            forcing = RainfallSeries.from_intervals(
                [row[i] for i, _ in rainfall_pairs], [row[d] for _, d in rainfall_pairs],
                intensity_units=rainfall_intensity_units, duration_units=rainfall_duration_units,
            )
        except (TypeError, ValueError) as error:
            raise ValueError(f"Slope unit {unit_id!r}: invalid rainfall intervals: {error}") from error
    if not isinstance(forcing, RainfallSeries):
        raise TypeError(f"Slope unit {unit_id!r}: supply rainfall in its row, material, or run argument.")
    plane = None
    if schema.depth_to_failure_plane_m and schema.depth_to_failure_plane_m in row and pd.notna(row[schema.depth_to_failure_plane_m]):
        plane = _finite_float(row, schema.depth_to_failure_plane_m, unit_id=unit_id)
    result = run_single_column(
        hydro=hydro, strength=strength, boundary=boundary, rainfall=forcing,
        output_times=output_times, depth_to_failure_plane=plane,
        storage_time_step=storage_time_step, antecedent=antecedent,
    )
    result.metadata.update(material_id=row[schema.material_id],
                           effective_hydraulic_properties={field: getattr(hydro, field) for field in _HYDRO_OVERRIDES.values()},
                           effective_strength_properties={field: getattr(strength, field) for field in _STRENGTH_OVERRIDES.values()})
    return unit_id, result


def run_slope_units(
    slope_units: pd.DataFrame,
    materials: Mapping | pd.DataFrame,
    *,
    output_times: np.ndarray,
    rainfall: RainfallSeries | None = None,
    schema: SlopeUnitSchema | None = None,
    column_mapping: Mapping[str, str] | None = None,
    default_boundary: BasalBoundary | None = None,
    n_roots: int | None = None,
    n_depths: int | None = None,
    n_jobs: int = -1,
    storage_time_step: float | None = None,
    antecedent: AntecedentRainfall | None = None,
    equation_set: str | None = None,
    rainfall_intensity_units: str = "mm/h",
    rainfall_duration_units: str = "h",
) -> SlopeUnitsResult:
    """Evaluate full-column responses for spatially variable slope units.

    Required SU columns are su_id, material_id, slope_angle_deg,
    maximum_depth_m, and initial_water_table_depth_m. column_mapping maps
    expected names to input names, e.g. {"su_id": "OBJECTID",
    "slope_angle_deg": "SLOPE"}. Unlisted names retain their defaults.
    Optional property and int_N/dur_N rainfall columns can also be mapped.
    Mapping changes names only, not units. Input attributes, geometry and
    column names are preserved in results and exports. Use either
    column_mapping or schema, not both.
    materials maps each ID to HydraulicProperties and StrengthProperties
    objects under hydro and strength. A material table is also accepted and
    converted with materials_from_dataframe. Geometry always comes from the SU.

    Optional non-null property columns override material defaults. Rainfall
    precedence is SU int_N/dur_N pairs, material rainfall, then run rainfall.
    Defaults for the pairs are mm/h and hours; durations are accumulated into
    breakpoints automatically. A row with all rainfall-pair values missing
    inherits forcing; a partly specified history is rejected.

    Optional failure_depth_m selects the pressure/FS time-series plane; omission
    selects the maximum depth. Full hydraulic and stability profiles are always
    returned. n_roots, n_depths and equation_set optionally override every
    material's corresponding hydraulic setting; otherwise templates retain
    their values. Output times (hours) and storage_time_step follow run_single_column.
    """
    if not isinstance(slope_units, pd.DataFrame) or slope_units.empty:
        raise ValueError("slope_units must be a non-empty pandas or GeoPandas DataFrame.")
    if rainfall is not None and not isinstance(rainfall, RainfallSeries):
        raise TypeError("rainfall must be a RainfallSeries or None.")
    if default_boundary is not None and not isinstance(default_boundary, BasalBoundary):
        raise TypeError("default_boundary must be a BasalBoundary or None.")
    original_units = slope_units
    if column_mapping is not None:
        if schema is not None:
            raise ValueError("Use either column_mapping or schema, not both.")
        if not isinstance(column_mapping, Mapping):
            raise TypeError("column_mapping must map expected names to input column names.")
        defaults = SlopeUnitSchema()
        allowed = set(vars(defaults).values()) | set(_HYDRO_OVERRIDES) | set(_STRENGTH_OVERRIDES)
        allowed |= {"saturated_hydraulic_conductivity", "alpha", "long_term_infiltration",
                    "cohesion", "unit_weight_soil", "unit_weight_water"}
        for expected, source in column_mapping.items():
            if not isinstance(expected, str) or (expected not in allowed and not re.fullmatch(r"(?:int|dur)_\d+", expected)):
                raise ValueError(f"Unknown expected column name: {expected!r}.")
            if not isinstance(source, str) or not source:
                raise ValueError("Mapped input names must be non-empty strings.")
            if source not in slope_units.columns:
                raise KeyError(f"Mapped input column not found: {source!r}.")
        if len(set(column_mapping.values())) != len(column_mapping):
            raise ValueError("Each mapped input column must represent only one expected column.")
        rename_columns = {source: expected for expected, source in column_mapping.items()}
        renamed = [rename_columns.get(c, c) for c in slope_units.columns]
        if len(set(renamed)) != len(renamed):
            raise ValueError("column_mapping creates duplicate expected column names.")
        slope_units = slope_units.rename(columns=rename_columns)
        result_schema = replace(defaults, **{
            field: column_mapping.get(name, name) for field, name in vars(defaults).items()
        })
    else:
        result_schema = schema or SlopeUnitSchema()
    schema = schema or SlopeUnitSchema()
    missing = sorted(set(schema.required_columns) - set(slope_units.columns))
    if missing:
        raise KeyError(f"Missing required slope-unit columns: {missing}")
    if slope_units[schema.unit_id].isna().any() or slope_units[schema.unit_id].duplicated().any():
        raise ValueError(f"{schema.unit_id!r} must contain unique, non-null identifiers.")
    if slope_units[schema.material_id].isna().any():
        raise ValueError("material_id must be non-null in every slope unit.")
    if isinstance(materials, pd.DataFrame):
        materials = materials_from_dataframe(materials)
    if not isinstance(materials, Mapping) or not materials:
        raise ValueError("materials must be a non-empty mapping or material table.")
    for material_id, material in materials.items():
        if (not isinstance(material, Mapping) or not isinstance(material.get("hydro"), HydraulicProperties)
                or not isinstance(material.get("strength"), StrengthProperties)):
            raise TypeError(f"Material {material_id!r} requires hydro and strength parameter objects.")
        if "boundary" in material and not isinstance(material["boundary"], BasalBoundary):
            raise TypeError(f"Material {material_id!r}: boundary must be a BasalBoundary.")
        if "rainfall" in material and not isinstance(material["rainfall"], RainfallSeries):
            raise TypeError(f"Material {material_id!r}: rainfall must be a RainfallSeries.")
    missing_ids = set(slope_units[schema.material_id]) - set(materials)
    if missing_ids:
        raise KeyError(f"No material properties found for material IDs: {sorted(missing_ids, key=str)}")
    for name, value in (("n_roots", n_roots), ("n_depths", n_depths)):
        if value is not None and (not isinstance(value, int) or value <= 0):
            raise ValueError(f"{name} must be a positive integer.")
    if not isinstance(n_jobs, int) or n_jobs == 0:
        raise ValueError("n_jobs must be a non-zero integer.")
    if equation_set is not None and equation_set not in ("baum_etal", "srivastava_yeh"):
        raise ValueError("equation_set must be 'baum_etal' or 'srivastava_yeh'.")
    rainfall_pairs = _rainfall_pairs(slope_units.columns)
    jobs = [delayed(_simulate_one)(
        row, materials[row[schema.material_id]], rainfall=rainfall, output_times=output_times,
        schema=schema, default_boundary=default_boundary or FinitePerchedBoundary(),
        n_roots=n_roots, n_depths=n_depths, equation_set=equation_set,
        storage_time_step=storage_time_step, antecedent=antecedent,
        rainfall_pairs=rainfall_pairs, rainfall_intensity_units=rainfall_intensity_units,
        rainfall_duration_units=rainfall_duration_units,
    ) for row in slope_units.to_dict(orient="records")]
    columns = dict(Parallel(n_jobs=n_jobs, prefer="processes")(jobs))
    return SlopeUnitsResult(slope_units=original_units.copy(), schema=result_schema, columns=columns)


def run_slope_units_from_gpkg(
    path: str | Path, materials: Mapping | pd.DataFrame, *,
    output_times: np.ndarray, rainfall: RainfallSeries | None = None,
    layer: str | None = None, schema: SlopeUnitSchema | None = None,
    column_mapping: Mapping[str, str] | None = None,
    default_boundary: BasalBoundary | None = None,
    n_roots: int | None = None, n_depths: int | None = None, n_jobs: int = -1,
    storage_time_step: float | None = None, antecedent: AntecedentRainfall | None = None,
    equation_set: str | None = None, rainfall_intensity_units: str = "mm/h",
    rainfall_duration_units: str = "h",
) -> SlopeUnitsResult:
    """Read a GeoPackage and run the material/override/forcing workflow.

    Parameters and precedence are the same as run_slope_units; layer selects
    the GeoPackage layer. column_mapping maps expected names to layer column
    names and is forwarded unchanged; schema is an alternative.
    """
    return run_slope_units(
        read_slope_units(path, layer=layer), materials, output_times=output_times,
        rainfall=rainfall, schema=schema, column_mapping=column_mapping,
        default_boundary=default_boundary,
        n_roots=n_roots, n_depths=n_depths, n_jobs=n_jobs,
        storage_time_step=storage_time_step, antecedent=antecedent,
        equation_set=equation_set, rainfall_intensity_units=rainfall_intensity_units,
        rainfall_duration_units=rainfall_duration_units,
    )
