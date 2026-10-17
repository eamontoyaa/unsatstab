# Inputs and outputs

Use this guide to configure a column or a slope-unit domain and read the results. Public times are in **hours**; geometry and material properties use SI units unless a constructor specifies otherwise. See the [tutorials](index.md#tutorials) for complete examples and the [API reference](api.md) for function signatures.

## Single-column parameter conventions

| Input | Meaning | Unit |
| --- | --- | --- |
| `slope_angle_deg` | Ground-surface slope | degrees |
| `maximum_depth` | Maximum vertical column depth | m |
| `initial_water_table_depth` | Initial water-table depth; defaults to the column base | m |
| `saturated_hydraulic_conductivity` | Saturated conductivity, $K_{\mathrm{s}}$ | m/s |
| `alpha` | Intrinsic Gardner parameter | 1/m |
| `theta_s`, `theta_r` | Saturated and residual water contents | m³/m³ |
| `long_term_infiltration` | Background flux defining the steady initial profile | m/s |
| `cohesion`, `friction_angle_deg` | Effective cohesion and friction angle | Pa, degrees |
| `unit_weight_soil`, `unit_weight_water` | Soil and water unit weights | N/m³ |
| `equation_set` | `"srivastava_yeh"` (default) or `"baum_etal"` | — |

Configure hydraulic parameters with `HydraulicProperties`, strength with `StrengthProperties`, and basal drainage with a [boundary object](api/boundaries.md). `srivastava_yeh` uses the intrinsic Gardner laws; `baum_etal` uses the TRIGRS-style conversions. Their differences are described in the [equation comparison](fundamentals.md#pressure-and-water-content-equation-sets).

### Column depth and observation plane

`maximum_depth` sets the reported profile extent. With a finite basal boundary, an initial water table below this depth is adjusted to the base; an infinite boundary retains the requested depth. Results record both requested and effective initial depths.

The optional `depth_to_failure_plane` argument selects the pressure and factor-of-safety time-series plane within the column. It defaults to the base and does not change the full profiles or hydraulic calculation.

```python
result = us.run_single_column(
    hydro=hydro, strength=strength, boundary=boundary, rainfall=rainfall,
    output_times=output_times_h,
    depth_to_failure_plane=0.8,  # vertical depth in metres
)
```

The analytical unsaturated domain extends to the effective initial water table. Water-table rise and pressure below it follow the [storage and reconstruction assumptions](fundamentals.md#basal-and-pressure-assumptions).

## Public time convention

Use hours for `output_times`, rainfall breakpoints, direct analytical evaluations, and `storage_time_step`. Returned times also use hours; the solver converts to SI seconds internally.

Rainfall constructors accept alternative units explicitly, for example `duration_units="days"`. Stored histories still use hours. `output_times` must start at zero and increase strictly; its last value sets the event simulation horizon.

## Rainfall forcing

`RainfallSeries` describes piecewise-constant gross rainfall. Convenience constructors default to **mm/h** and **hours**:

```python
rainfall = us.RainfallSeries.from_intervals(
    intensities=[0.0, 9.0, 3.0],  # mm/h
    durations=[3.0, 2.0, 1.0],    # hours
)
```

Include zero-intensity intervals for dry gaps. Each column admits $I(t)=\min[P(t),K_{\mathrm{s}}]$; rainfall above its local conductivity is discarded without surface storage or routing.

Intervals use `[start, end)`: the next rate applies at an interior breakpoint, and rainfall is zero at the final breakpoint. Plot instantaneous rates with `step(..., where="post")`.

## Antecedent rainfall

Pass an `AntecedentRainfall` object to simulate rainfall before the triggering event. Omit it for an event without antecedent forcing.

```python
antecedent = us.AntecedentRainfall.uniform(
    total_rainfall_mm=500,
    duration=90, duration_units="days",  # stored history uses hours
    drainage_efficiency=1.0,
)
result = us.run_single_column(
    hydro=hydro, strength=strength, boundary=boundary, rainfall=rainfall,
    output_times=output_times_h, antecedent=antecedent,
    storage_time_step=1.0,  # hours
)
```

For a measured history, use `AntecedentRainfall(series, drainage_efficiency=...)` with a `RainfallSeries` beginning at zero. The same object is accepted by spatial solvers. Each antecedent interval receives the local conductivity cap.

`drainage_efficiency=None` uses the event boundary's setting throughout. An explicit value in [0, 1] applies during the antecedent only. Efficiency multiplies drainage capacity, $q_{\mathrm{c}}=\eta K_{\mathrm{s}}\sin^2\delta$; a value of 1 does not impose a fixed water table.

Moisture and retained storage continue into the event. Returned times are event-relative, while water-table rise and cumulative storage reference the start of the complete history and can be positive at event zero. See the [antecedent tutorial](tutorials/02_antecedent_wetness.ipynb) and [physical assumptions](fundamentals.md#finite-antecedent-infiltration-history).

### Antecedent diagnostic fields

Metadata includes rainfall amounts, transfer, drainage, storage, and profile-water changes. All volume fields below are in millimetres of water; antecedent `flux_m_s` fields are period means.

| Metadata field | Meaning |
| --- | --- |
| `antecedent_effective_infiltration_mm` | Rainfall admitted after the local cap |
| `antecedent_basal_transfer_mm` | Integrated flux arriving at the base |
| `antecedent_basal_storage_mm` | Retained basal water |
| `antecedent_drainage_mm` | Drainage allocated by the storage balance |
| `antecedent_capacity_drainage_mm` | Integral of basal flux limited by drainage capacity |
| `antecedent_background_correction_mm` | Allocated minus capacity-limited drainage |
| `antecedent_profile_water_change_mm` | Reconstructed profile-water change |
| `antecedent_unsaturated_response_water_change_mm` | Profile-water change before reconstruction |
| `antecedent_water_balance_residual_mm` | Effective input minus allocated drainage and reconstructed profile-water change |
| `antecedent_storage_above_pore_capacity_mm` | Retained storage exceeding fillable pore volume |

Profile changes reference the start of the complete history. Basal storage is already represented in the reconstructed profile; do not add them as independent stores. See [water-balance checks](verification.md#antecedent-water-balance-checks).

## Slope-unit domain and materials {#spatial-tables}

`run_slope_units` represents each slope unit (SU) as an independent column. Supply a table or GeoDataFrame with these required fields:

| Column | Meaning | Unit |
| --- | --- | --- |
| `su_id` | Unique slope-unit identifier | — |
| `material_id` | Key in the material mapping | — |
| `slope_angle_deg` | Local slope | degrees |
| `maximum_depth_m` | Maximum vertical column depth | m |
| `initial_water_table_depth_m` | Requested initial water-table depth | m |

Optional `failure_depth_m` selects the observation plane; omission selects the base. Geometry comes from each SU. Material objects supply shared hydraulic and strength properties:

```python
materials = {
    "soil_a": {
        "hydro": us.HydraulicProperties(
            saturated_hydraulic_conductivity=1e-5,
            alpha=3.0, theta_s=0.41, theta_r=0.06,
            long_term_infiltration=1.5e-7,
        ),
        "strength": us.StrengthProperties(
            cohesion=3000.0, friction_angle_deg=33.6, unit_weight_soil=18000.0,
        ),
        "boundary": us.FinitePerchedBoundary(drainage_efficiency=0.1),
    },
}
```

Each material requires `hydro` and `strength`; `boundary` and `rainfall` are optional defaults. Each SU receives its own resolved configuration. Material DataFrames or CSV files are also supported; see [material-table fields](api/slope_units.md).

### Spatially varying property overrides

Non-null SU properties override the mapped material; missing values inherit it.

| SU column | Property |
| --- | --- |
| `saturated_hydraulic_conductivity_m_s` | Saturated conductivity |
| `alpha_m_inv` | Gardner parameter |
| `theta_s`, `theta_r` | Water-content limits |
| `long_term_infiltration_m_s` | Background flux |
| `cohesion_pa`, `friction_angle_deg` | Strength |
| `unit_weight_soil_n_m3`, `unit_weight_water_n_m3` | Unit weights |
| `boundary`, `drainage_efficiency` | Basal drainage |

Constructor names such as `alpha` and `cohesion` are also accepted in SI units. Use one name per property; conflicting aliases are rejected. Run-level `n_roots`, `n_depths`, and `equation_set` apply to every unit when supplied; otherwise material settings apply.

### Spatial rainfall histories

SU columns `int_0`, `dur_0`, `int_1`, `dur_1`, … define consecutive rainfall intervals in **mm/h** and **hours**. Indices start at zero without gaps; each intensity needs a positive duration. Use zero intensity for dry intervals. Alternative units are set with `rainfall_intensity_units` and `rainfall_duration_units`.

Rainfall precedence is **SU intervals → material rainfall → shared run rainfall**. Rows with no interval values inherit rainfall; partial histories are rejected. Output times are shared across units, and each history supplies zero rainfall after it ends.

```python
units["int_0"], units["dur_0"] = 0.0, 3.0
units["int_1"], units["dur_1"] = 9.0, 2.0
result = us.run_slope_units(
    units, materials, output_times=output_times_h, n_jobs=2,
)
```

### Column names from GIS files

Use `column_mapping` to map **expected names to input-file names**. List only fields needing a mapping:

```python
column_mapping = {
    "su_id": "OBJECTID",
    "material_id": "MATERIAL",
    "slope_angle_deg": "SLOPE",
    "maximum_depth_m": "DEPTH",
    "initial_water_table_depth_m": "WT_DEPTH",
}
result = us.run_slope_units(
    slope_units, materials, rainfall=rainfall,
    output_times=output_times_h, column_mapping=column_mapping,
)
```

`run_slope_units_from_gpkg` accepts the same mapping, including optional properties and rainfall fields. Mapping changes names only; values must use the documented units. `SlopeUnitSchema` offers an alternative for geometry and boundary fields; supply either `schema` or `column_mapping`.

Input attributes, geometry, and CRS are retained in spatial results and exports. `SlopeUnitsResult.summary()` reports minimum FS over the full column and at the selected plane. Use `result_for(su_id)` for a column result and `to_geopackage()` for a spatial export.

## Column outputs

| Accessor | Content |
| --- | --- |
| `time_dataframe()` | Time series at the selected plane, fluxes, storage, and water-table rise |
| `profile_dataframe()` | Depth profiles for each output time |
| `arrays()` | NumPy arrays |
| `metadata` | Configuration, effective initial water table, and diagnostics |

`rainfall_flux` is gross rainfall, `infiltration_flux` is the admitted rate, and `rainfall_excess_flux` is the rejected difference. `basal_flux` is analytical transfer at the effective initial water table; `excess_flux` is undrained basal transfer. See [result fields](api/results.md).

The selected-plane FS differs from the minimum over the whole column. $\mathrm{FS}\le1$ marks the model's deterministic failure threshold, rather than a failure probability.

## Numerical controls

| Control | Default | Role |
| --- | --- | --- |
| `n_roots` | 250 | Analytical series terms; `srivastava_yeh` resolves extra terms at short times when needed |
| `n_depths` | 250 | Depth intervals for profiles and storage calculations |
| `small_time_threshold` | 1e-4 | Dimensionless cutoff for the `baum_etal` early-time approximation |
| `output_times` | Required | Saved time samples and initial storage integration grid |
| `storage_time_step` | None; 1 h with antecedent forcing | Maximum internal storage integration interval, in hours |
| `max_factor_of_safety` | 10 | Upper FS output cap |
| `n_jobs` | -1 | Spatial worker count; 1 runs serially |

### Output sampling and storage integration

`output_times` controls saved results. `storage_time_step` refines internal storage integration independently: hourly outputs with `storage_time_step=5 / 60` use intervals of at most five minutes while saving hourly profiles. Rainfall breakpoints and drainage transitions are included automatically.

Without an explicit `storage_time_step`, event-only simulations use output times and rainfall breakpoints; simulations with antecedent forcing use a maximum interval of 1 h. Finer integration can change storage, water-table rise, and FS. Finer output sampling is needed to resolve the first sampled failure time.

Compare root counts, depth resolutions, and storage intervals until the quantities of interest stabilize. Numerical refinement does not remove the model's [water-table reconstruction assumptions](fundamentals.md#basal-and-pressure-assumptions). Details of the series and early-time response are in [analytical infiltration](fundamentals.md#evaluating-the-response-at-short-times).
