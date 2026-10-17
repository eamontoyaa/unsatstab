import geopandas as gpd
import numpy as np
import pandas as pd
import pytest
from shapely.geometry import Polygon

import unsatstab as us


def _slope_units():
    return pd.DataFrame({
        "su_id": ["below_wt", "above_wt"],
        "material_id": ["volcanic_soil", "volcanic_soil"],
        "slope_angle_deg": [40.0, 45.0],
        "maximum_depth_m": [1.01, 0.83],
        "initial_water_table_depth_m": [0.75, 1.2],
        "long_term_infiltration_m_s": [1.5e-7, 1.5e-7],
        "boundary": ["finite", "infinite"],
    })


def _materials():
    return pd.DataFrame({
        "material_id": ["volcanic_soil"],
        "saturated_hydraulic_conductivity_m_s": [1.0e-5],
        "alpha_m_inv": [3.0], "theta_s": [0.41], "theta_r": [0.06],
        "cohesion_pa": [3000.0], "friction_angle_deg": [33.6],
        "unit_weight_soil_n_m3": [18000.0],
    })


def test_slope_units_reuse_materials_and_support_parallel_execution():
    result = us.run_slope_units(
        _slope_units(), _materials(),
        rainfall=us.RainfallSeries.constant(4.5, 24.0),
        output_times=(np.arange(5.0) * 3600.0) / 3600.0,
        n_roots=40, n_depths=30, n_jobs=2,
    )
    summary = result.summary()
    assert set(result.columns) == {"below_wt", "above_wt"}
    assert set(summary["su_id"]) == {"below_wt", "above_wt"}
    assert set(summary["boundary_model"]) == {"finite_perched", "infinite_drainage"}
    assert np.all(np.isfinite(summary["minimum_factor_of_safety"]))
    assert result.result_for("below_wt").metadata["initial_water_table_depth_m"] == 0.75
    assert result.result_for("above_wt").metadata["initial_water_table_depth_m"] == 1.2


@pytest.mark.parametrize('equation_set',['baum_etal','srivastava_yeh'])
def test_equation_set_reaches_spatial_columns(equation_set):
    result=us.run_slope_units(_slope_units(),_materials(),
        rainfall=us.RainfallSeries.constant(4.5,24.),output_times=(np.array([0.])) / 3600.0,
        n_roots=20,n_depths=20,n_jobs=1,equation_set=equation_set)
    assert all(c.metadata['equation_set']==equation_set for c in result.columns.values())


def test_slope_units_reject_missing_material_and_duplicate_unit_ids():
    with pytest.raises(KeyError, match="material IDs"):
        us.run_slope_units(
            _slope_units().assign(material_id="missing"), _materials(),
            rainfall=us.RainfallSeries.constant(4.5, 24.0), output_times=(np.arange(2.0) * 3600.0) / 3600.0,
            n_jobs=1,
        )
    duplicated = _slope_units()
    duplicated.loc[1, "su_id"] = "below_wt"
    with pytest.raises(ValueError, match="unique"):
        us.run_slope_units(
            duplicated, _materials(), rainfall=us.RainfallSeries.constant(4.5, 24.0),
            output_times=(np.arange(2.0) * 3600.0) / 3600.0, n_jobs=1,
        )


def test_geopackage_round_trip_preserves_geometry_and_summary(tmp_path):
    units = gpd.GeoDataFrame(
        _slope_units().iloc[[0]].copy(),
        geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])],
        crs="EPSG:4326",
    )
    result = us.run_slope_units(
        units, _materials(), rainfall=us.RainfallSeries.constant(4.5, 24.0),
        output_times=(np.arange(3.0) * 3600.0) / 3600.0, n_roots=20, n_depths=20, n_jobs=1,
    )
    path = tmp_path / "results.gpkg"
    result.to_geopackage(path)
    recovered = us.read_slope_units(path, layer="unsatstab_results")
    assert recovered.crs == units.crs
    assert recovered.loc[0, "su_id"] == "below_wt"
    assert "minimum_factor_of_safety" in recovered
    simulated = us.run_slope_units_from_gpkg(
        path,_materials(),layer="unsatstab_results",
        rainfall=us.RainfallSeries.constant(4.5,24),output_times=(np.array([0.,3600.])) / 3600.0,
        antecedent=us.AntecedentRainfall.uniform(24,(1) * 24.0,drainage_efficiency=1.),
        n_roots=20,n_depths=20,n_jobs=1,
    ).result_for("below_wt")
    assert simulated.metadata["antecedent_rainfall_mm"] == pytest.approx(24.)
    assert simulated.metadata["antecedent_drainage_efficiency"] == 1.
    np.testing.assert_array_equal((simulated.time) * 3600.0,np.array([0.,3600.]))


def test_spatial_antecedent_history_uses_independent_quadrature():
    units = _slope_units().iloc[[0]].copy()
    history = us.RainfallSeries.constant(4.5,24).with_antecedent_rainfall(500,(90) * 24.0)
    times = np.array([0.,90*86400.,90*86400.+24*3600.])
    result = us.run_slope_units(
        units,_materials(),rainfall=history,output_times=(times) / 3600.0,
        n_roots=40,n_depths=30,n_jobs=1,storage_time_step=(3600.) / 3600.0,
    ).result_for("below_wt")
    direct = us.run_single_column(
        hydro=us.HydraulicProperties(slope_angle_deg=40,maximum_depth=1.01,saturated_hydraulic_conductivity=1e-5,alpha=3,theta_s=.41,theta_r=.06,long_term_infiltration=1.5e-7,
            initial_water_table_depth=.75,n_roots=40,n_depths=30),
        strength=us.StrengthProperties(3000,33.6,18000),
        boundary=us.FinitePerchedBoundary(),rainfall=history,output_times=(times) / 3600.0,
        storage_time_step=(3600.) / 3600.0,
    )
    assert len(result.profiles) == 3
    np.testing.assert_allclose(result.cumulative_excess,direct.cumulative_excess)
    np.testing.assert_allclose(result.factor_of_safety_at_failure_plane,direct.factor_of_safety_at_failure_plane)


def test_spatial_two_stage_drainage():
    units = _slope_units().iloc[[0]].copy()
    units["long_term_infiltration_m_s"] = 0.
    result = us.run_slope_units(
        units,_materials(),rainfall=us.RainfallSeries.constant(4.5,24),
        output_times=(np.array([0.,24*3600.])) / 3600.0,
        default_boundary=us.FinitePerchedBoundary(drainage_efficiency=.015),
        n_roots=40,n_depths=30,n_jobs=1,storage_time_step=(3600) / 3600.0,
        antecedent=us.AntecedentRainfall.uniform(1000,(90) * 24.0,drainage_efficiency=1.),
    ).result_for("below_wt")
    assert result.water_table_depth[0] == .75
    assert result.cumulative_excess[0] == 0
    assert result.metadata["antecedent_profile_water_change_mm"] > 0
    assert result.cumulative_excess[-1] > 0
    np.testing.assert_array_equal((result.time) * 3600.0,np.array([0.,86400.]))


def test_shared_antecedent_uses_each_units_capacity_and_ks_in_parallel():
    units = _slope_units()
    units.loc[1,"material_id"] = "less_permeable"
    units["drainage_efficiency"] = [.015,1.]
    materials = pd.concat([
        _materials(),
        _materials().assign(material_id="less_permeable",saturated_hydraulic_conductivity_m_s=1e-6),
    ],ignore_index=True)
    previous = us.RainfallSeries.from_intensities(np.array([0.,50.,0.]),6,time_units="h")
    antecedent = us.AntecedentRainfall(previous)
    rainfall = us.RainfallSeries.constant(4.5,24)
    times = np.array([0.,3600.])
    spatial = us.run_slope_units(
        units,materials,rainfall=rainfall,output_times=(times) / 3600.0,antecedent=antecedent,
        n_roots=40,n_depths=30,n_jobs=2,storage_time_step=(1800.) / 3600.0,
    )
    for i, row in units.iterrows():
        material = materials.iloc[i]
        boundary = us.FinitePerchedBoundary(drainage_efficiency=.015) if i==0 else us.InfiniteDrainageBoundary()
        direct = us.run_single_column(
            hydro=us.HydraulicProperties(
                slope_angle_deg=row.slope_angle_deg,maximum_depth=row.maximum_depth_m,
                saturated_hydraulic_conductivity=material.saturated_hydraulic_conductivity_m_s,alpha=3,theta_s=.41,theta_r=.06,long_term_infiltration=1.5e-7,
                initial_water_table_depth=row.initial_water_table_depth_m,n_roots=40,n_depths=30,
            ),
            rainfall=rainfall,strength=us.StrengthProperties(3000,33.6,18000),boundary=boundary,
            output_times=(times) / 3600.0,antecedent=antecedent,storage_time_step=(1800.) / 3600.0,
        )
        result = spatial.result_for(row.su_id)
        np.testing.assert_allclose(result.factor_of_safety_at_failure_plane,direct.factor_of_safety_at_failure_plane)
        np.testing.assert_allclose(result.cumulative_excess,direct.cumulative_excess)
        assert result.metadata["antecedent_drainage_efficiency"] == boundary.drainage_efficiency
        assert result.metadata["antecedent_effective_infiltration_mm"] == pytest.approx(216 if i==0 else 21.6)


def _object_materials():
    return {"volcanic_soil": {
        "hydro": us.HydraulicProperties(saturated_hydraulic_conductivity=1e-5,
            alpha=3, theta_s=.41, theta_r=.06, long_term_infiltration=1.5e-7,
            n_roots=30, n_depths=30),
        "strength": us.StrengthProperties(3000, 33.6, 18000),
        "boundary": us.FinitePerchedBoundary(drainage_efficiency=.015),
        "rainfall": us.RainfallSeries.constant(4.5, 24),
    }}


def test_five_required_spatial_columns_inherit_typed_materials_without_mutating_them():
    units = _slope_units()[list(us.SlopeUnitSchema().required_columns)]
    materials = _object_materials()
    result = us.run_slope_units(units, materials, output_times=(np.array([0, 3600])) / 3600.0, n_jobs=2)
    for _, row in units.iterrows():
        r = result.result_for(row.su_id)
        assert r.arrays()["z"][-1] == row.maximum_depth_m
        assert r.metadata["initial_water_table_depth_m"] == min(row.initial_water_table_depth_m, row.maximum_depth_m)
        assert r.metadata["effective_hydraulic_properties"]["long_term_infiltration"] == 1.5e-7
    assert materials["volcanic_soil"]["hydro"].maximum_depth is None
    assert materials["volcanic_soil"]["hydro"].slope_angle_deg is None


def test_local_properties_and_intensity_duration_pairs_override_material_and_shared_rain():
    units = _slope_units()[list(us.SlopeUnitSchema().required_columns)].copy()
    units["alpha"] = [2., np.nan]
    units["cohesion"] = [1000., np.nan]
    units["failure_depth_m"] = [.5, .5]
    units["int_0"] = [0., np.nan]; units["dur_0"] = [1., np.nan]
    units["int_1"] = [10., np.nan]; units["dur_1"] = [2., np.nan]
    m = _object_materials()
    result = us.run_slope_units(units, m, rainfall=us.RainfallSeries.constant(1, 3),
                               output_times=(np.array([0, 3600, 3*3600])) / 3600.0, n_jobs=1)
    local = result.result_for("below_wt")
    np.testing.assert_allclose(local.rainfall_flux*3.6e6, [0, 10, 0])
    assert local.metadata["effective_hydraulic_properties"]["alpha"] == 2
    assert local.metadata["effective_strength_properties"]["cohesion"] == 1000
    assert local.metadata["failure_plane_depth_m"] == .5
    inherited = result.result_for("above_wt")
    assert inherited.rainfall_flux[0]*3.6e6 == pytest.approx(4.5)
    assert inherited.metadata["effective_hydraulic_properties"]["alpha"] == 3
    assert m["volcanic_soil"]["hydro"].alpha == 3
    broken = units.copy(); broken.loc[0, "dur_1"] = np.nan
    with pytest.raises(ValueError, match="below_wt"):
        us.run_slope_units(broken, m, output_times=(np.array([0])) / 3600.0, n_jobs=1)
    for column in ("int_1", "dur_1"):
        with pytest.raises(ValueError, match="consecutive"):
            us.run_slope_units(units.drop(columns=column), m, output_times=(np.array([0])) / 3600.0, n_jobs=1)


@pytest.mark.parametrize("column,property_name,value,object_key", [
    ("saturated_hydraulic_conductivity_m_s", "saturated_hydraulic_conductivity", 2e-5, "effective_hydraulic_properties"),
    ("theta_s", "theta_s", .45, "effective_hydraulic_properties"),
    ("theta_r", "theta_r", .05, "effective_hydraulic_properties"),
    ("long_term_infiltration_m_s", "long_term_infiltration", 0., "effective_hydraulic_properties"),
    ("friction_angle_deg", "friction_angle_deg", 32., "effective_strength_properties"),
    ("unit_weight_soil", "unit_weight_soil", 20000., "effective_strength_properties"),
])
def test_optional_soil_columns_override_inherited_properties(column, property_name, value, object_key):
    units = _slope_units()[list(us.SlopeUnitSchema().required_columns)].iloc[:1].copy()
    units[column] = value
    result = us.run_slope_units(units, _object_materials(), output_times=(np.array([0.])) / 3600.0, n_jobs=1)
    assert result.result_for("below_wt").metadata[object_key][property_name] == value


def test_conflicting_optional_aliases_are_rejected():
    units = _slope_units()[list(us.SlopeUnitSchema().required_columns)].iloc[:1].copy()
    units["alpha"] = 2.; units["alpha_m_inv"] = 3.
    with pytest.raises(ValueError, match="conflicting aliases"):
        us.run_slope_units(units, _object_materials(), output_times=(np.array([0.])) / 3600.0, n_jobs=1)


def test_required_material_identifier_cannot_be_missing():
    units = _slope_units()[list(us.SlopeUnitSchema().required_columns)].iloc[:1].copy()
    units["material_id"] = None
    with pytest.raises(ValueError, match="material_id must be non-null"):
        us.run_slope_units(units, {None: _object_materials()["volcanic_soil"]},
                           output_times=(np.array([0.])) / 3600.0, n_jobs=1)


def test_gis_mapping_matches_canonical_inputs_and_preserves_export(tmp_path):
    units = gpd.GeoDataFrame(_slope_units().iloc[:1].copy(),
        geometry=[Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])], crs="EPSG:4326")
    units["cohesion_pa"] = 1200.
    units["int_0"], units["dur_0"] = 3., 2.
    mapping = {"su_id": "OBJECTID", "material_id": "MATERIAL",
               "slope_angle_deg": "SLOPE", "maximum_depth_m": "DEPTH",
               "initial_water_table_depth_m": "WT_DEPTH", "cohesion_pa": "COHESION",
               "int_0": "RAIN_1", "dur_0": "HOURS_1"}
    gis = units.rename(columns=mapping)
    before = gis.copy(deep=True)
    kwargs = dict(output_times=np.array([0., 1., 2.]), n_jobs=1, n_roots=20, n_depths=20)
    expected = us.run_slope_units(units, _materials(), **kwargs)
    path = tmp_path / "gis.gpkg"
    gis.to_file(path, layer="units", driver="GPKG")
    actual = us.run_slope_units_from_gpkg(path, _materials(), layer="units",
                                        column_mapping=mapping, **kwargs)
    np.testing.assert_allclose(actual.result_for("below_wt").arrays()["factor_of_safety"],
                               expected.result_for("below_wt").arrays()["factor_of_safety"])
    pd.testing.assert_frame_equal(gis, before)
    assert actual.slope_units.crs == gis.crs
    assert set(actual.slope_units.columns) == set(gis.columns)
    pd.testing.assert_frame_equal(actual.slope_units[gis.columns], gis, check_dtype=False)
    assert actual.summary().loc[0, "OBJECTID"] == "below_wt"
    output = tmp_path / "output.gpkg"
    actual.to_geopackage(output)
    recovered = us.read_slope_units(output)
    assert recovered.loc[0, "OBJECTID"] == "below_wt"
    assert "COHESION" in recovered and "cohesion_pa" not in recovered


@pytest.mark.parametrize("mapping,error,match", [
    ({"unknown": "su_id"}, ValueError, "Unknown expected"),
    ({"su_id": "MISSING"}, KeyError, "not found"),
    ({"su_id": "material_id", "material_id": "material_id"}, ValueError, "only one"),
    ({"su_id": "material_id"}, ValueError, "duplicate expected"),
    ({"su_id": None}, ValueError, "non-empty strings"),
])
def test_invalid_gis_mappings_are_rejected(mapping, error, match):
    with pytest.raises(error, match=match):
        us.run_slope_units(_slope_units(), _materials(), output_times=np.array([0.]),
                           column_mapping=mapping, n_jobs=1)


def test_gis_mapping_and_schema_are_alternatives():
    with pytest.raises(ValueError, match="either column_mapping or schema"):
        us.run_slope_units(_slope_units(), _materials(), output_times=np.array([0.]),
                           column_mapping={}, schema=us.SlopeUnitSchema(), n_jobs=1)
