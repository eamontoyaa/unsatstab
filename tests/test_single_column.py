import numpy as np
import pytest
import unsatstab as us


def test_single_column_two_boundaries():
    hydro = us.HydraulicProperties(
        slope_angle_deg=40.0,
        maximum_depth=1.01,
        saturated_hydraulic_conductivity=1.0e-5,
        alpha=3.0,
        theta_s=0.41,
        theta_r=0.06,
        long_term_infiltration=1.5e-7,
        n_roots=40,
        n_depths=50,
    )
    rainfall = us.RainfallSeries.constant(4.5, 24.0)
    strength = us.StrengthProperties(
        cohesion=3.0e3,
        friction_angle_deg=33.6,
        unit_weight_soil=19.0e3,
    )
    times = np.linspace(0.0, 24.0 * 3600.0, 25)

    finite = us.run_single_column(
        hydro=hydro,
        rainfall=rainfall,
        strength=strength,
        boundary=us.FinitePerchedBoundary(),
        output_times=(times) / 3600.0,
    )
    infinite = us.run_single_column(
        hydro=hydro,
        rainfall=rainfall,
        strength=strength,
        boundary=us.InfiniteDrainageBoundary(),
        output_times=(times) / 3600.0,
    )

    assert finite.time_dataframe().shape[0] == len(times)
    assert infinite.time_dataframe().shape[0] == len(times)
    assert finite.arrays()["factor_of_safety"].shape == (len(times), hydro.n_depths + 1)
    assert infinite.arrays()["factor_of_safety"].shape == (len(times), hydro.n_depths + 1)


@pytest.mark.parametrize('water_table_depth', [1.3, 1.0, .75, None])
@pytest.mark.parametrize('equation_set',['baum_etal','srivastava_yeh'])
def test_failure_plane_initial_pressure_above_at_and_below_water_table(water_table_depth,equation_set):
    hydro = us.HydraulicProperties(
        slope_angle_deg=40., maximum_depth=1.5, saturated_hydraulic_conductivity=1e-5, alpha=3., theta_s=.41, theta_r=.06,
        initial_water_table_depth=water_table_depth, n_roots=40, n_depths=40,
        equation_set=equation_set,
    )
    result = us.run_single_column(
        hydro=hydro, rainfall=us.RainfallSeries.constant(0., 24.),
        strength=us.StrengthProperties(3000., 33.6, 18000.),
        boundary=us.FinitePerchedBoundary(), depth_to_failure_plane=1.0, output_times=(np.array([0., 3600.])) / 3600.0,
    )
    depth = hydro.maximum_depth if water_table_depth is None else water_table_depth
    # The consistent closure preserves the cos² gradient above and below
    # the water table. Legacy pressure above it has a cos gradient.
    gap = 1.0 - depth
    gradient = (np.cos(np.deg2rad(40.)) if gap < 0 and equation_set=='baum_etal'
                else np.cos(np.deg2rad(40.))**2)
    np.testing.assert_allclose(result.pressure_head_at_failure_plane, gradient * gap, atol=1e-12)
    np.testing.assert_array_equal(result.water_table_depth, np.full(2, depth))
    np.testing.assert_array_equal(result.water_table_rise, np.zeros(2))
    assert result.profiles[0]['z'][-1] == hydro.maximum_depth


@pytest.mark.parametrize("equation_set", ["baum_etal", "srivastava_yeh"])
def test_storage_refinement_does_not_require_extra_output_profiles(equation_set):
    hydro = us.HydraulicProperties(
        slope_angle_deg=40.0, maximum_depth=1.0, saturated_hydraulic_conductivity=1e-5, alpha=3.0, theta_s=0.41, theta_r=0.06,
        initial_water_table_depth=1.3, n_roots=40, n_depths=50,
        equation_set=equation_set,
    )
    rainfall = us.RainfallSeries(
        breakpoints=np.array([0.0, 2.0, 5.0, 8.0, 12.0]) * 3600,
        fluxes=np.array([1.0, 12.0, 1.0, 6.0]) / 3.6e6,
    )
    kwargs = dict(
        hydro=hydro, rainfall=rainfall,
        strength=us.StrengthProperties(3000.0, 33.6, 18000.0),
        boundary=us.FinitePerchedBoundary(drainage_efficiency=0.015),
        storage_time_step=(600.0) / 3600.0,
    )
    sparse = us.run_single_column(**kwargs, output_times=(np.array([0.0, 12 * 3600.0])) / 3600.0)
    hourly = us.run_single_column(**kwargs, output_times=(np.arange(13) * 3600.0) / 3600.0)
    assert len(sparse.profiles) == 2
    assert len(hourly.profiles) == 13
    np.testing.assert_array_equal((sparse.time) * 3600.0, np.array([0.0, 12 * 3600.0]))
    assert sparse.cumulative_excess[-1] > 0.0
    for field in ("cumulative_excess", "water_table_rise",
                  "pressure_head_at_failure_plane", "factor_of_safety_at_failure_plane"):
        np.testing.assert_allclose(getattr(sparse, field), getattr(hourly, field)[[0, -1]],
                                   rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("equation_set", ["baum_etal", "srivastava_yeh"])
def test_finite_base_clamps_deep_requested_nf_and_accumulates_from_base(equation_set):
    from dataclasses import replace
    h = us.HydraulicProperties(slope_angle_deg=40, maximum_depth=1,
        initial_water_table_depth=2, saturated_hydraulic_conductivity=1e-5,
        alpha=3, theta_s=.41, theta_r=.06, equation_set=equation_set,
        n_roots=40, n_depths=40)
    args = dict(rainfall=us.RainfallSeries.constant(4.5, 24),
        strength=us.StrengthProperties(3000, 33.6, 18000),
        boundary=us.FinitePerchedBoundary(drainage_efficiency=0),
        output_times=(np.array([0, 6*3600, 24*3600])) / 3600.0, storage_time_step=(600) / 3600.0)
    actual = us.run_single_column(hydro=h, **args)
    at_base = us.run_single_column(hydro=replace(h, initial_water_table_depth=1), **args)
    for key in actual.arrays():
        np.testing.assert_allclose(actual.arrays()[key], at_base.arrays()[key])
    assert actual.water_table_depth[0] == 1
    assert actual.pressure_head_at_failure_plane[0] == 0
    assert actual.cumulative_excess[0] == 0
    assert actual.water_table_rise[-1] > 0
    assert actual.metadata["requested_initial_water_table_depth_m"] == 2
    assert actual.metadata["initial_water_table_depth_was_clamped"]
    assert h.initial_water_table_depth == 2


def test_observation_plane_does_not_change_full_column_or_storage():
    h = us.HydraulicProperties(slope_angle_deg=40, maximum_depth=2,
        initial_water_table_depth=.8, saturated_hydraulic_conductivity=1e-5,
        alpha=3, theta_s=.41, theta_r=.06, n_roots=40, n_depths=40)
    args = dict(hydro=h, rainfall=us.RainfallSeries.constant(4.5, 24),
        strength=us.StrengthProperties(3000, 33.6, 18000),
        boundary=us.FinitePerchedBoundary(), output_times=(np.array([0, 12*3600])) / 3600.0)
    base = us.run_single_column(**args)
    selected = us.run_single_column(**args, depth_to_failure_plane=.5)
    for key in ("basal_flux", "cumulative_excess", "water_table_rise",
                "pressure_head_total", "theta", "hydraulic_conductivity", "factor_of_safety"):
        np.testing.assert_allclose(base.arrays()[key], selected.arrays()[key])
    assert selected.arrays()["z"][-1] == 2
    assert not np.allclose(base.pressure_head_at_failure_plane, selected.pressure_head_at_failure_plane)
    for plane in (0, 2.1, np.nan):
        with pytest.raises(ValueError, match="maximum_depth"):
            us.run_single_column(**args, depth_to_failure_plane=plane)


def test_infinite_base_retains_deeper_nf_but_limits_reported_profiles():
    h = us.HydraulicProperties(slope_angle_deg=40, maximum_depth=1,
        initial_water_table_depth=2, saturated_hydraulic_conductivity=1e-5,
        alpha=3, theta_s=.41, theta_r=.06, n_roots=40, n_depths=40)
    r = us.run_single_column(hydro=h, rainfall=us.RainfallSeries.constant(0, 24),
        strength=us.StrengthProperties(3000, 33.6, 18000),
        boundary=us.InfiniteDrainageBoundary(), output_times=(np.array([0, 3600])) / 3600.0)
    assert r.water_table_depth[0] == 2
    assert r.arrays()["z"][-1] == 1
    assert r.pressure_head_at_failure_plane[0] < 0
    assert not r.metadata["initial_water_table_depth_was_clamped"]


def test_hour_api_preserves_existing_physical_storm_results():
    from unsatstab.analytical import basal_flux_at_time
    h = us.HydraulicProperties(slope_angle_deg=40, maximum_depth=1,
        initial_water_table_depth=1, saturated_hydraulic_conductivity=1e-5,
        alpha=3, theta_s=.41, theta_r=.06, long_term_infiltration=1.5e-7,
        equation_set="baum_etal", n_roots=120, n_depths=120)
    rain = us.RainfallSeries.constant(4.5, 24)
    times_h = np.arange(97) / 4
    r = us.run_single_column(hydro=h, rainfall=rain,
        strength=us.StrengthProperties(3000, 33.6, 18000, max_factor_of_safety=2),
        boundary=us.FinitePerchedBoundary(), output_times=times_h,
        storage_time_step=5/60)
    np.testing.assert_array_equal(r.time, times_h)
    np.testing.assert_array_equal(r.time_dataframe().time_h, times_h)
    assert r.water_table_rise.max() == pytest.approx(.505691719, abs=1e-8)
    assert r.factor_of_safety_at_failure_plane.min() == pytest.approx(.9176313286, abs=1e-8)
    crossing = np.flatnonzero(r.factor_of_safety_at_failure_plane <= 1)[0]
    assert r.time[crossing] == 13.75
    assert basal_flux_at_time(24, h, rain) == pytest.approx(r.basal_flux[-1])
