from dataclasses import replace

import numpy as np
import pytest

import unsatstab as us
from unsatstab.analytical import initial_profile, profile_at_time, basal_flux_at_time
from unsatstab.solver import _run_history


@pytest.fixture
def hydro():
    return us.HydraulicProperties(
        slope_angle_deg=40, maximum_depth=1.3, saturated_hydraulic_conductivity=1e-5, alpha=3, theta_s=.41, theta_r=.06, long_term_infiltration=1.5e-7,
        initial_water_table_depth=1.3, n_roots=80, n_depths=80,
    )


def test_antecedent_volume_duration_and_event_transition():
    event = us.RainfallSeries.from_intensities(np.array([4.5, 2.]), 1, time_units="h")
    history = event.with_antecedent_rainfall(500, (90) * 24.0)
    end = 90 * 86400.
    assert history.fluxes[0] * 3.6e6 == pytest.approx(500 / (90 * 24))
    assert history.fluxes[0] * end == pytest.approx(.5)
    assert history.flux_at((end - 1) / 3600.0) == history.fluxes[0]
    assert history.flux_at((end) / 3600.0) == event.flux_at((0) / 3600.0)
    np.testing.assert_array_equal(history.breakpoints[1:] * 3600.0 - end, event.breakpoints * 3600.0)
    assert history.flux_at(history.breakpoints[-1]) == 0


def test_antecedent_preserves_delayed_event_gap():
    event = us.RainfallSeries((np.array([3600., 7200.])) / 3600.0, np.array([1e-6]))
    history = event.with_antecedent_rainfall(500, (90) * 24.0)
    assert history.flux_at((90 * 86400.) / 3600.0) == 0
    assert history.flux_at((90 * 86400. + 3600) / 3600.0) == 1e-6


@pytest.mark.parametrize("total,days", [(-1,90),(np.nan,90),(500,0),(500,-1),(500,np.inf)])
def test_invalid_antecedent(total, days):
    with pytest.raises(ValueError):
        us.RainfallSeries.constant(4.5,24).with_antecedent_rainfall(total,(days) * 24.0)


def test_interval_splitting_does_not_change_analytical_state(hydro):
    one = us.RainfallSeries.constant(4.5,24)
    split = us.RainfallSeries.from_intensities(np.array([4.5,4.5]),12,time_units="h")
    for t in [0, 1, 12*3600, 12*3600+1, 24*3600, 25*3600, 90*86400]:
        assert basal_flux_at_time((t) / 3600.0,hydro,one) == pytest.approx(basal_flux_at_time((t) / 3600.0,hydro,split))
        for key, value in profile_at_time((t) / 3600.0,hydro,one).items():
            np.testing.assert_allclose(value,profile_at_time((t) / 3600.0,hydro,split)[key],atol=1e-12)


def test_matching_equilibrium_rate_does_not_wet_profile(hydro):
    rain = us.RainfallSeries.constant(hydro.long_term_infiltration,90*86400,
                                    intensity_units="m/s",duration_units="s")
    for t in [1., 30*86400., 90*86400.]:
        assert basal_flux_at_time((t) / 3600.0,hydro,rain) == pytest.approx(hydro.long_term_infiltration)
        for key,value in initial_profile(hydro).items():
            np.testing.assert_allclose(profile_at_time((t) / 3600.0,hydro,rain)[key],value,atol=1e-12)


def test_long_antecedent_reaches_correct_steady_state_and_transition_is_continuous(hydro):
    event = us.RainfallSeries.constant(4.5,24)
    history = event.with_antecedent_rainfall(500,(90) * 24.0)
    transition = 90*86400.
    expected = initial_profile(replace(hydro,long_term_infiltration=history.fluxes[0]))
    state = profile_at_time((transition) / 3600.0,hydro,history)
    for key,value in expected.items():
        np.testing.assert_allclose(state[key],value,atol=1e-12)
        # A flux step changes surface head as sqrt(elapsed time), not linearly.
        np.testing.assert_allclose(profile_at_time((transition+1e-8) / 3600.0,hydro,history)[key],state[key],atol=1e-6)
    assert basal_flux_at_time((transition) / 3600.0,hydro,history) == pytest.approx(history.fluxes[0])
    assert basal_flux_at_time((transition+1e-6) / 3600.0,hydro,history) == pytest.approx(history.fluxes[0])


def test_drying_after_event_recovers_zero_flux_steady_profile(hydro):
    rain = us.RainfallSeries.constant(4.5,24)
    time = 90*86400.
    assert basal_flux_at_time((time) / 3600.0,hydro,rain) == pytest.approx(0,abs=1e-18)
    for key,value in initial_profile(replace(hydro,long_term_infiltration=0.)).items():
        np.testing.assert_allclose(profile_at_time((time) / 3600.0,hydro,rain)[key],value,atol=1e-12)


def test_event_result_matches_continuous_history(hydro):
    event = us.RainfallSeries.constant(4.5,24)
    times = np.arange(25)*3600.
    kwargs = dict(hydro=hydro,strength=us.StrengthProperties(3000,33.6,18000),
                  boundary=us.FinitePerchedBoundary(),storage_time_step=(3600.) / 3600.0)
    result = us.run_single_column(**kwargs, rainfall=event, output_times=(times) / 3600.0, antecedent=us.AntecedentRainfall.uniform(500, (90) * 24.0))
    continuous = us.run_single_column(
        **kwargs,rainfall=event.with_antecedent_rainfall(500,(90) * 24.0),output_times=(np.r_[0,90*86400.+times]) / 3600.0,
    )
    for key,value in result.arrays().items():
        if key == "time":
            np.testing.assert_array_equal(value,times / 3600.0)
        elif key == "z":
            np.testing.assert_array_equal(value,continuous.arrays()[key])
        else:
            np.testing.assert_allclose(value,continuous.arrays()[key][1:])
    assert len(result.profiles) == len(times)
    assert result.metadata["event_initial_water_table_depth_m"] == result.water_table_depth[0]


def test_dry_history_creates_no_negative_storage_debt(hydro):
    event = us.RainfallSeries.constant(4.5,24)
    result = us.run_single_column(
        hydro=hydro,
        rainfall=event,
        strength=us.StrengthProperties(3000, 33.6, 18000),
        boundary=us.FinitePerchedBoundary(),
        output_times=(np.array([0.0, 24 * 3600.0])) / 3600.0,
        antecedent=us.AntecedentRainfall.uniform(0, (90) * 24.0),
    )
    assert result.cumulative_excess[0] == 0
    assert result.cumulative_excess[-1] > 0


def test_antecedent_changes_failure_time_and_quadrature_converges(hydro):
    hydro = replace(hydro,long_term_infiltration=0.)
    kwargs = dict(hydro=hydro,rainfall=us.RainfallSeries.constant(4.5,24),
                  strength=us.StrengthProperties(3000,33.6,18000),
                  boundary=us.FinitePerchedBoundary(drainage_efficiency=.02),
                  output_times=(np.arange(97)*900.) / 3600.0)
    dry = us.run_single_column(**kwargs, antecedent=us.AntecedentRainfall.uniform(0, (90) * 24.0))
    wet = us.run_single_column(**kwargs, antecedent=us.AntecedentRainfall.uniform(500, (90) * 24.0))
    refined = us.run_single_column(**kwargs, antecedent=us.AntecedentRainfall.uniform(500, (90) * 24.0), storage_time_step=(1800) / 3600.0)
    assert wet.profiles[0]["volumetric_water_content"][0] > dry.profiles[0]["volumetric_water_content"][0]
    assert wet.factor_of_safety_at_failure_plane[-1] < dry.factor_of_safety_at_failure_plane[-1]
    wet_failure = (wet.time * 3600.0)[wet.factor_of_safety_at_failure_plane <= 1]
    dry_failure = (dry.time * 3600.0)[dry.factor_of_safety_at_failure_plane <= 1]
    assert len(wet_failure) > 0
    assert len(dry_failure) == 0 or wet_failure[0] < dry_failure[0]
    np.testing.assert_allclose(wet.cumulative_excess,refined.cumulative_excess,atol=1e-5)
    np.testing.assert_allclose(wet.factor_of_safety_at_failure_plane,refined.factor_of_safety_at_failure_plane,atol=1e-4)


@pytest.mark.parametrize("total", [0,500,1000])
def test_efficient_antecedent_preserves_nf_for_moderate_rainfall(hydro,total):
    hydro = replace(hydro,long_term_infiltration=0.)
    event = us.RainfallSeries.constant(4.5,24)
    result = us.run_single_column(
        hydro=hydro,
        rainfall=event,
        strength=us.StrengthProperties(3000, 33.6, 18000),
        boundary=us.FinitePerchedBoundary(drainage_efficiency=0.015),
        output_times=(np.array([0.0, 3600.0])) / 3600.0,
        antecedent=us.AntecedentRainfall.uniform(total, (90) * 24.0, drainage_efficiency=1.0),
    )
    expected = profile_at_time((90*86400.) / 3600.0,hydro,event.with_antecedent_rainfall(total,(90) * 24.0))
    state = result.profiles[0]
    for key in ["hydraulic_conductivity","volumetric_water_content"]:
        np.testing.assert_allclose(state[key],expected[key])
    np.testing.assert_allclose(state["pressure_head_total"],expected["pressure_head"])
    assert result.water_table_depth[0] == hydro.unsaturated_domain_depth
    assert result.water_table_rise[0] == 0
    assert result.cumulative_excess[0] == 0
    assert result.metadata["storage_reference"] == "antecedent_start"
    assert result.metadata["antecedent_effective_infiltration_mm"] == pytest.approx(total)
    expected_water = np.trapezoid(
        expected["volumetric_water_content"]-initial_profile(hydro)["volumetric_water_content"],
        state["z"],
    )*1000
    assert result.metadata["antecedent_profile_water_change_mm"] == pytest.approx(expected_water)


def test_two_stage_drainage_carries_transient_memory_into_event(hydro):
    hydro = replace(hydro,long_term_infiltration=0.)
    event = us.RainfallSeries.constant(4.5,24)
    times = np.array([0.,3600.,24*3600.])
    result = us.run_single_column(
        hydro=hydro,
        rainfall=event,
        strength=us.StrengthProperties(3000, 33.6, 18000),
        boundary=us.FinitePerchedBoundary(),
        output_times=(times) / 3600.0,
        antecedent=us.AntecedentRainfall.uniform(500, (0.5) * 24.0, drainage_efficiency=0.7),
    )
    full = _run_history(
        hydro=hydro,rainfall=event.with_antecedent_rainfall(500,(.5) * 24.0),
        strength=us.StrengthProperties(3000,33.6,18000),boundary=us.FinitePerchedBoundary(),
        output_times=(np.r_[0.,.5*86400+times]) / 3600.0,antecedent_drainage_efficiency=.7,
        antecedent_end_time_h=(.5*86400) / 3600.0,
        storage_time_step=(3600) / 3600.0,
    )
    for key in ["basal_flux","cumulative_excess","water_table_depth","factor_of_safety"]:
        np.testing.assert_allclose(result.arrays()[key],full.arrays()[key][1:])
    # A newly imposed steady state would erase the antecedent transient.
    steady = initial_profile(replace(hydro,long_term_infiltration=min(
        500/1000/(.5*86400),hydro.saturated_hydraulic_conductivity,
    )))
    transient = profile_at_time((.5*86400) / 3600.0, hydro, event.with_antecedent_rainfall(500,(.5) * 24.0).capped(
        hydro.saturated_hydraulic_conductivity,
    ))
    np.testing.assert_allclose(result.profiles[0]["pressure_head_unsaturated"], transient["pressure_head"])
    assert not np.allclose(transient["volumetric_water_content"], steady["volumetric_water_content"])


def test_antecedent_duration_changes_storage_and_surface_cap_is_reported(hydro):
    hydro = replace(hydro,long_term_infiltration=0.)
    kwargs = dict(hydro=hydro,rainfall=us.RainfallSeries.constant(4.5,24),
        strength=us.StrengthProperties(3000,33.6,18000),boundary=us.FinitePerchedBoundary(),
        output_times=(np.array([0.])) / 3600.0)
    long = us.run_single_column(**kwargs, antecedent=us.AntecedentRainfall.uniform(500, (90) * 24.0, drainage_efficiency=1.0))
    short = us.run_single_column(**kwargs, antecedent=us.AntecedentRainfall.uniform(500, (1) * 24.0, drainage_efficiency=1.0))
    capped = us.run_single_column(
        **kwargs,
        antecedent=us.AntecedentRainfall.uniform(500, (0.1) * 24.0, drainage_efficiency=1.0),
        storage_time_step=(60) / 3600.0,
    )
    assert short.metadata["antecedent_profile_water_change_mm"] > long.metadata["antecedent_profile_water_change_mm"]
    assert short.profiles[0]["hydraulic_conductivity"][0] > long.profiles[0]["hydraulic_conductivity"][0]
    assert capped.metadata["antecedent_effective_infiltration_mm"] == pytest.approx(86.4)
    assert capped.metadata["antecedent_rainfall_excess_mm"] == pytest.approx(413.6)
    assert capped.metadata["antecedent_infiltration_mm"] == pytest.approx(86.4)
    assert long.metadata["antecedent_drainage_mm"] > 490
    # Expose the closure residual instead of replacing analytical drainage by
    # rainfall minus retained profile water to conceal the approximation.
    assert long.metadata["antecedent_water_balance_residual_mm"] == pytest.approx(
        long.metadata["antecedent_effective_infiltration_mm"]
        - long.metadata["antecedent_drainage_mm"]
        - long.metadata["antecedent_profile_water_change_mm"]
    )
    assert capped.water_table_depth[0] <= hydro.unsaturated_domain_depth


def test_antecedent_efficiency_controls_storage_and_defaults_to_event_boundary(hydro):
    kwargs = dict(hydro=replace(hydro,long_term_infiltration=0.),
        rainfall=us.RainfallSeries.constant(4.5,24),strength=us.StrengthProperties(3000,33.6,18000),
        boundary=us.FinitePerchedBoundary(drainage_efficiency=.015),output_times=(np.array([0.,3600.])) / 3600.0)
    results = [us.run_single_column(**kwargs, antecedent=us.AntecedentRainfall.uniform(500, (90) * 24.0, drainage_efficiency=eta))
               for eta in [0.,.01,.015,.02,1.]]
    storage = np.array([r.cumulative_excess[0] for r in results])
    nf = np.array([r.water_table_depth[0] for r in results])
    assert np.all(np.diff(storage) <= 1e-12)
    assert np.all(np.diff(nf) >= -1e-12)
    assert storage[0] > storage[2] > storage[-1]
    assert storage[-1] == 0
    assert results[0].water_table_depth[0] == 0
    assert results[0].metadata["antecedent_storage_above_pore_capacity_mm"] > 0
    assert results[-1].metadata["antecedent_storage_above_pore_capacity_mm"] == 0
    assert results[0].metadata["antecedent_drainage_mm"] == pytest.approx(0,abs=1e-8)
    for result in results:
        assert result.metadata["antecedent_basal_storage_mm"] == result.cumulative_excess[0]*1000
        assert result.metadata["storage_reference"] == "antecedent_start"
        assert result.metadata["antecedent_basal_transfer_mm"] == pytest.approx(
            result.metadata["antecedent_drainage_mm"] + result.metadata["antecedent_basal_storage_mm"],
        )
        assert np.isfinite(result.metadata["antecedent_water_balance_residual_mm"])
        assert result.drainage_capacity[0] == kwargs["boundary"].drainage_capacity(kwargs["hydro"])
    default = us.run_single_column(**kwargs, antecedent=us.AntecedentRainfall.uniform(500, (90) * 24.0))
    np.testing.assert_allclose(default.cumulative_excess,results[2].cumulative_excess)
    assert default.metadata["antecedent_drainage_efficiency"] == .015


def test_unit_efficiency_does_not_enforce_fixed_nf_under_intense_rain(hydro):
    result = us.run_single_column(
        hydro=replace(hydro, long_term_infiltration=0.0),
        rainfall=us.RainfallSeries.constant(4.5, 24),
        strength=us.StrengthProperties(3000, 33.6, 18000),
        boundary=us.FinitePerchedBoundary(),
        output_times=(np.array([0.0])) / 3600.0,
        storage_time_step=(900.0) / 3600.0,
        antecedent=us.AntecedentRainfall.uniform(500, (1) * 24.0, drainage_efficiency=1.0),
    )
    assert result.cumulative_excess[0] > 0
    assert result.water_table_depth[0] < hydro.unsaturated_domain_depth


def test_drainage_switch_has_no_endpoint_leak_into_antecedent_storage(hydro):
    hydro = replace(hydro,long_term_infiltration=0.)
    result = us.run_single_column(
        hydro=hydro,
        rainfall=us.RainfallSeries.constant(4.5, 24),
        strength=us.StrengthProperties(3000, 33.6, 18000),
        boundary=us.FinitePerchedBoundary(drainage_efficiency=0.0),
        output_times=(np.array([0.0, 3600.0])) / 3600.0,
        storage_time_step=(90 * 86400.0) / 3600.0,
        antecedent=us.AntecedentRainfall.uniform(500, (90) * 24.0, drainage_efficiency=1.0),
    )
    # Even with coarse quadrature, the final antecedent endpoint uses antecedent
    # capacity. A right-sided switch there would inject spurious storage at t=0.
    assert result.cumulative_excess[0] == 0
    assert result.cumulative_excess[1] > 0
    assert result.drainage_capacity[0] == 0


@pytest.mark.parametrize("eta", [-.1,1.1,np.nan,np.inf])
def test_invalid_antecedent_efficiency(hydro,eta):
    with pytest.raises(ValueError,match="drainage_efficiency"):
        us.run_single_column(
            hydro=hydro,
            rainfall=us.RainfallSeries.constant(4.5, 24),
            strength=us.StrengthProperties(3000, 33.6, 18000),
            boundary=us.FinitePerchedBoundary(),
            output_times=(np.array([0.0])) / 3600.0,
            antecedent=us.AntecedentRainfall.uniform(500, (90) * 24.0, drainage_efficiency=eta),
        )


@pytest.mark.parametrize("rainfall", [None, 500, "rainfall"])
def test_invalid_antecedent_series_type(rainfall):
    with pytest.raises(TypeError, match="RainfallSeries"):
        us.AntecedentRainfall(rainfall)


def test_antecedent_series_must_start_at_zero():
    with pytest.raises(ValueError, match="begin at zero"):
        us.AntecedentRainfall(us.RainfallSeries((np.array([1.,2.])) / 3600.0,np.array([0.])))


@pytest.mark.parametrize("antecedent", [500, {}, us.RainfallSeries.constant(1,1)])
def test_solver_requires_antecedent_object(hydro,antecedent):
    with pytest.raises(TypeError, match="AntecedentRainfall"):
        us.run_single_column(
            hydro=hydro,rainfall=us.RainfallSeries.constant(4.5,24),
            strength=us.StrengthProperties(3000,33.6,18000),
            boundary=us.FinitePerchedBoundary(),output_times=(np.array([0.,3600.])) / 3600.0,
            antecedent=antecedent,
        )


@pytest.mark.parametrize("times", [[],[1],[0,0],[0,np.nan],[[0,1]]])
def test_invalid_event_times(hydro,times):
    with pytest.raises(ValueError):
        us.run_single_column(
            hydro=hydro,
            rainfall=us.RainfallSeries.constant(4.5, 24),
            strength=us.StrengthProperties(3000, 33.6, 18000),
            boundary=us.FinitePerchedBoundary(),
            output_times=(np.array(times)) / 3600.0,
            antecedent=us.AntecedentRainfall.uniform(500, (90) * 24.0),
        )


@pytest.mark.parametrize("step", [0,-1,np.nan,np.inf])
def test_invalid_storage_step(hydro,step):
    with pytest.raises(ValueError):
        us.run_single_column(
            hydro=hydro,rainfall=us.RainfallSeries.constant(4.5,24),
            strength=us.StrengthProperties(3000,33.6,18000),boundary=us.FinitePerchedBoundary(),
            output_times=(np.array([0.,3600.])) / 3600.0,storage_time_step=(step) / 3600.0,
        )


def test_optional_antecedent_none_preserves_event_results(hydro):
    kwargs = dict(hydro=hydro,rainfall=us.RainfallSeries.constant(4.5,24),
        strength=us.StrengthProperties(3000,33.6,18000),boundary=us.FinitePerchedBoundary(),
        output_times=(np.array([0.,3600.,86400.])) / 3600.0)
    ordinary = us.run_single_column(**kwargs)
    explicit = us.run_single_column(**kwargs,antecedent=None)
    reference = _run_history(**kwargs)
    for key, value in ordinary.arrays().items():
        np.testing.assert_array_equal(value,explicit.arrays()[key])
        np.testing.assert_array_equal(value,reference.arrays()[key])


def test_custom_antecedent_hyetograph_keeps_timing_and_local_cap(hydro):
    hydro = replace(hydro,long_term_infiltration=0.)
    # Different intervals, including dry days and a pulse above local Ks.
    previous = us.RainfallSeries.from_intensities(np.array([0.,50.,2.,0.]),6,time_units="h")
    antecedent = us.AntecedentRainfall(previous)
    event = us.RainfallSeries.constant(4.5,24)
    times = np.array([0.,3600.,86400.])
    kwargs = dict(hydro=hydro,strength=us.StrengthProperties(3000,33.6,18000),
        boundary=us.FinitePerchedBoundary(),storage_time_step=(1800.) / 3600.0)
    result = us.run_single_column(**kwargs,rainfall=event,output_times=(times) / 3600.0,antecedent=antecedent)
    history = antecedent.prepend_to(event)
    reference = us.run_single_column(**kwargs,rainfall=history,output_times=(np.r_[0.,86400.+times]) / 3600.0)
    for key, value in result.arrays().items():
        if key == "time": np.testing.assert_array_equal(value,times / 3600.0)
        elif key == "z": np.testing.assert_array_equal(value,reference.arrays()[key])
        else: np.testing.assert_allclose(value,reference.arrays()[key][1:])
    assert antecedent.total_rainfall_mm == pytest.approx(312.)
    assert antecedent.duration_days == 1.
    assert result.metadata["antecedent_effective_infiltration_mm"] == pytest.approx(228.)
    assert result.metadata["antecedent_rainfall_excess_mm"] == pytest.approx(84.)
    assert result.metadata["antecedent_infiltration_flux_m_s"] == pytest.approx(.228/86400)
    expected = np.trapezoid(
        result.profiles[0]["volumetric_water_content"]-reference.profiles[0]["volumetric_water_content"],
        result.profiles[0]["z"],
    ) * 1000
    assert result.metadata["antecedent_profile_water_change_mm"] == pytest.approx(expected)


@pytest.mark.parametrize("step", [0,-1,np.nan,np.inf])
def test_invalid_storage_step_with_antecedent(hydro,step):
    with pytest.raises(ValueError, match="storage_time_step"):
        us.run_single_column(
            hydro=hydro,rainfall=us.RainfallSeries.constant(4.5,24),
            strength=us.StrengthProperties(3000,33.6,18000),boundary=us.FinitePerchedBoundary(),
            output_times=(np.array([0.,3600.])) / 3600.0,storage_time_step=(step) / 3600.0,
            antecedent=us.AntecedentRainfall.uniform(500,(90) * 24.0),
        )
