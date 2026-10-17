import numpy as np

import unsatstab as us
from unsatstab.analytical import gardner_water_content
from unsatstab.roots import characteristic_equation, compute_roots


def _strength():
    return us.StrengthProperties(
        cohesion=3_000.0,
        friction_angle_deg=33.6,
        unit_weight_soil=18_000.0,
    )


def test_characteristic_roots_satisfy_fortran_equation():
    alpha_transformed = 3.0 * np.cos(np.deg2rad(40.0)) ** 2
    roots = compute_roots(alpha_transformed, 1.01, 20)
    residual = [characteristic_equation(root, alpha_transformed * 1.01) for root in roots]
    assert np.max(np.abs(residual)) < 1.0e-7
    assert np.all(np.diff(roots) > 0.0)


def test_gardner_retention_curve_uses_intrinsic_alpha_for_baum_figure_4():
    theta = gardner_water_content(
        np.array([0.0, -1.0 / 3.5, -1.0]), alpha=3.5, theta_s=0.41, theta_r=0.05
    )
    np.testing.assert_allclose(theta[0], 0.41)
    np.testing.assert_allclose(theta[1], 0.05 + 0.36 / np.e)
    assert 0.05 < theta[-1] < theta[1]


def test_drainage_capacity_matches_trigrs_default_flow_direction():
    hydro = us.HydraulicProperties(slope_angle_deg=40.0, maximum_depth=1.01, saturated_hydraulic_conductivity=1.0e-5, alpha=3.0, theta_s=0.41, theta_r=0.06, long_term_infiltration=1.5e-7)
    expected = 0.1 * hydro.saturated_hydraulic_conductivity * np.sin(hydro.slope_angle_rad) ** 2
    np.testing.assert_allclose(
        us.FinitePerchedBoundary().drainage_capacity(hydro), expected, rtol=1.0e-12
    )


def test_qva_infinite_boundary_reference_case_remains_stable():
    hydro = us.HydraulicProperties(
        slope_angle_deg=45.0,
        maximum_depth=0.83,
        saturated_hydraulic_conductivity=5.0e-5,
        alpha=3.5,
        theta_s=0.41,
        theta_r=0.06,
        long_term_infiltration=1.5e-7,
        n_roots=80,
        n_depths=80,
        equation_set="baum_etal",
    )
    result = us.run_single_column(
        hydro=hydro,
        rainfall=us.RainfallSeries.constant(4.5, 24.0),
        strength=_strength(),
        boundary=us.InfiniteDrainageBoundary(),
        output_times=(np.arange(25, dtype=float) * 3600.0) / 3600.0,
    )
    np.testing.assert_allclose(result.factor_of_safety_at_failure_plane[[0, -1]], 1.06, atol=0.02)
    assert result.water_table_rise.max() == 0.0


def test_qtb_finite_boundary_tracks_baum_2010_figure_10():
    hydro = us.HydraulicProperties(
        slope_angle_deg=40.0,
        maximum_depth=1.01,
        saturated_hydraulic_conductivity=1.0e-5,
        alpha=3.0,
        theta_s=0.41,
        theta_r=0.06,
        long_term_infiltration=1.5e-7,
        n_roots=120,
        n_depths=120,
        equation_set="baum_etal",
    )
    result = us.run_single_column(
        hydro=hydro,
        rainfall=us.RainfallSeries.constant(4.5, 24.0),
        strength=_strength(),
        boundary=us.FinitePerchedBoundary(),
        output_times=(np.arange(25, dtype=float) * 3600.0) / 3600.0,
    )

    fs = result.factor_of_safety_at_failure_plane
    np.testing.assert_allclose(fs[0], 1.13, atol=0.02)
    np.testing.assert_allclose(fs[16], 1.00, atol=0.04)
    np.testing.assert_allclose(fs[24], 0.94, atol=0.04)
    first_failure_hour = int(np.flatnonzero(fs <= 1.0)[0])
    assert 13 <= first_failure_hour <= 17
