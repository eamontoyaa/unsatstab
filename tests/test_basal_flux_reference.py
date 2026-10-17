"""Basal-flux checks against original, compiled USGS TRIGRS 2.1 routines."""
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest

import unsatstab as us
from unsatstab.analytical import basal_flux_at_time

REFERENCE = json.loads((Path(__file__).parent / 'fixtures' / 'trigrs_2_1_basal_flux.json').read_text())


def setup_case(case, n_roots=250):
    hydro = us.HydraulicProperties(slope_angle_deg=0., maximum_depth=1., saturated_hydraulic_conductivity=1., alpha=case['alpha_depth'], theta_s=1., theta_r=0., n_roots=n_roots, equation_set='baum_etal')
    rainfall = us.RainfallSeries(
        np.array(case['breakpoints_t_star']) / hydro.time_factor,
        np.array(case['rates']),
    )
    return hydro, rainfall


@pytest.mark.parametrize('case', REFERENCE['cases'],
                         ids=lambda c: f"{c['name']}-aD={c['alpha_depth']}")
def test_basal_flux_matches_compiled_trigrs_series(case):
    hydro, rainfall = setup_case(case)
    actual = [basal_flux_at_time(t / hydro.time_factor, hydro, rainfall)
              for t in case['times_t_star']]
    # The original flux routine terminates its series at 1e-6 and uses some
    # single-precision inputs. Compare absolute error relative to the forcing;
    # do not require relative accuracy when the physical response is near zero.
    np.testing.assert_allclose(actual, case['series_audit_basal_flux'],
                               rtol=1e-5, atol=8e-6 * max(case['rates']))
    # Comparison with the original, unmodified smt=0.1 configuration is only
    # appropriate after that early-time branch has ended for every active step.
    breaks = np.array(case['breakpoints_t_star'])
    mask = [all(t - b >= .1 for b in breaks if t > b)
            and t >= .1 for t in case['times_t_star']]
    np.testing.assert_allclose(np.array(actual)[mask],
                               np.array(case['original_basal_flux'])[mask],
                               rtol=1e-5, atol=8e-6 * max(case['rates']))


def test_shallow_column_does_not_delay_basal_response_until_time_threshold():
    case = next(c for c in REFERENCE['cases']
                if c['name'] == 'constant' and c['alpha_depth'] == .01)
    hydro, rainfall = setup_case(case)
    index = case['times_t_star'].index(9.9e-5)
    flux = basal_flux_at_time(9.9e-5 / hydro.time_factor, hydro, rainfall)
    # The old scalar cutoff returned zero here; the independently executed
    # eigenfunction solution already transfers most of the imposed 0.25 flux.
    assert case['original_basal_flux'][index] == 0.
    assert flux > .2
    np.testing.assert_allclose(flux, case['series_audit_basal_flux'][index], atol=2e-6)


@pytest.mark.parametrize('alpha_depth', [.01, .1, .5, 1., 2.])
def test_basal_branch_transition_is_continuous_and_root_convergent(alpha_depth):
    hydro = us.HydraulicProperties(slope_angle_deg=0., maximum_depth=1., saturated_hydraulic_conductivity=1., alpha=alpha_depth, theta_s=1., theta_r=0., equation_set='baum_etal')
    rainfall = us.RainfallSeries(np.array([0., 100.]) / hydro.time_factor, np.array([.25]))
    cutoff = hydro.small_time_threshold * min(1., alpha_depth**2)
    times_star = np.array([.999 * cutoff, 1.001 * cutoff, 2 * cutoff, .1, 1.])
    values = []
    for count in [250, 500, 1000]:
        h = replace(hydro, n_roots=count)
        values.append(np.array([basal_flux_at_time(t / h.time_factor, h, rainfall)
                                for t in times_star]))
    assert np.isfinite(values).all()
    assert np.abs(values[0][1] - values[0][0]) < 1e-9
    np.testing.assert_allclose(values[0], values[2], rtol=1e-8, atol=1e-9)
    np.testing.assert_allclose(values[1], values[2], rtol=1e-8, atol=1e-9)


def test_early_basal_flux_keeps_equilibrium_background():
    hydro = us.HydraulicProperties(slope_angle_deg=40., maximum_depth=1.01, saturated_hydraulic_conductivity=1e-5, alpha=3., theta_s=.41, theta_r=.06, long_term_infiltration=1.5e-7, equation_set='baum_etal')
    rainfall = us.RainfallSeries.constant(4.5, 24.)
    time = .5 * hydro.small_time_threshold / hydro.time_factor
    np.testing.assert_allclose(basal_flux_at_time(time, hydro, rainfall),
                               hydro.long_term_infiltration, rtol=0, atol=1e-20)
