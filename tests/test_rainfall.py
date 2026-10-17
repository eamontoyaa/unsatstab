import numpy as np
import pytest

import unsatstab as us


def test_rainfall_half_open_intervals_and_total():
    rain = us.RainfallSeries.from_intensities(np.array([4.5, 9.]), 12, time_units="h")
    assert rain.flux_at((-1) / 3600.0) == 0
    assert rain.flux_at((0) / 3600.0) == rain.fluxes[0]
    assert rain.flux_at((12 * 3600 - 1e-6) / 3600.0) == rain.fluxes[0]
    assert rain.flux_at((12 * 3600) / 3600.0) == rain.fluxes[1]
    assert rain.flux_at((24 * 3600) / 3600.0) == 0
    assert rain.flux_at((25 * 3600) / 3600.0) == 0
    constant = us.RainfallSeries.constant(4.5, 24)
    assert np.sum(constant.fluxes * np.diff(constant.breakpoints)) * 3600.0 == pytest.approx(0.108)


def test_intensity_duration_pairs_preserve_variable_event_and_dry_gaps():
    rates = np.array([0, 3.5, 6, 10, 13, 9, 5, 3.5, 1.5, 4, 6.5, 10, 7.5, 1])
    durations = np.array([3, 2, 1, 1, .5, 1, 2, 4, 2, 2.5, 1, 1, 2, 1])
    rain = us.RainfallSeries.from_intervals(rates, durations)
    np.testing.assert_allclose(rain.breakpoints, np.r_[0, np.cumsum(durations)])
    assert rain.duration == 24
    assert np.dot(np.diff(rain.breakpoints), rain.fluxes) * 3.6e6 == pytest.approx(108)
    assert rain.flux_at((3 * 3600 - 1) / 3600.0) == 0
    assert rain.flux_at((3 * 3600) / 3600.0) * 3.6e6 == pytest.approx(3.5)
    assert rain.flux_at((24 * 3600) / 3600.0) == 0
    minutes = us.RainfallSeries.from_intervals(rates / 3.6e6, durations * 60,
                                              intensity_units="m/s", duration_units="min")
    np.testing.assert_allclose(minutes.breakpoints, rain.breakpoints)
    np.testing.assert_allclose(minutes.fluxes, rain.fluxes)


@pytest.mark.parametrize("rates,durations", [
    ([], []), ([1, 2], [1]), ([[1]], [1]), ([1], [[1]]),
    ([-1], [1]), ([np.nan], [1]), ([1], [0]), ([1], [-1]), ([1], [np.inf]),
])
def test_intensity_duration_pairs_reject_invalid_intervals(rates, durations):
    with pytest.raises(ValueError):
        us.RainfallSeries.from_intervals(rates, durations)


@pytest.mark.parametrize("units", [{"intensity_units": "cm/h"}, {"duration_units": "days"}])
def test_intensity_duration_pairs_reject_unknown_units(units):
    with pytest.raises(ValueError):
        us.RainfallSeries.from_intervals([1], [1], **units)


def test_all_public_time_vectors_and_durations_default_to_hours():
    direct = us.RainfallSeries(np.array([0., 12., 24.]), np.array([4.5, 9.]) / 3.6e6)
    regular = us.RainfallSeries.from_intensities(np.array([4.5, 9.]), 12.)
    intervals = us.RainfallSeries.from_intervals([4.5, 9.], [12., 12.])
    for rain in (direct, regular, intervals):
        np.testing.assert_array_equal(rain.breakpoints, [0., 12., 24.])
        assert rain.duration == 24.
        assert rain.flux_at(12.) * 3.6e6 == pytest.approx(9.)
        assert rain.flux_at(24.) == 0.
    seconds = us.RainfallSeries.constant(4.5, 86400., duration_units="s")
    np.testing.assert_array_equal(seconds.breakpoints, [0., 24.])
    minutes = us.RainfallSeries.constant(4.5, 1440., duration_units="min")
    np.testing.assert_array_equal(minutes.breakpoints, [0., 24.])
    antecedent = us.AntecedentRainfall.uniform(108., 24.)
    days = us.AntecedentRainfall.uniform(108., 1., duration_units="days")
    assert antecedent.duration == days.duration == 24.
    assert antecedent.total_rainfall_mm == pytest.approx(108.)
    np.testing.assert_allclose(antecedent.rainfall.fluxes, days.rainfall.fluxes)
