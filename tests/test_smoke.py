import numpy as np
import unsatstab as us


def test_single_cell_smoke():
    hydro = us.HydraulicProperties(slope_angle_deg=40.0, maximum_depth=1.01, saturated_hydraulic_conductivity=1.0e-5, alpha=3.0, theta_s=0.41, theta_r=0.06, long_term_infiltration=1.5e-7, n_roots=40, n_depths=50)
    rainfall = us.RainfallSeries.constant(4.5, 24.0)
    strength = us.StrengthProperties(3.0e3, 33.6, 19.0e3)
    result = us.run_single_column(
        hydro=hydro,
        rainfall=rainfall,
        strength=strength,
        boundary=us.FinitePerchedBoundary(),
        output_times=(np.linspace(0, 24 * 3600, 25)) / 3600.0,
    )
    assert len(result.time_dataframe()) == 25
    assert result.arrays()["factor_of_safety"].shape == (25, 51)
