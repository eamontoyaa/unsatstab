from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest
from scipy.integrate import simpson

import unsatstab as us
from unsatstab.analytical import (
    pressure_from_conductivity, water_content_from_pressure, initial_profile,
    profile_at_time, basal_flux_at_time, _conductivity_step_response,
    _half_space_step,
)
from unsatstab.roots import compute_roots

REFERENCE=json.loads((Path(__file__).parent/'fixtures/srivastava_yeh_1991_example1.json').read_text())


@pytest.mark.parametrize('slope',[0.,30.,40.,60.])
def test_intrinsic_gardner_round_trip_and_legacy_pressure(slope):
    h=us.HydraulicProperties(slope_angle_deg=slope,maximum_depth=1.,saturated_hydraulic_conductivity=1e-5,alpha=3.,theta_s=.41,theta_r=.06,equation_set='srivastava_yeh')
    ratios=np.array([1e-3,.1,.5,1.])
    psi=pressure_from_conductivity(ratios*h.saturated_hydraulic_conductivity,hydro=h)
    np.testing.assert_allclose(np.exp(h.alpha*psi),ratios)
    theta=water_content_from_pressure(psi,hydro=h)
    np.testing.assert_allclose((theta-h.theta_r)/(h.theta_s-h.theta_r),ratios)
    legacy=replace(h,equation_set='baum_etal')
    legacy_psi=pressure_from_conductivity(ratios*h.saturated_hydraulic_conductivity,hydro=legacy)
    np.testing.assert_allclose(legacy_psi,psi/np.cos(h.slope_angle_rad))


@pytest.mark.parametrize('background',[0.,1e-7,1e-6])
def test_steady_pressure_gradient_preserves_background_darcy_flux(background):
    h=us.HydraulicProperties(slope_angle_deg=40.,maximum_depth=1.,saturated_hydraulic_conductivity=1e-5,alpha=3.,theta_s=.41,theta_r=.06,long_term_infiltration=background,
                            equation_set='srivastava_yeh',n_depths=2000)
    p=initial_profile(h)
    gradient=np.gradient(p['pressure_head'],p['z'],edge_order=2)
    q=p['hydraulic_conductivity']*(1-gradient/np.cos(h.slope_angle_rad)**2)
    np.testing.assert_allclose(q,background,atol=2e-11)
    saturated_q=h.saturated_hydraulic_conductivity*(1-h.beta/np.cos(h.slope_angle_rad)**2)
    np.testing.assert_allclose(saturated_q,background,atol=1e-20)
    # The capacity is an external boundary assumption, unchanged by the selector.
    boundary=us.FinitePerchedBoundary()
    assert boundary.drainage_capacity(h)==boundary.drainage_capacity(replace(h,equation_set='baum_etal'))


@pytest.mark.parametrize('case',REFERENCE['cases'],ids=lambda c:f"{c['name']}-alpha={c['alpha_m_inv']}")
def test_example1_against_independent_finite_difference_solution(case):
    ks=.01/3600.
    h=us.HydraulicProperties(slope_angle_deg=0.,maximum_depth=1.,saturated_hydraulic_conductivity=ks,alpha=case['alpha_m_inv'],theta_s=case['theta_s'],theta_r=case['theta_r'],
                            long_term_infiltration=case['initial_rate_cm_h']*ks,n_roots=500,n_depths=200,
                            equation_set='srivastava_yeh')
    rain=us.RainfallSeries((np.array([0.,101*3600])) / 3600.0,np.array([case['applied_rate_cm_h']*ks]))
    roots=compute_roots(h.transformed_alpha,1.,h.n_roots)
    q=[]
    for i,t in enumerate(REFERENCE['times_hours']):
        p=profile_at_time((t*3600) / 3600.0,h,rain,roots=roots)
        for field,key,tolerance in [('conductivity','hydraulic_conductivity',2e-5*ks),
                                    ('pressure_head','pressure_head',2e-5),
                                    ('volumetric_water_content','volumetric_water_content',1e-5)]:
            values=np.interp(REFERENCE['sample_depths_m'],p['z'],p[key])
            np.testing.assert_allclose(values,case[field][i],rtol=1e-4,atol=tolerance)
        q.append(basal_flux_at_time((t*3600) / 3600.0,h,rain,roots=roots))
    np.testing.assert_allclose(q,case['basal_flux_m_s'],rtol=1e-3,atol=1e-4*ks)


@pytest.mark.parametrize('tstar',[1e-7,1e-5,1e-3])
def test_early_profile_satisfies_imposed_surface_flux(tstar):
    h=us.HydraulicProperties(slope_angle_deg=0.,maximum_depth=1.,saturated_hydraulic_conductivity=1.,alpha=10.,theta_s=.4,theta_r=.06,equation_set='srivastava_yeh')
    roots=compute_roots(10.,1.,250)
    dz=np.sqrt(tstar)/h.alpha*1e-3
    values=[_conductivity_step_response(tstar/h.time_factor,z,step_flux=.25,
             hydro=h,depth=1.,roots=roots) for z in [0.,dz,2*dz]]
    gradient=(-3*values[0]+4*values[1]-values[2])/(2*dz)
    np.testing.assert_allclose(values[0]-gradient/h.alpha,.25,rtol=2e-6,atol=1e-8)


def test_fixed_domain_water_balance_for_inclined_consistent_flow():
    h=us.HydraulicProperties(slope_angle_deg=40.,maximum_depth=1.,saturated_hydraulic_conductivity=.01/3600,alpha=1.,theta_s=.45,theta_r=.2,long_term_infiltration=.1*.01/3600,
                            equation_set='srivastava_yeh',n_depths=1000)
    r=us.RainfallSeries((np.array([0.,101*3600])) / 3600.0,np.array([.9*.01/3600]))
    roots=compute_roots(h.transformed_alpha,1.,250)
    times=np.linspace(0.,10*3600,401)
    outflow=simpson([basal_flux_at_time((t) / 3600.0,h,r,roots=roots) for t in times],x=times)
    start=profile_at_time((0.) / 3600.0,h,r,roots=roots);end=profile_at_time((times[-1]) / 3600.0,h,r,roots=roots)
    storage=np.trapezoid(end['volumetric_water_content']-start['volumetric_water_content'],end['z'])
    np.testing.assert_allclose(r.fluxes[0]*times[-1]-outflow,storage,atol=2e-6)


@pytest.mark.parametrize('tstar',[1e-7,1e-5])
def test_short_time_limit_satisfies_diffusion_advection_equation(tstar):
    x=2*np.sqrt(tstar)
    dx=np.sqrt(tstar)*1e-3
    dt=tstar*1e-3
    center=_half_space_step(tstar,x)
    lower=_half_space_step(tstar,x-dx);upper=_half_space_step(tstar,x+dx)
    temporal=(_half_space_step(tstar+dt,x)-_half_space_step(tstar-dt,x))/(2*dt)
    spatial=(upper-2*center+lower)/dx**2-(upper-lower)/(2*dx)
    np.testing.assert_allclose(temporal,spatial,rtol=2e-5)


def test_invalid_equation_set_is_rejected():
    with pytest.raises(ValueError,match='equation_set'):
        us.HydraulicProperties(slope_angle_deg=40.,maximum_depth=1.,saturated_hydraulic_conductivity=1e-5,alpha=3.,theta_s=.41,theta_r=.06,equation_set='unknown')


def test_consistent_retention_reduces_antecedent_fixed_domain_balance_error():
    residuals={}
    for equation_set in ['baum_etal','srivastava_yeh']:
        h=us.HydraulicProperties(slope_angle_deg=40.,maximum_depth=1.3,saturated_hydraulic_conductivity=1e-5,alpha=3.,theta_s=.41,theta_r=.06,initial_water_table_depth=1.3,
                                n_roots=80,n_depths=80,equation_set=equation_set)
        result=us.run_single_column(hydro=h,rainfall=us.RainfallSeries.constant(4.5,24.),
            strength=us.StrengthProperties(3000.,33.6,18000.),boundary=us.FinitePerchedBoundary(),
            output_times=(np.array([0.])) / 3600.0,antecedent=us.AntecedentRainfall.uniform(500.,(90.) * 24.0,drainage_efficiency=1.))
        assert result.water_table_rise[0]==0.
        residuals[equation_set]=abs(result.metadata['antecedent_water_balance_residual_mm'])
    assert residuals['srivastava_yeh'] < .001  # numerical resolution, in mm
    assert residuals['baum_etal'] > .01
