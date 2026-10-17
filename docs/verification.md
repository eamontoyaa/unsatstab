# Numerical verification

Verification checks the implementation against its equations, published results, and independent calculations. Validation evaluates model predictions against observations from a real application.

## Reference comparison notebooks

1. [Baum and Srivastava–Yeh equation sets](verification/01_baum_etal_vs_srivastava_yeh.ipynb): compares pressure and moisture conversions, water-table rise, stability, and early-time response under identical inputs.
2. [Gardner retention-curve comparison](verification/02_baum_etal_retention_curve.ipynb): compares the constitutive relation with approximate points from Baum et al. (2010), Figure 4, alongside Figure 9 of Baum et al. (2008).
3. [Space-time profile verification](verification/03_baum_etal_space_time.ipynb): compares Qva and Qtb pressure-head and factor-of-safety profiles with Baum et al. (2010), Figure 10, including numerical checkpoints.
4. [Srivastava–Yeh homogeneous soils](verification/04_srivastava_yeh_homogeneous_soils.ipynb): reproduces wetting and drainage for two soils (Figures 1–4) and their basal outflow (Figure 5), with an independent finite-difference comparison.

Each notebook contains its inputs, computed figures, reference comparisons, and interpretation. Readings from published figures are approximate; the Srivastava–Yeh numerical errors are measured against the independent solver. The Qva/Qtb cases use `baum_etal`; the homogeneous-soil cases use `srivastava_yeh`. See the [equation comparison](fundamentals.md#pressure-and-water-content-equation-sets).

## Automated checks

Tests cover analytical profiles and basal flux, rainfall conversions and infiltration limits, antecedent histories, storage integration, column geometry, material mapping, and agreement between spatial and individual-column runs.

Basal-flux tests also use reference outputs from the TRIGRS v2.1.0 Fortran routines under controlled forcing and early-time settings. Reproduction scripts and source checksums are in `research/trigrs-reference/`; the independent Srivastava–Yeh solver is defined directly in the [homogeneous-soil notebook](verification/04_srivastava_yeh_homogeneous_soils.ipynb#independent-numerical-comparison).

### Antecedent water-balance checks

Diagnostics compare effective rainfall with allocated drainage and reconstructed profile-water change. The residual reports numerical and model approximations; exact conservation is not imposed. Basal storage is represented in the reconstructed profile and must not be added as a separate water store. Storage above the fillable pore volume is flagged separately from rainfall rejected by the infiltration cap. See [diagnostic fields](inputs-outputs.md#antecedent-diagnostic-fields) and [water-table assumptions](fundamentals.md#water-table-rise).

## Scope

These checks support the selected equations and reference cases, rather than field predictive skill or equivalence with the complete TRIGRS model. The [Cerro Bravo application](validation/01_ruiz_etal_case_study.ipynb) demonstrates deterministic scenarios using real slope-unit geometry; it does not compare predictions with an observed landslide inventory. Full-model validation requires such observational evidence.
