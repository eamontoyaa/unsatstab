# unsatstab

![unsatstab logo: rainfall, an inclined soil column, a saturated zone, and basal drainage](assets/unsatstab_logo.svg){ .unsatstab-logo }

**Transient unsaturated infiltration and spatial slope-unit stability analysis.**

`unsatstab` is a Python package for physics-based modelling of rainfall-induced slope instability. It couples a one-dimensional analytical solution for transient unsaturated infiltration with a storage-based water-table model and an infinite-slope factor-of-safety formulation. The same modelling framework can be applied to individual slope columns or extended to spatial analyses over geomorphological slope units.

<div class="unsatstab-intro-end"></div>

## Installation

### Requirements

Python 3.12 or newer. NumPy, SciPy, pandas, Matplotlib, GeoPandas, and Joblib are among the dependencies installed with the package.

### Install from GitHub

```bash
python -m pip install 'git+https://github.com/eamontoyaa/unsatstab.git'
```

## Explore the documentation

- [Theory and foundations](fundamentals.md): model equations, assumptions, equation provenance, and references.
- [Inputs and outputs](inputs-outputs.md): column parameters, rainfall histories, slope-unit domains, results, and numerical controls.
- [Numerical verification](verification.md): reference comparisons and the scope of the evidence.
- [Cerro Bravo study-site application](validation/01_ruiz_etal_case_study.ipynb): deterministic scenarios for the Cerro Bravo case study and the scope of the validation evidence.
- [API reference](api.md): Python functions and classes, organized by module.

## Tutorials

- [Single-column response and boundaries](tutorials/01_single_column.ipynb)
- [Antecedent rainfall history](tutorials/02_antecedent_wetness.ipynb)
- [Multiple columns over slope units](tutorials/03_multiple_columns_slope_units.ipynb)

## Citing unsatstab

If you use `unsatstab` in a scientific publication or in a technical report, cite the associated methodological work:

> Montoya-Araque, E. A., Ruiz, D. F., Castillo-Cardona, E., & Gamboa, M. (202X). *unsatstab: a physics-based model for transient unsaturated infiltration and rainfall-induced slope stability at the slope-unit scale* [journal paper under preparation].

> Montoya-Araque, E. A., Ruiz, D. F., Castillo-Cardona, E., & Gamboa, M. (2026). *Fragility Curves for the Probabilistic Analysis of Rainfall-Induced Landslides through Physics-Based Numerical Modeling*. Manuscript prepared for the XIX Colombian Geotechnical Congress.

For reproducibility, also cite the exact `unsatstab` version or repository commit, along with parameter tables, rainfall-to-infiltration assumptions, and numerical controls.

## Authors

- [Exneyder A. Montoya-Araque](https://www.researchgate.net/profile/Exneyder-A-Montoya-Araque)
- [Daniel F. Ruiz](https://www.researchgate.net/profile/Daniel-F-Ruiz-2)
- [Emanuel Castillo-Cardona](https://www.researchgate.net/profile/Emanuel-Castillo-Cardona)

Copyright © 2026 Exneyder A. Montoya-Araque, Daniel F. Ruiz, Emanuel Castillo-Cardona, Universidad EAFIT.
Distributed under the MIT License.
