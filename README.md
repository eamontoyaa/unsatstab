# unsatstab

<img src="https://raw.githubusercontent.com/eamontoyaa/unsatstab/main/docs/assets/unsatstab_logo.svg" alt="unsatstab logo: rainfall, an inclined soil column, a saturated zone, and basal drainage" align="right" width="110">

**Transient unsaturated infiltration and spatial slope-unit stability analysis.**

`unsatstab` is a Python package for physics-based modelling of rainfall-induced slope instability. It couples a one-dimensional analytical solution for transient unsaturated infiltration with a storage-based water-table model and an infinite-slope factor-of-safety formulation. The same modelling framework can be applied to individual slope columns or extended to spatial analyses over geomorphological slope units.

[Documentation](https://eamontoyaa.github.io/unsatstab/)

## Installation

Python 3.12 or newer is required. Install from the repository:

```bash
python -m pip install 'git+https://github.com/eamontoyaa/unsatstab.git'
```

## Citing unsatstab

If you use `unsatstab` in a scientific publication or in a technical report, cite the associated methodological work:

> Montoya-Araque, E. A., Ruiz, D. F., Castillo-Cardona, E., & Gamboa, M. (202X). *unsatstab: a physics-based model for transient unsaturated infiltration and rainfall-induced slope stability at the slope-unit scale* [journal paper under preparation].

> Montoya-Araque, E. A., Ruiz, D. F., Castillo-Cardona, E., & Gamboa, M. (2026). *Fragility Curves for the Probabilistic Analysis of Rainfall-Induced Landslides through Physics-Based Numerical Modeling*. Manuscript prepared for the XIX Colombian Geotechnical Congress.

For reproducibility, also cite the exact `unsatstab` version or repository commit, along with parameter tables, rainfall-to-infiltration assumptions, and numerical controls.

## License

The Python implementation is distributed under the [MIT License](https://github.com/eamontoyaa/unsatstab/blob/main/LICENSE). Its scientific foundations are cited in the [documentation](https://eamontoyaa.github.io/unsatstab/fundamentals/#references).

## Authors

- [Exneyder A. Montoya-Araque](https://www.researchgate.net/profile/Exneyder-A-Montoya-Araque)
- [Daniel F. Ruiz](https://www.researchgate.net/profile/Daniel-F-Ruiz-2)
- [Emanuel Castillo-Cardona](https://www.researchgate.net/profile/Emanuel-Castillo-Cardona)

Copyright © 2026 Exneyder A. Montoya-Araque, Daniel F. Ruiz, Emanuel Castillo-Cardona, Universidad EAFIT.
