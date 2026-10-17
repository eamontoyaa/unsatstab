# API reference

Functions and classes are documented by their defining module. Common entry points are also available from `import unsatstab as us`.

- [Parameters (`unsatstab.parameters`)](api/parameters.md): Hydraulic and strength properties, initial conditions, and model resolution.
- [Rainfall (`unsatstab.rainfall`)](api/rainfall.md): Piecewise-constant rainfall forcing and optional antecedent histories.
- [Basal boundaries (`unsatstab.boundaries`)](api/boundaries.md): Drainage capacity and basal boundary assumptions.
- [Column solver (`unsatstab.solver`)](api/solver.md): Single-column infiltration, storage, and stability simulation.
- [Column results (`unsatstab.results`)](api/results.md): Time series, depth profiles, and exported column results.
- [Spatial slope units (`unsatstab.slope_units`)](api/slope_units.md): Spatial schemas, input readers, distributed simulations, and result exports.
- [Analytical flow (`unsatstab.analytical`)](api/analytical.md): Gardner relations and the analytical unsaturated-flow response.
- [Characteristic roots (`unsatstab.roots`)](api/roots.md): Root equations and numerical solution of the characteristic spectrum.
- [Water table (`unsatstab.water_table`)](api/water_table.md): Fillable pore volume, retained storage, and water-table rise.
- [Slope stability (`unsatstab.stability`)](api/stability.md): Infinite-slope factor of safety and effective-stress calculations.

Verified source equation numbers and implementation adaptations are summarized in [equation provenance](fundamentals.md#equation-provenance).
