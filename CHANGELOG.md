# Changelog

## 1.0.0

- Preserve Colab's NumPy 2.1.3 when it satisfies the model, and keep notebook-only plotting dependencies out of the runtime installation.

- Provide analytical unsaturated infiltration with the default Srivastava–Yeh constitutive conversion and an explicit Baum et al. comparison set.
- Compute retained water, estimated water-table rise, pressure, and infinite-slope stability throughout each column.
- Separate maximum column depth from the optional observation plane; initialize a deeper requested water table at a finite base with zero initial water-column height.
- Use hours for public time inputs, outputs, rainfall breakpoints, and storage integration controls.
- Support antecedent rainfall histories with continuous state and optional stage-specific drainage.
- Resolve spatial columns from material mappings, optional local property overrides, and intensity–duration rainfall histories.
- Include three tutorials and four numerical-verification notebooks, plus a separate study-site application, with documented assumptions, equation provenance, and independent reference comparisons.
- Standardize documentation and figure styles; simplify the repository README and remove superseded guides and unused integration code.
- Accept GIS column mapping while retaining original attributes in spatial exports; compare both equation sets in a dedicated notebook.
- Use per-unit mean geometry and two geozones from SU_tutorial.gpkg in the spatial tutorial; relocate the previous example to the Ruiz et al. study-site application.

This version prepares the release artifacts; publication to PyPI is a separate step.
