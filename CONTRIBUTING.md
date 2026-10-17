# Contributing

Use the shared `~/.venvs/unsatstab-env` development environment. For the local documentation preview and static site commands, see [docs/README.md](docs/README.md). Activate the environment before syncing development tools:

```bash
source "$HOME/.venvs/unsatstab-env/bin/activate"
uv sync --active --locked --all-groups
```

The public source contains the `unsatstab` model, spatial slope-unit support, tests, and English documentation notebooks. Keep changes to the equations, assumptions, units, and numerical controls explicit. For model changes, add a regression case for meaningful behavior and discuss effects on the Qva/Qtb benchmark comparisons.

Put executable user examples in Jupyter notebooks under `docs/tutorials/`, reference comparisons under `docs/verification/`, and study-site applications under `docs/validation/`. Keep their input data alongside the corresponding notebooks. Unpublished manuscripts remain under `research/` and are excluded from the package distributions; the wheel contains only the Python implementation.

Before a release, build the package distributions with `uv build`, inspect their contents, and build the HTML documentation with `uv run --active --locked mkdocs build --strict`.

Package publication to PyPI is manual; there is no automatic release-upload workflow. Documentation deployment remains automatic on updates to `main`; see [the GitHub Pages setup](docs/README.md#publish-the-documentation-on-github-pages).
