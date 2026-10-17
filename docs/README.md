# Building the documentation locally

This guide uses the shared development environment at `~/.venvs/unsatstab-uv`. Activate it before running project commands. `uv sync --active` installs the locked documentation dependencies into the active environment and does not create a `.venv` in the repository.

From the repository root:

```bash
source "$HOME/.venvs/unsatstab-uv/bin/activate"
uv sync --active --locked --group docs
```

## Live preview

```bash
uv run --active --locked mkdocs serve
```

Open the local address printed by MkDocs, usually <http://127.0.0.1:8000/>. The preview rebuilds when documentation sources change. Stop it with `Ctrl+C`.

## Build and serve the static site

Build the HTML into `site/`:

```bash
uv run --active --locked mkdocs build --strict
```

To serve that built output locally, start a static file server from the repository root:

```bash
python -m http.server 8000 --directory site
```

Then open <http://127.0.0.1:8000/>. Stop the server with `Ctrl+C`. The `site/` directory is generated output; edit the Markdown, notebooks, and assets under `docs/` instead.

## Publish the documentation on GitHub Pages

MkDocs builds the HTML; GitHub Pages hosts it. The repository's
`.github/workflows/docs.yml` connects these steps. It builds with
`mkdocs build --strict` and deploys the generated `site/` directory whenever
changes reach `main`. Pull requests build the documentation without deploying it.

Configure the repository once under **Settings → Pages → Build and deployment**:
select **GitHub Actions** as the source. Then edit the documentation, rerun any
changed notebooks to save their outputs, inspect the local preview, and commit
and push the source changes to `main`. The **Documentation** workflow publishes
the updated site at <https://eamontoyaa.github.io/unsatstab/>. Check its status in
the repository's **Actions** tab. It can also be started manually with
**Run workflow**; select `main` for the published documentation.

Do not commit `site/`; it is generated during the workflow. This configuration
uses an Actions artifact deployment rather than a separate `gh-pages` branch,
so it does not need `mkdocs gh-deploy`.

Package releases and documentation deployment are independent. Publishing a
GitHub release does not upload the package to PyPI; package publication is manual.

## Writing notebooks

Put parameter values directly in constructors when they are not reused. For later plot labels or checks, read geometry and material values from the configured objects instead of duplicating them in separate variables. Keep named histories, shared settings, and derived quantities when they are reused across experiments.

Comment every unit conversion with its source and destination units, including conversions for plots and summary tables. Keep rainfall intensities and durations explicit, and verify the integrated amount when comparing histories with equal totals. After changing an experiment, rerun it and update any numerical interpretation in Markdown cells.

Number notebooks independently within `tutorials/`, `verification/`, and `validation/` using `01_`, `02_`, and so on, following their order in `mkdocs.yml`. Update navigation and cross-references when renaming a notebook.

Embed published reference figures in dedicated code cells with `IPython.display.Image` and a caption, immediately after the matching computed figure. Use the same display width and panel order for aligned comparisons. Add the `remove-input` cell tag and `jupyter: {source_hidden: true}` metadata: MkDocs hides the source while retaining outputs, and JupyterLab starts with the source collapsed. Other notebook frontends may handle collapsed input differently. Execute the notebook and save its image outputs before building the documentation.

Generated plots are stored in notebook outputs and rendered from those outputs in the documentation. Do not export duplicate figure files with `savefig`. Keep published reference images and branding assets in `docs/assets/`; notebooks may read those files for comparison.

## Writing equations

Use `$...$` for inline mathematics and `$$...$$` on separate lines for display equations, in both Markdown files and notebook Markdown cells. Leave a blank line before and after a display block:

```markdown
The drainage efficiency is $\eta$.

$$
q_{\mathrm{c}} = \eta K_{\mathrm{s}}\sin^2\delta.
$$
```

This convention works in Jupyter and the notebook-to-HTML conversion used by the site. `assets/mathjax.js` enables the same delimiters and processes both MkDocs `arithmatex` wrappers and notebook Markdown containers. Keep code examples and Matplotlib label strings unchanged; they are not documentation Markdown equations.

Avoid `ignoreHtmlClass: ".*"`: it suppresses nested notebook paragraphs even when the outer Markdown container is enabled. Ignore only explicitly marked `mathjax_ignore` or `tex2jax_ignore` elements. Verify browser-rendered equations after changing this configuration; a successful HTML build alone does not establish successful typesetting.

Use upright descriptive subscripts, for example `$\theta_{\mathrm{s}}$`, `$K_{\mathrm{s}}$`, `$q_{\mathrm{b}}$` and `$I_{\mathrm{ZLT}}$`. Write the safety factor as `$\mathrm{FS}$`. Keep variable subscripts and iteration indices italic, for example `$S_i$`, `$t_n$` and `$r_m$`. Put verified source equation numbers and implementation notes in Markdown footnotes; the theme displays their superscripts with hover previews. Keep explanations needed to understand a formula in the main text. Distinguish published relationships, algebraic rewrites, and package definitions; use the [equation provenance table](fundamentals.md#equation-provenance).
