"""Check metadata restrictions enforced by the PyPI upload API."""

from pathlib import Path
import tomllib


def test_project_url_labels_fit_pypi_limit():
    project = tomllib.loads(
        (Path(__file__).resolve().parents[1] / "pyproject.toml").read_text()
    )["project"]
    for label in project.get("urls", {}):
        assert len(label) <= 32, f"PyPI rejects Project-URL labels longer than 32 characters: {label!r}"
