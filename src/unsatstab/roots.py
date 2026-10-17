"""Characteristic roots for the analytical unsaturated solution."""

from __future__ import annotations

from math import pi, sin, cos
import numpy as np
from scipy.optimize import brentq


def characteristic_equation(root: float, alpha_depth: float) -> float:
    """Evaluate the pole-free rewrite of Baum et al. (2010), Eq. (6):
    sin(alpha_depth*r) + 2*r*cos(alpha_depth*r), where alpha_depth = a*z_w0.
    """
    return sin(alpha_depth * root) + 2.0 * root * cos(alpha_depth * root)


def compute_roots(
    transformed_alpha: float,
    depth: float,
    n_roots: int,
    *,
    eps: float = 1.0e-10,
) -> np.ndarray:
    """Compute characteristic roots for the analytical series."""

    if transformed_alpha <= 0.0:
        raise ValueError("transformed_alpha must be positive.")
    if depth <= 0.0:
        raise ValueError("depth must be positive.")
    if n_roots <= 0:
        raise ValueError("n_roots must be positive.")

    alpha_depth = transformed_alpha * depth
    roots = np.empty(n_roots, dtype=float)

    for n in range(1, n_roots + 1):
        lo = (n - 1) * pi / alpha_depth + eps
        hi = n * pi / alpha_depth - eps
        roots[n - 1] = brentq(
            characteristic_equation,
            lo,
            hi,
            args=(alpha_depth,),
            xtol=eps,
            rtol=eps,
        )

    return roots
