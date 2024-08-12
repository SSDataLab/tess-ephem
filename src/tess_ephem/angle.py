"""Helper functions to deal with spherical angles."""

from typing import Callable, Sequence

import numpy as np
from numpy import pi
from scipy.interpolate import CubicSpline


def enforce_positive_angle(angle: float) -> float:
    """Ensures `angle` is in the range [0, 2*pi[.

    Parameters
    ----------
    angle : float
        An angle given in radians.

    Returns
    -------
    angle : float
        Same angle guaranteed to be in the range [0, 2*pi[
    """
    return (angle + 2 * pi) % (2 * pi)


def create_angle_interpolator(
    x: Sequence[float], angle: Sequence[float], enforce_positive: bool = False
) -> Callable[[float], float]:
    """Returns a function to interpolate angles as a function of an independent variable.

    The function returned deals appropriately with the boundary at 360-0 degrees
    by interpolating the sine and cosine of the angle separately.

    Parameters
    ----------
    x : Sequence[float]
        Independent variable.
    angle : Sequence[float]
        Angle in degrees (i.e., the dependent variable)/
    enforce_positive : bool
        If `True`, ensures angle is in the range [0, 2*pi[

    Returns
    -------
    interpf : Callable
        Interpolation function
    """
    sin = np.sin(np.radians(angle))
    cos = np.cos(np.radians(angle))
    cs_sin = CubicSpline(x, sin)
    cs_cos = CubicSpline(x, cos)

    def interpf(x):
        if enforce_positive:
            return np.degrees(enforce_positive_angle(np.arctan2(cs_sin(x), cs_cos(x))))
        return np.degrees(np.arctan2(cs_sin(x), cs_cos(x)))

    return interpf
