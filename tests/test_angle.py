import numpy as np
import pytest
from numpy import pi
from numpy.testing import assert_almost_equal, assert_array_equal

from tess_ephem.angle import create_angle_interpolator, enforce_positive_angle


@pytest.mark.parametrize(
    "angle,expected",
    [(0, 0), (-pi, pi), (2 * pi, 0), (-2 * pi, 0), (3 * pi, pi), (-3 * pi, pi)],
)
def test_enforce_positive_angle(angle, expected):
    """Does `enforce_positive_angle` return the expected result?"""
    assert enforce_positive_angle(angle) == expected


def test_enforce_positive_angle_vectorized():
    """Is `enforce_positive_angle` vectorized?"""
    assert all(
        enforce_positive_angle(np.array([0, 2 * pi, 4 * pi])) == np.array([0, 0, 0])
    )


def test_angle_interpolator():
    x = [-10, 10]
    angle = [-10 * pi, 10 * pi]
    f = create_angle_interpolator(x, angle)
    testvalues = np.array([-10, 0, 10])
    assert_array_equal(f(testvalues), testvalues * pi)


@pytest.mark.parametrize("x,angle", [([1, 2, 3], [0, 0, 0]), ([1, 2, 3], [1, 2, 3])])
def test_angle_interpolator_roundtrip(x, angle):
    f = create_angle_interpolator(x, angle)
    assert_almost_equal(f(x), angle)
