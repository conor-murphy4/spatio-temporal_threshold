import unittest

import numpy as np

from seismod import utils


class TestMomentMagnitudeTransformation(unittest.TestCase):
    def test_when_running_magnitude_transformation_with_defaults_in_series_then_output_is_equal_to_input(self):
        mags = np.linspace(1, 9, 100)
        out = utils.mom2mag(utils.mag2mom(mags))
        np.testing.assert_allclose(out, mags)

    def test_when_running_magnitude_transformation_with_given_parameters_in_series_then_output_is_equal_to_input(self):
        mags = np.linspace(1, 9, 100)
        out = utils.mom2mag(utils.mag2mom(mags, d=1, c=1), d=1, c=1)
        np.testing.assert_allclose(out, mags)

    def test_when_running_moment_transformation_with_defaults_in_series_then_output_is_equal_to_input(self):
        moms = np.linspace(1e10, 1e20, 100)
        out = utils.mag2mom(utils.mom2mag(moms))
        np.testing.assert_allclose(out, moms)

    def test_when_running_moment_transformation_with_given_parameters_in_series_then_output_is_equal_to_input(self):
        moms = np.linspace(1e10, 1e20, 100)
        out = utils.mag2mom(utils.mom2mag(moms, d=1, c=1), d=1, c=1)
        np.testing.assert_allclose(out, moms)


if __name__ == '__main__':
    unittest.main()
