import unittest
from unittest import mock

import numpy as np

from seismod import utils, distributions


class TestPowerLawExponentialTaperRandomValueSampling(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.random_caller = mock.patch('numpy.random.random', autospec=True)

    def test_when_sampling_without_size_then_sample_is_scalar(self):
        beta, zeta, mmin = .7, 1e-5, 1.5
        distr = distributions.PowerLawExponentialTaper.dist(beta=beta, zeta=zeta, Mmin=mmin)
        sample = distr.rvs(beta=beta, zeta=zeta, Mmin=mmin)
        self.assertTupleEqual(sample.shape, ())

    def test_when_sampling_with_size_then_sample_has_corresponding_shape(self):
        beta, zeta, mmin = .7, 1e-5, 1.5
        distr = distributions.PowerLawExponentialTaper.dist(beta=beta, zeta=zeta, Mmin=mmin)
        sample = distr.rvs(size=(2, 3, 5), beta=beta, zeta=zeta, Mmin=mmin)
        self.assertTupleEqual(sample.shape, (2, 3, 5))

    def test_when_generating_rvs_with_parameters_then_parameter_values_are_used_for_calculation(self):
        beta, zeta = .7, 1e-5
        mmin, mmax = 1.5, 7
        expected_mags = np.linspace(mmin, mmax * .8, 100, endpoint=False)
        sratio = utils.mag2mom(expected_mags) / utils.mag2mom(mmin)
        expected_exc_prob = sratio ** -beta * np.exp(-zeta * (sratio - 1))
        with self.random_caller as mock_random:
            mock_random.return_value = 1 - expected_exc_prob
            catalogue = distributions.PowerLawExponentialTaper.rvs(beta=beta, zeta=zeta, Mmin=mmin)
        np.testing.assert_allclose(catalogue, expected_mags)


class TestPowerLawExponentialTaperLogLikelihood(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.random_caller = mock.patch('numpy.random.random', autospec=True)

    def test_when_computing_logp_without_catalogues_then_logp_has_no_shape(self):
        beta, zeta = .7, 1e-5
        mmin, mmax = 1.5, 7
        distr = distributions.PowerLawExponentialTaper.dist(beta=beta, zeta=zeta, Mmin=mmin)
        events = np.random.random(50) * (mmax - mmin) + mmin
        logp = distr.logp(events)
        self.assertTupleEqual(logp.eval().shape, ())

    def test_when_computing_logp_with_catalogues_then_logp_shape_matches_catalogues(self):
        beta, zeta = .7, 1e-5
        mmin, mmax = 1.5, 7
        distr = distributions.PowerLawExponentialTaper.dist(beta=beta, zeta=zeta, Mmin=mmin)
        events = np.random.random((50, 10)) * (mmax - mmin) + mmin
        logp = distr.logp(events)
        self.assertTupleEqual(logp.eval().shape, (10,))

    def test_when_calculating_likelihood_with_logp_then_results_are_equal_to_calculation_with_loglike(self):
        distr = distributions.PowerLawExponentialTaper.dist(beta=.7, zeta=1e-2, Mmin=1.5)
        with self.random_caller as mock_random:
            mock_random.return_value = np.linspace(0, 1, 1000, endpoint=False)
            catalogue = distr.rvs(beta=.7, zeta=1e-2, Mmin=1.5)
            from_logp = distr.logp(catalogue).eval()
            from_loglike = distr.loglike(catalogue, beta=.7, zeta=1e-2, Mmin=1.5).sum()
        self.assertAlmostEqual(from_loglike, from_logp)


if __name__ == '__main__':
    unittest.main()
