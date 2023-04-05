# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
"""Distributions available for models (based on PyMC3 distributions).

Distributions present:

* :class:`~seismod.distributions.PowerLawExponentialTaper`
* :class:`~seismod.distributions.SpaceTime`
"""
import numpy as np
import scipy.special

from . import utils
from ._bases.distributions import BaseDistribution


class PowerLawExponentialTaper(BaseDistribution):
    r"""Power law with exponential taper.

    The probability density function of this distribution is:

    .. math::
       f(s \mid \beta, \zeta, s_m) =
           \frac{1}{s_m}\left(\beta + \zeta \frac{s}{s_m}\right) \left(\frac{s}{s_m}\right)^{-\beta - 1}
           e^{-\zeta (s/s_m - 1)}

    with support over :math:`s \in [s_m, \infty)`, where :math:`s_m` corresponds to the minimum seismic moment.

    **Important:** The values for :math:`s` and :math:`s_m` must always  be given in terms of magnitude rather
    than seismic moment. The transformation is handled internally using :func:`seismod.utils.mag2mom`.
    """
    parameters = ('beta', 'zeta', 'Mmin')

    @staticmethod
    def _extract_parameters(**kwargs):
        beta, zeta, mmin = kwargs.pop('beta'), kwargs.pop('zeta'), kwargs.pop('Mmin')
        if kwargs:
            raise ValueError('Unexpected arguments: {}'.format(', '.join(kwargs)))
        return beta, zeta, mmin

    @classmethod
    def rvs(cls, size=None, rng=None, **kwargs):
        r"""Draw random magnitude values from the distribution.

        The parameters can be scalars or arrays, in which case the sampling will be performed for each parameter
        in the array (as long as all parameters are broadcastable with each other).

        Uses the same probability function defined above, unless a maximum magnitude is provided. If that's the case
        then the distribution is constrained according based on the exceedance probability:

        .. math::
           \mathrm{Pr}(\geq s \mid s_m \leq s \leq s_M) =
               \frac{P(\geq s \mid s \geq s_m) - P(\geq s_M \mid s_M \geq s_m)}{1 - P(\geq s_M \mid s_M \geq s_m)},

        where :math:`s_M` is the maximum seismic moment and :math:`P(\geq s \mid s \geq s_m)` is the survival function.

        Args:
            size (int, tuple[int], None): Number of events (returns scalar if omitted)
            rng (np.random.Generator, None): Random number generator used to generate the samples

        Keyword Args:
            beta (float, np.ndarray[float]): :math:`\beta` value(s) of the distribution
            zeta (float, np.ndarray[float]): :math:`\zeta` value(s) of the distribution
            Mmin (float, np.ndarray[float]): Minimum magnitude value(s) of the distribution
            Mmax (float, np.ndarray[float]): Maximum magnitude possible (optional)

        Returns:
            float, np.ndarray[float]: Randomly sampled magnitude values
        """
        mmax = kwargs.pop('Mmax', None)
        beta, zeta, mmin = cls._extract_parameters(**kwargs)
        exc_prob = 1 - rng.random(size=size) if rng else 1 - np.random.random(size=size)
        smin = utils.mag2mom(mmin)
        if mmax is not None:
            with np.errstate(under='ignore'):
                smax = utils.mag2mom(mmax)
                max_unconstrained = (smax / smin) ** -beta * np.exp(-zeta * (smax / smin - 1))
                exc_prob += max_unconstrained * (1 - exc_prob)
        sratio = np.power(exc_prob, -1 / beta)
        if zeta != 0:
            zbratio = zeta / beta
            sratio = scipy.special.lambertw(zbratio * np.exp(zbratio) * sratio).real / zbratio
        return utils.mom2mag(smin * sratio)

    @classmethod
    def loglike(cls, value, **kwargs):
        r"""Computes the log-likelihood of the distribution for the given magnitude values.

        For :math:`n_e` events, this is calculated using,

        .. math::
           \ell = \sum_{i=1}^{n_e} \left(\log\left(\beta_i + \zeta_i \frac{s_i}{s_m}\right) - (1 + \beta_i)
               \log\left(\frac{s_i}{s_m}\right) - \zeta_i \left(\frac{s_i}{s_m} - 1\right) - \log(s_m)\right)

        Args:
            value (float, np.ndarray[float], theano.tensor.TensorVariable): Magnitude values
            **kwargs: Distribution parameters at which log-likelihood will be calculated.

        Returns:
            float, np.ndarray[float], theano.tensor.TensorVariable: Total log-likelihood
        """
        beta, zeta, mmin = cls._extract_parameters(**kwargs)
        s, smin = utils.mag2mom(value), utils.mag2mom(mmin)
        logpi = np.log(beta + zeta * s / smin) - (1 + beta) * np.log(s / smin) - zeta * (s / smin - 1) - np.log(smin)
        return logpi.sum(axis=0)

    @classmethod
    def sf(cls, value, **kwargs):
        r"""Computes the survival function (exceedance probability) for the given magnitude values.

        The function is given by,

        .. math::
           P(\geq s \mid s \geq s_m) = \left(\frac{s}{s_m}\right)^{-\beta} e^{-\zeta\left(\frac{s}{s_m} -1\right)}.

        Args:
            value (float, np.ndarray[float], theano.tensor.TensorVariable): Magnitude values
            **kwargs: Distribution parameters at which log-likelihood will be calculated.

        Returns:
            float, np.ndarray[float], theano.tensor.TensorVariable: Total log-likelihood
        """
        beta, zeta, mmin = cls._extract_parameters(**kwargs)
        s, smin = utils.mag2mom(value), utils.mag2mom(mmin)
        with np.errstate(under='ignore'):
            return (s / smin) ** -beta * np.exp(-zeta * (s / smin - 1))


class SpaceTime(BaseDistribution):
    r"""Space-time dependent Poisson process distribution.

    The probability density function of this distribution is:

    .. math::
       f\left(\lambda(\mathbf{x}, t) \mid \mathrm{d}\lambda(\mathbf{x}', t_i), \Lambda_\mathrm{ETAS}\right) =
           \lambda(\mathbf{x}, t) \exp\left(-\sum_{\mathbf{x}'}\sum_{i=0}^{n_t - 1}\mathrm{d}\lambda(\mathbf{x}', t_i)
           -\Lambda_\mathrm{ETAS}\right),

    with support over :math:`\lambda(\mathbf{x}, t) \in (0, \infty)`, where :math:`n_t` is the number of epochs.
    """
    parameters = ('dlambda', 'etas')

    @staticmethod
    def _extract_parameters(**kwargs):
        dlambda, etas = kwargs.pop('dlambda'), kwargs.pop('etas', 0)
        if kwargs:
            raise ValueError('Unexpected arguments: {}'.format(', '.join(kwargs)))
        return dlambda, etas

    @classmethod
    def rvs(cls, size=None, rng=None, **kwargs):
        r"""Draw random event locations and times (as indices) from the distribution.

        Only the parameter ``dlambda`` is required (``etas`` will be ignored if provided).

        The probability density function for each epoch and grid node is given by,

        .. math::
           f(\mathrm{x}, t_i) = \mathrm{d}\lambda(\mathbf{x}, t_i).

        The distribution *within* each epoch and grid node is considered to be uniform.

        Args:
            size (int, tuple[int], None): Number of events (returns scalar if omitted)
            rng (np.random.Generator, None): Random number generator used to generate the samples

        Keyword Args:
            dlambda (np.ndarray[float]): :math:`\mathrm{d}\lambda(\mathbf{x}, t_i)`. Shape ``(nt, nx, ny)``

        Returns:
            tuple[float, np.ndarray[float]]: Time, Easting and Northing grid indices for the sampled events. The
            indices themselves are floats, representing the fractional displacement within the index (e.g., an index
            of 4.25, corresponds to the fourth index and a factor of 0.25 towards the fifth)
        """
        dlambda, _ = cls._extract_parameters(**kwargs)
        try:
            size = tuple(size)
        except TypeError:
            size = (size,) if size else ()
        nt, nx, ny, nv = dlambda.shape
        itime, ix, iy = np.empty(size + (nv,)), np.empty(size + (nv,)), np.empty(size + (nv,))
        for ic in range(nv):
            density = dlambda[..., ic].ravel() / np.sum(dlambda[..., ic])
            i = rng.choice(nt * nx * ny, size, p=density) if rng else np.random.choice(nt * nx * ny, size, p=density)
            i = i.astype('f8')
            itime[..., ic], ix[..., ic], iy[..., ic] = i // (nx * ny), i % (nx * ny) // ny, i % (nx * ny) % ny
            itime[..., ic] += rng.random(size=size) if rng else np.random.random(size=size)
            ix[..., ic] += rng.uniform(-.5, .5, size=size) if rng else np.random.uniform(-.5, .5, size=size)
            iy[..., ic] += rng.uniform(-.5, .5, size=size) if rng else np.random.uniform(-.5, .5, size=size)
        if nv == 1:
            return itime.squeeze(axis=-1), ix.squeeze(axis=-1), iy.squeeze(axis=-1)
        return itime, ix, iy

    @classmethod
    def loglike(cls, value, **kwargs):
        r"""Computes the log-likelihood of the distribution for the given :math:`\lambda` values.

        For :math:`n_e` events, this is calculated using,

        .. math::
           \ell = \sum_{i=1}^{n_e} \log\left(\lambda(\mathbf{x}_i, t_i)\right) -
               \sum_{\mathbf{x}}\sum_{j=0}^{n_t - 1}\mathrm{d}\lambda(\mathbf{x}, t_j) -\Lambda_\mathrm{ETAS}

        Args:
            value (float, np.ndarray[float], theano.tensor.TensorVariable): Magnitude values
            **kwargs: Distribution parameters at which log-likelihood will be calculated.

        Returns:
            float, np.ndarray[float], theano.tensor.TensorVariable: Total log-likelihood
        """
        dlambda, etas = cls._extract_parameters(**kwargs)
        return np.log(value[value > 0]).sum(axis=0) - (dlambda.sum(axis=(0, 1, 2)) + etas).squeeze()

    @classmethod
    def nevents(cls, rng=None, **kwargs):
        r"""Randomly draw the number of events from the distribution.

        Only the parameter ``dlambda`` is required (``etas`` will be ignored if provided). The number of events
        is drawn from a Poisson distribution with mean:

        .. math::
           \Lambda = \sum_{\mathbf{x}}\sum_{i=0}^{n_t - 1}\mathrm{d}\lambda(\mathbf{x}, t_i).

        Args:
            rng (np.random.Generator, None): Random number generator used for the draw

        Keyword Args:
            dlambda (np.ndarray[float]): :math:`\mathrm{d}\lambda(\mathbf{x}, t_i)`. Shape ``(nt, nx, ny)``

        Returns:
            int: Number of events
        """
        dlambda, _ = cls._extract_parameters(**kwargs)
        return rng.poisson(dlambda.sum(axis=(0, 1, 2))) if rng else np.random.poisson(dlambda.sum(axis=(0, 1, 2)))
