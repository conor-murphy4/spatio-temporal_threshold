# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
import numpy as np
import theano.tensor

from .._bases.models import StressInvariantActivityRateModel, StressDependentActivityRateModel


class UniformRate(StressInvariantActivityRateModel):
    r"""Model with stress dependent activity rate.

    The variables controlling this model are:

    * ``Mmin`` (:math:`M_m`): Minimum magnitude.
    * ``beta0`` (:math:`\beta_0`, [log(m\ :sup:`-2` day\ :sup:`-1`\ )]): Background event rate exponent.

    The parameters of the observed distribution are computed following:

    .. math::
       \mathrm{d} \lambda(\mathbf{x}, t_i) &= e^{\beta_0} \left(t_{i + 1} - t_i\right) \mathrm{d}S, \\
       \Lambda_\mathrm{ETAS} &= 0,

    where :math:`t_i` is the :math:`i`-th epoch (the difference is computed in days) and :math:`\mathrm{d}S` is the
    grid node area (in meters squared). Note that :math:`\mathrm{d} \lambda(\mathbf{x}, t_i)` is set to zero for
    any epoch and location where there is no depletion (either because of lack or data or due to pressure recovery).

    The value of the :math:`i`-th event, used to calculate the log-likelihood, is computed using,

    .. math::
       \lambda(\mathbf{x}_i, t_i) = e^{\beta_0}
    """
    _vars = StressInvariantActivityRateModel._vars + ('beta0',)

    def __init__(self, thinsheet, label='uni.uni1'):
        super().__init__(thinsheet, label)

    def _observations(self, catalogue, var):
        b0 = var['beta0']
        return np.exp(b0) * np.ones(catalogue.nevents)

    def _generate_parameter(self, name, var, covs):
        if name == 'dlambda':
            b0 = var['beta0']
            uniform = covs['uniform_grid'][..., None]
            return np.exp(b0) * uniform * self._thinsheet.ds
        if name == 'etas':
            return 0
        return super()._generate_parameter(name, var, covs)


class ExponentialTrendCoulomb(StressDependentActivityRateModel):
    r"""Model with stress dependent activity rate.

    The variables controlling this model are:

    * ``Mmin`` (:math:`M_m`): Minimum magnitude.
    * ``beta0`` (:math:`\beta_0`, [log(m\ :sup:`-3`\ )]): Background event rate exponent.
    * ``beta1`` (:math:`\beta_1`, [MPa\ :sup:`-1`\ ]): Characteristic exponential increase in rate with Coulomb stress.
    * ``beta2`` (:math:`\beta_2`, [m]): Topographic gradient grid smoothing length-scale.
    * ``beta3`` (:math:`\beta_3`): Maximum throw-to-thickness ratio used when computing topographic gradient grid.
    * ``beta4`` (:math:`\beta_4`, [log\ :sub:`10`\ (MPa)]): Skeleton modulus exponent; :math:`H_s = 10^{\beta_4}`.

    The parameters of the observed distribution are computed following:

    .. math::
       \mathrm{d} \lambda(\mathbf{x}, t_i) &= h e^{\beta_0} \left(e^{\beta_1 \Delta C(\mathbf{x}, t_{i + 1})} -
       e^{\beta_1 \Delta C(\mathbf{x}, t_i)}\right) \mathrm{d}S, \\
       \Lambda_\mathrm{ETAS} &= 0,

    where :math:`\Delta C(\mathbf{x}, t_i)` is the incremental Coulomb stress grid on the :math:`i`-th epoch,
    :math:`\mathrm{d}S` is the grid node area (in meters squared) and :math:`h` is the thickness grid (in meters).

    The value of the :math:`i`-th event, used to calculate the log-likelihood, is computed using,

    .. math::
       \lambda(\mathbf{x}_i, t_i) = h\beta_1\left.\frac{\partial \Delta C}{\partial t}\right|_{\mathbf{x}_i, t_i}
       e^{\beta_0 + \beta_1 \Delta C(\mathbf{x}_i, t_i)}.
    """
    _vars = StressDependentActivityRateModel._vars + ('beta0', 'beta1')

    def __init__(self, thinsheet, label='ets.etc'):
        super().__init__(thinsheet, label)

    def _observations(self, catalogue, var):
        b0, b1 = var['beta0'], var['beta1']
        b2, b3, b4 = var['beta2'], var['beta3'], var['beta4']
        ethickness = self._thinsheet.thickness_event(catalogue.location)
        estress = self._thinsheet.stress_event(catalogue.time, catalogue.location, b2, b3, b4, nu=0)
        estress_dt = self._thinsheet.stress_event(catalogue.time, catalogue.location, b2, b3, b4, nu=1)
        return ethickness * b1 * estress_dt * np.exp(b0 + b1 * estress)

    def _generate_parameter(self, name, var, covs):
        if name == 'dlambda':
            b0, b1 = var['beta0'], var['beta1']
            stress, thickness = covs['stress_grid'][..., None], covs['thickness_grid'][..., None]
            return thickness * (np.exp(b0 + b1 * stress[1:]) - np.exp(b0 + b1 * stress[:-1])) * self._thinsheet.ds
        if name == 'etas':
            return 0
        return super()._generate_parameter(name, var, covs)


class ExponentialTrendCoulombEpidemicTypeAftershockSequence(ExponentialTrendCoulomb):
    r"""Model with stress dependent activity rate and epidemic type aftershock sequence (ETAS).

    The variables controlling this model are:

    * ``Mmin`` (:math:`M_m`): Minimum magnitude.
    * ``beta0`` (:math:`\beta_0`, [log(m\ :sup:`-3`\ )]): Background event rate exponent.
    * ``beta1`` (:math:`\beta_1`, [MPa\ :sup:`-1`\ ]): Characteristic exponential increase in rate with Coulomb stress.
    * ``beta2`` (:math:`\beta_2`, [m]): Topographic gradient grid smoothing length-scale.
    * ``beta3`` (:math:`\beta_3`): Maximum throw-to-thickness ratio used when computing topographic gradient grid.
    * ``beta4`` (:math:`\beta_4`, [log\ :sub:`10`\ (MPa)]): Skeleton modulus exponent; :math:`H_s = 10^{\beta_4}`.
    * ``etas_k`` (:math:`k`): Aftershock trigger function weight.
    * ``etas_a`` (:math:`a`): Aftershock productivity parameter.
    * ``etas_p`` (:math:`p`): Aftershock temporal clustering exponent.
    * ``etas_q`` (:math:`q`): Aftershock spatial clustering exponent.
    * ``etas_c`` (:math:`c`, [days]): Aftershock characteristic inter-event time.
    * ``etas_d`` (:math:`d`, [m\ :sup:`2`\ ]): Aftershock characteristic inter-event distance squared.

    The parameters of the observed distribution are computed following:

    .. math::
       \mathrm{d} \lambda(\mathbf{x}, t_i) &= h e^{\beta_0} \left(e^{\beta_1 \Delta C(\mathbf{x}, t_{i + 1})} -
       e^{\beta_1 \Delta C(\mathbf{x}, t_i)}\right) \mathrm{d}S, \\
       \Lambda_\mathrm{ETAS} &= k \sum_{j = 1}^{n_e} e^{a(M_j - M_m)},

    where :math:`\Delta C(\mathbf{x}, t_i)` is the incremental Coulomb stress grid on the :math:`i`-th epoch,
    :math:`\mathrm{d}S` is the grid node area (in meters squared), :math:`h` is the thickness grid (in meters),
    :math:`n_e` is the number of events and :math:`M_j` is the magnitude of the :math:`j`-th event.

    The value of the :math:`i`-th event, used to calculate the log-likelihood, is computed using,

    .. math::
       \lambda(\mathbf{x}_i, t_i) = h\beta_1\left.\frac{\partial \Delta C}{\partial t}\right|_{\mathbf{x}_i, t_i}
       e^{\beta_0 + \beta_1 \Delta C(\mathbf{x}_i, t_i)} + k \frac{p-1}{c} \frac{q-1}{\pi d}
       \sum_{j < i} e^{a(M_j - M_m)} \left(\frac{\Delta t_{ij}}{c} + 1\right)^{-p}
       \left(\frac{\Delta\mathbf{x}_{ij}^2}{d} + 1\right)^{-q}.
    """
    _vars = ExponentialTrendCoulomb._vars + ('etas_k', 'etas_a', 'etas_p', 'etas_q', 'etas_c', 'etas_d')

    def __init__(self, thinsheet, label='ets.etc.etas'):
        super().__init__(thinsheet, label)

    @staticmethod
    def __trigger(catalogue, var):
        dm = catalogue.magnitude - var['Mmin']
        dt = catalogue.dt
        dr2 = catalogue.dr ** 2

        p, q, c, d = var['etas_p'], var['etas_q'], var['etas_c'], var['etas_d']
        temporal = (p - 1.) * c ** (p - 1) * theano.tensor.switch(dt > 0, (dt + c) ** -p, 0)
        spatial = (q - 1.) * d ** (q - 1) * theano.tensor.switch(dt > 0, (dr2 + d) ** -q, 0) / np.pi

        k, a = var['etas_k'], var['etas_a']
        return k * (np.exp(a * dm) * spatial * temporal).sum(axis=-1)

    def _observations(self, catalogue, var):
        return super()._observations(catalogue, var) + self.__trigger(catalogue, var)

    def _generate_parameter(self, name, var, covs):
        if name == 'etas':
            dm = covs['magnitude_event'][..., None] - var['Mmin']
            k, a = var['etas_k'], var['etas_a']
            return k * np.exp(a * dm).sum(axis=0)
        return super()._generate_parameter(name, var, covs)
