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

from .. import utils
from .._bases.models import StressInvariantMagnitudeModel, StressDependentMagnitudeModel


class UniformBetaZeta(StressInvariantMagnitudeModel):
    r"""Stress independent magnitude model.

    The moment distribution is a power-law with an exponential taper. The variables controlling this model are:

    * ``Mmin`` (:math:`M_m`): Minimum magnitude.
    * ``beta`` (:math:`\beta`): Power-law exponent.
    * ``zeta`` (:math:`\zeta`): Exponential taper characteristic exponent.

    The variables match exactly with the observed distribution parameters.
    """
    _vars = StressInvariantMagnitudeModel._vars + ('beta', 'zeta')

    def __init__(self, thinsheet, label='uni.beta.zeta'):
        super().__init__(thinsheet, label)

    def _generate_parameter(self, name, var, covs):
        if name == 'Mmin':
            return var[name]
        if name == 'beta':
            return var[name]
        if name == 'zeta':
            return var[name]
        return super()._generate_parameter(name, var, covs)


class UniformBetaCornerMagnitude(StressInvariantMagnitudeModel):
    r"""Stress independent magnitude model.

    The moment distribution is a power-law with an exponential taper. The variables controlling this model are:

    * ``Mmin`` (:math:`M_m`): Minimum magnitude.
    * ``beta`` (:math:`\beta`): Power-law exponent.
    * ``Mc`` (:math:`M_c`): Corner magnitude.

    The first two variables map match exactly the corresponding observed distribution parameters. The last one
    is obtained using:

    .. math::
       \zeta = \frac{\mathcal{M}_m}{\mathcal{M}_c},

    where :math:`\mathcal{M}_m` and :math:`\mathcal{M}_c` are the seismic moments corresponding to :math:`M_m` and
    :math:`M_c`, respectively.
    """
    _vars = StressInvariantMagnitudeModel._vars + ('beta', 'Mc')

    def __init__(self, thinsheet, label='uni.beta.Mc'):
        super().__init__(thinsheet, label)

    def _generate_parameter(self, name, var, covs):
        if name == 'Mmin':
            return var[name]
        if name == 'beta':
            return var[name]
        if name == 'zeta':
            return utils.mag2mom(var['Mmin']) / utils.mag2mom(var['Mc'])
        return super()._generate_parameter(name, var, covs)


class UniformBetaCornerMoment(StressInvariantMagnitudeModel):
    r"""Stress independent magnitude model.

    The moment distribution is a power-law with an exponential taper. The variables controlling this model are:

    * ``Mmin`` (:math:`M_m`): Minimum magnitude.
    * ``beta`` (:math:`\beta`): Power-law exponent.
    * ``Sc`` (:math:`\mathcal{M}_c`, [Nm]): Corner seismic moment.

    The first two variables map match exactly the corresponding observed distribution parameters. The last one
    is obtained using:

    .. math::
       \zeta = \frac{\mathcal{M}_m}{\mathcal{M}_c},

    where :math:`\mathcal{M}_m` is the seismic moments corresponding to :math:`M_m`.
    """
    _vars = StressInvariantMagnitudeModel._vars + ('beta', 'Sc')

    def __init__(self, thinsheet, label='uni.beta.Sc'):
        super().__init__(thinsheet, label)

    def _generate_parameter(self, name, var, covs):
        if name == 'Mmin':
            return var[name]
        if name == 'beta':
            return var[name]
        if name == 'zeta':
            return utils.mag2mom(var['Mmin']) / var['Sc']
        return super()._generate_parameter(name, var, covs)


class InversePowerLawBeta(StressDependentMagnitudeModel):
    r"""Model with stress dependent :math:`\beta`.

    The moment distribution is a power-law. The variables controlling this model are:

    * ``Mmin`` (:math:`M_m`): Minimum magnitude.
    * ``beta2`` (:math:`\beta_2`, [m]): Topographic gradient grid smoothing length-scale.
    * ``beta3`` (:math:`\beta_3`): Maximum throw-to-thickness ratio used when computing topographic gradient grid.
    * ``beta4`` (:math:`\beta_4`, [log\ :sub:`10`\ (MPa)]): Skeleton modulus exponent; :math:`H_s = 10^{\beta_4}`.
    * ``theta0`` (:math:`\theta_0`): Minimum value for :math:`\beta`.
    * ``theta1`` (:math:`\theta_1`, [MPa]): Coulomb stress at which :math:`\beta` becomes stress dependent.
    * ``theta2`` (:math:`\theta_2`, [MPa]): Stress normalization factor.
    * ``theta3`` (:math:`\theta_3`): Power-law exponent on stress dependency.

    The remaining parameters of the observed distribution are computed following:

    .. math::
       \beta_i &= \theta_0 + \left(\frac{\Delta C_i - \theta_1}{\theta_2}\right)^{\theta_3}, \\
       \zeta &= 0,

    where :math:`\Delta C_i` is the incremental Coulomb stress for the :math:`i`-th event.

    If the computed :math:`\beta` is higher than one, it will be set to one.

    If :math:`\Delta C_i \leq \theta_1`, :math:`\beta` will be set to :math:`\theta_0`.
    """
    _vars = StressDependentMagnitudeModel._vars + ('theta0', 'theta1', 'theta2', 'theta3')

    def __init__(self, thinsheet, label='ets.ipc'):
        super().__init__(thinsheet, label)

    def _generate_parameter(self, name, var, covs):
        if name == 'Mmin':
            return var[name]
        if name == 'beta':
            t0, t1, t2, t3 = var['theta0'], var['theta1'], var['theta2'], var['theta3']
            estress = covs['stress_event']
            return theano.tensor.minimum(t0 + theano.tensor.switch(estress > t1, ((estress - t1) / t2) ** -t3, 0.), 1.)
        if name == 'zeta':
            return 0
        return super()._generate_parameter(name, var, covs)


class HyperbolicTangentBeta(StressDependentMagnitudeModel):
    r"""Model with stress dependent :math:`\beta`.

    The moment distribution is a power-law. The variables controlling this model are:

    * ``Mmin`` (:math:`M_m`): Minimum magnitude.
    * ``beta2`` (:math:`\beta_2`, [m]): Topographic gradient grid smoothing length-scale.
    * ``beta3`` (:math:`\beta_3`): Maximum throw-to-thickness ratio used when computing topographic gradient grid.
    * ``beta4`` (:math:`\beta_4`, [log\ :sub:`10`\ (MPa)]): Skeleton modulus exponent; :math:`H_s = 10^{\beta_4}`.
    * ``theta0`` (:math:`\theta_0`): Base value for :math:`\beta`.
    * ``theta1`` (:math:`\theta_1`): Weight factor for hyperbolic tangent contribution.
    * ``theta2`` (:math:`\theta_2`, [MPa\ :sup:`-1`\ ]): Stress variability rate for :math:`\beta`.
    * ``theta3`` (:math:`\theta_3`): Normalized stress offset.

    The remaining parameters of the observed distribution are computed following:

    .. math::
       \beta_i &= \theta_0 + \theta_1 (1 - \tanh(\theta_2 \Delta C_i - \theta_3)), \\
       \zeta &= 0,

    where :math:`\Delta C_i` is the incremental Coulomb stress for the :math:`i`-th event.
    """
    _vars = StressDependentMagnitudeModel._vars + ('theta0', 'theta1', 'theta2', 'theta3')

    def __init__(self, thinsheet, label='ets.htb'):
        super().__init__(thinsheet, label)

    def _generate_parameter(self, name, var, covs):
        if name == 'Mmin':
            return var[name]
        if name == 'beta':
            t0, t1, t2, t3 = var['theta0'], var['theta1'], var['theta2'], var['theta3']
            estress = covs['stress_event']
            return t0 + t1 * (1 - np.tanh(t2 * estress - t3))
        if name == 'zeta':
            return 0
        return super()._generate_parameter(name, var, covs)


class CriticalPointScalingZeta(StressDependentMagnitudeModel):
    r"""Model with stress dependent :math:`\zeta`.

    The moment distribution is a power-law with an exponential taper. The variables controlling this model are:

    * ``Mmin`` (:math:`M_m`): Minimum magnitude.
    * ``beta2`` (:math:`\beta_2`, [m]): Topographic gradient grid smoothing length-scale.
    * ``beta3`` (:math:`\beta_3`): Maximum throw-to-thickness ratio used when computing topographic gradient grid.
    * ``beta4`` (:math:`\beta_4`, [log\ :sub:`10`\ (MPa)]): Skeleton modulus exponent; :math:`H_s = 10^{\beta_4}`.
    * ``theta0`` (:math:`\theta_0`): Value for :math:`\beta`.
    * ``theta1`` (:math:`\theta_1`, [MPa\ :math:`^{-\theta_2}`\ ]): Weight factor for critical point scaling.
    * ``theta2`` (:math:`\theta_2`): Critical point scaling exponent.
    * ``theta3`` (:math:`\theta_3`, [MPa]): Stress critical point.

    The remaining parameters of the observed distribution are computed following:

    .. math::
       \beta_i &= \theta_0, \\
       \zeta_i &= \theta_1 (\theta_3 - \Delta C_i)^{\theta_2},

    where :math:`\Delta C_i` is the incremental Coulomb stress for the :math:`i`-th event.

    If :math:`\Delta C_i \geq \theta_3`, :math:`\zeta` will be set to zero.
    """
    _vars = StressDependentMagnitudeModel._vars + ('theta0', 'theta1', 'theta2', 'theta3')

    def __init__(self, thinsheet, label='ets.cps'):
        super().__init__(thinsheet, label)

    def _generate_parameter(self, name, var, covs):
        if name == 'Mmin':
            return var[name]
        if name == 'beta':
            t0 = var['theta0']
            return t0
        if name == 'zeta':
            t1, t2, t3 = var['theta1'], var['theta2'], var['theta3']
            estress = covs['stress_event']
            return theano.tensor.switch(estress < t3, t1 * (t3 - estress) ** t2, 0.)
        return super()._generate_parameter(name, var, covs)


class ExponentialTrendZeta(StressDependentMagnitudeModel):
    r"""Model with stress dependent :math:`\zeta`.

    The moment distribution is a power-law with an exponential taper. The variables controlling this model are:

    * ``Mmin`` (:math:`M_m`): Minimum magnitude.
    * ``beta2`` (:math:`\beta_2`, [m]): Topographic gradient grid smoothing length-scale.
    * ``beta3`` (:math:`\beta_3`): Maximum throw-to-thickness ratio used when computing topographic gradient grid.
    * ``beta4`` (:math:`\beta_4`, [log\ :sub:`10`\ (MPa)]): Skeleton modulus exponent; :math:`H_s = 10^{\beta_4}`.
    * ``theta0`` (:math:`\theta_0`): Value for :math:`\beta`.
    * ``theta1`` (:math:`\theta_1`): Base value for :math:`\zeta`.
    * ``theta2`` (:math:`\theta_2`, [MPa\ :sup:`-1`\ ]): Stress decay rate for :math:`\zeta`.

    The remaining parameters of the observed distribution are computed following:

    .. math::
       \beta_i &= \theta_0, \\
       \zeta_i &= \theta_1 e^{-\theta_2 \Delta C_i},

    where :math:`\Delta C_i` is the incremental Coulomb stress for the :math:`i`-th event.
    """
    _vars = StressDependentMagnitudeModel._vars + ('theta0', 'theta1', 'theta2')

    def __init__(self, thinsheet, label='ets.etz'):
        super().__init__(thinsheet, label)

    def _generate_parameter(self, name, var, covs):
        if name == 'Mmin':
            return var[name]
        if name == 'beta':
            t0 = var['theta0']
            return t0
        if name == 'zeta':
            t1, t2 = var['theta1'], var['theta2']
            estress = covs['stress_event']
            return t1 * np.exp(-t2 * estress)
        return super()._generate_parameter(name, var, covs)


class HyperbolicTangentBetaExponentialTrendZeta(StressDependentMagnitudeModel):
    r"""Model with stress dependent :math:`\beta` and :math:`\zeta`.

    The moment distribution is a power-law with an exponential taper. The variables controlling this model are:

    * ``Mmin`` (:math:`M_m`): Minimum magnitude.
    * ``beta2`` (:math:`\beta_2`, [m]): Topographic gradient grid smoothing length-scale.
    * ``beta3`` (:math:`\beta_3`): Maximum throw-to-thickness ratio used when computing topographic gradient grid.
    * ``beta4`` (:math:`\beta_4`, [log\ :sub:`10`\ (MPa)]): Skeleton modulus exponent; :math:`H_s = 10^{\beta_4}`.
    * ``theta0`` (:math:`\theta_0`): Base value for :math:`\beta`.
    * ``theta1`` (:math:`\theta_1`): Weight factor for hyperbolic tangent contribution to :math:`\beta`.
    * ``theta2`` (:math:`\theta_2`, [MPa\ :sup:`-1`\ ]): Stress variability rate for :math:`\beta`.
    * ``theta3`` (:math:`\theta_3`): Base value for :math:`\zeta`.
    * ``theta4`` (:math:`\theta_4`, [MPa\ :sup:`-1`\ ]): Stress decay rate for :math:`\zeta`

    The remaining parameters of the observed distribution are computed following:

    .. math::
       \beta_i &= \theta_0 + \theta_1 (1 - \tanh(\theta_2 \Delta C_i)), \\
       \zeta_i &= \theta_3 e^{-\theta_4 \Delta C_i},

    where :math:`\Delta C_i` is the incremental Coulomb stress for the :math:`i`-th event.
    """
    _vars = StressDependentMagnitudeModel._vars + ('theta0', 'theta1', 'theta2', 'theta3', 'theta4')

    def __init__(self, thinsheet, label='ets.htb.etz'):
        super().__init__(thinsheet, label)

    def _generate_parameter(self, name, var, covs):
        if name == 'Mmin':
            return var[name]
        estress = covs['stress_event']
        if name == 'beta':
            t0, t1, t2 = var['theta0'], var['theta1'], var['theta2']
            return t0 + t1 * (1 - np.tanh(t2 * estress))
        if name == 'zeta':
            t3, t4 = var['theta3'], var['theta4']
            return t3 * np.exp(-t4 * estress)
        return super()._generate_parameter(name, var, covs)
