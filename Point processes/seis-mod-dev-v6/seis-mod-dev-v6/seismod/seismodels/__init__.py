# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
"""Contains classes that handle different magnitude and activity rate models.

Each class is capable of parameter sampling, evaluation of out of sample log-likelihood and simulate either
time and location of events (for activity rate models) or event magnitudes (for magnitude models) for a given set
of parameters.

Moment-magnitude models present:

* :class:`~seismod.seismodels.UniformBetaZeta`
* :class:`~seismod.seismodels.UniformBetaCornerMoment`
* :class:`~seismod.seismodels.UniformBetaCornerMagnitude`
* :class:`~seismod.seismodels.InversePowerLawBeta`
* :class:`~seismod.seismodels.HyperbolicTangentBeta`
* :class:`~seismod.seismodels.CriticalPointScalingZeta`
* :class:`~seismod.seismodels.ExponentialTrendZeta`
* :class:`~seismod.seismodels.HyperbolicTangentBetaExponentialTrendZeta`

Activity rate models present:

* :class:`~seismod.seismodels.UniformRate`
* :class:`~seismod.seismodels.ExponentialTrendCoulomb`
* :class:`~seismod.seismodels.ExponentialTrendCoulombEpidemicTypeAftershockSequence`
"""
from ._magnitude import UniformBetaZeta, UniformBetaCornerMoment, UniformBetaCornerMagnitude
from ._magnitude import InversePowerLawBeta, HyperbolicTangentBeta, CriticalPointScalingZeta, ExponentialTrendZeta
from ._magnitude import HyperbolicTangentBetaExponentialTrendZeta
from ._rate import UniformRate, ExponentialTrendCoulomb, ExponentialTrendCoulombEpidemicTypeAftershockSequence

__all__ = [
    'UniformBetaZeta',
    'UniformBetaCornerMoment',
    'UniformBetaCornerMagnitude',
    'InversePowerLawBeta',
    'HyperbolicTangentBeta',
    'CriticalPointScalingZeta',
    'ExponentialTrendZeta',
    'HyperbolicTangentBetaExponentialTrendZeta',
    'UniformRate',
    'ExponentialTrendCoulomb',
    'ExponentialTrendCoulombEpidemicTypeAftershockSequence'
]
