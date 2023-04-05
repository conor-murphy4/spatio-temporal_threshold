# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
import abc

import pymc3 as pm
import theano.tensor


class BaseDistribution(pm.Continuous, abc.ABC):
    """Base class for distributions used within model.

    Expected parameter names must be provided as a tuple class variable in sub-classes.
    """
    parameters = ()  # MUST be filled in derived classes

    def __init__(self, **kwargs):
        """The object must be initialized with the expected parameters listed in
        :py:attr:`parameters`.

        Args:
            **kwargs: Distribution parameters and optional additional arguments (details for extra arguments `here
                <https://docs.pymc.io/api/distributions/utilities.html#pymc3.distributions.Distribution>`_)
        """
        for param in self.parameters:
            setattr(self, param, theano.tensor.as_tensor_variable(kwargs.pop(param)))
        super().__init__(**kwargs)

    def logp(self, value):
        """Calculate log-likelihood of the distribution at specified value(s).

        Args:
            value (float, np.ndarray[float]): Value(s) for which log-likelihood is calculated

        Returns:
            theano.tensor.TensorVariable: Total log-likelihood for the values
        """
        parameters = {name: getattr(self, name) for name in self.parameters}
        return self.loglike(value, **parameters)

    @classmethod
    @abc.abstractmethod
    def rvs(cls, size=None, rng=None, **kwargs):
        """Draw random values from the distribution.

        The parameters can be scalars or arrays, in which case the sampling will be performed for each parameter
        in the array (as long as all parameters are broadcastable with each other). This is set to facilitate sampling
        of full sets of parameters (such as traces).

        Args:
            size (int, tuple[int], None): Shape of random sample (returns one sample if not specified)
            rng (np.random.Generator, None): Random number generator used to generate the samples
            **kwargs: Distribution parameters

        Returns:
            float, np.ndarray[float]: Randomly sampled values
        """
        pass

    @classmethod
    @abc.abstractmethod
    def loglike(cls, value, **kwargs):
        """Calculate total log-likelihood of the distribution at specified value(s).

        The parameters can be scalars or arrays, in which case the sampling will be performed for each parameter
        in the array (as long as all parameters are broadcastable with each other and the values). This is set to
        facilitate computing likelihood for multiple catalogues simultaneously.

        Args:
            value (float, np.ndarray[float]): Value(s) for which log-likelihood is calculated
            **kwargs: Distribution parameter(s) at which log-likelihood will be calculated.

        Returns:
            float, np.ndarray[float]: Total log-likelihood for the input sample and parameters
        """
        pass
