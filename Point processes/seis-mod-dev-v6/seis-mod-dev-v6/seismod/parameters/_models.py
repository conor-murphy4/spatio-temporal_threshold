# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
from .. import seismodels
from .._bases import parameters, models


class ModelParameters(parameters.Parameters):
    """Class to read model specific parameters from a JSON file."""
    _def_params = dict(
        label='',
        type='',
        dof=dict()
    )

    def __init__(self, path, Mmin=1.5, **kwargs):
        """Accepted keys in JSON file are:

        * ``label``: Identifier string for the model (it is only used for identification)
        * ``type``: Model type (one of :mod:`seismod.seismodels`)
        * ``dof``: Definition of variables used in model (depends on ``type`` and must not include ``Mmin``). The
          variables themselves are defined using:

          * ``distr``: Variable distribution. One of ``constant``, ``uniform``, ``uniform-log`` or ``exponential``
          * ``value``: Variable value (required in ``constant``)
          * ``lower``: Lower bound of distribution (required in ``uniform``, ``uniform-log`` and ``exponential``)
          * ``upper``: Upper bound of distribution (required in ``uniform`` and ``uniform-log``)
          * ``mean``: Mean of the distribution (required in ``exponential``)

        Args:
            path (str): Path to JSON with input parameters
            Mmin (int): Minimum magnitude allowed for events in model
            **kwargs: Keyword arguments can be used to override parameters from input file
        """
        super().__init__(path, **kwargs)
        if 'Mmin' in self.dof:
            raise ValueError('Model parameters cannot explicitly set Mmin')
        self.dof['Mmin'] = {
            models.Options.PRIOR: models.Priors.CONST,
            models.Options.VALUE: Mmin
        }

    @property
    def is_magnitude(self):
        """bool: Whether model is a magnitude model"""
        return issubclass(getattr(seismodels, self._access('type')), models.StressInvariantMagnitudeModel)

    @property
    def is_invariant(self):
        """bool: Whether model is stress invariant"""
        if self.is_magnitude:
            return not issubclass(getattr(seismodels, self._access('type')), models.StressDependentMagnitudeModel)
        else:
            return not issubclass(getattr(seismodels, self._access('type')), models.StressDependentActivityRateModel)

    @property
    def dof(self):
        """dict[str, dict[str, float]]: Dictionary with variable definitions"""
        return self._access('dof')

    @property
    def values(self):
        """dict[str, float]: Dictionary with variable values (only functional if ``values`` is defined for all)"""
        return {name: params[models.Options.VALUE] for name, params in self.dof.items()}

    def model_generator(self, thinsheet):
        """Creates a model generator instance based on the parameters.

        Args:
            thinsheet (seismod.geomodels.ThinSheetModel): Thin-sheet model used for computing covariates

        Returns:
            seismod._bases.models._BaseModel: Model generator
        """
        return getattr(seismodels, self._access('type'))(thinsheet, label=self._access('label'))
