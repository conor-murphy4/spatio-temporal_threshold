# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
import logging

import numpy as np
from methodtools import lru_cache

from ._faults import FaultModel
from ._reservoirs import ReservoirGrid, ReservoirGridEpochs

logger = logging.getLogger('geomodels')


class ThinSheetModel:
    """Class for handling thin-sheet models and computing covariates."""

    def __init__(self, paths, kernels=None, sigmas=None):
        r"""Initialize and pre-process reservoir grids and fault data used in the model.

        Note that all grids and fault data are, in principle, optional. However, lacking one or more of them will
        make it impossible to extract certain covariates. This can be the intended behavior, depending on the model
        that uses the thin-sheet model.

        Reservoir grids can, optionally, be median filtered and/or smoothed with a Gaussian. For the latter, undefined
        grid values will be filled in using :meth:`~seismod.geomodels.ReservoirGrid.snap` before smoothing. Finally,
        the pressure grid is smoothed using :meth:`~seismod.geomodels.ReservoirGrid.adapted_gaussian_filter`, while
        all other grids use :meth:`~seismod.geomodels.ReservoirGrid.gaussian_filter`.

        Recognized inputs in ``paths`` are:

        * ``pressure``: Reservoir pressure data file for all epochs [bar].
        * ``compressibility``: Reservoir compressibility data file [MPa\ :sup:`-1`\ ].
        * ``thickness``: Reservoir thickness data file [m].
        * ``faults``: Fault data file (see :class:`~seismod.geomodels.FaultModel`). Does not support median filtering
          nor Gaussian smoothing.

        Args:
            paths (dict[str, str]): Dictionary with paths to grids/fault data
            kernels (dict[str, int, None], None): Kernel size for median filtering reservoir grids
            sigmas (dict[str, float, None], None): Smoothing length-scale for Gaussian filtering of reservoir grids
        """
        logger.info('Generating thinsheet model')

        sigmas = sigmas or dict()
        kernels = kernels or dict()

        self.__depletion = self.__generate_depletion_from_pressure(paths, kernels, sigmas)
        self.__thickness = self.__generate_standard_grid('thickness', paths, kernels, sigmas, fill=1)
        self.__compressibility = self.__generate_standard_grid('compressibility', paths, kernels, sigmas)

        ReservoirGrid.match_nodes(self.__depletion, self.__thickness, self.__compressibility)

        self.__fault_model = None
        if paths.get('faults'):
            self.__fault_model = FaultModel.from_csv(paths['faults'], self.__depletion.xnodes, self.__depletion.ynodes)

    @property
    def ds(self):
        r"""float: Grid node area [m\ :sup:`2`\ ]"""
        return self.__depletion.dx * self.__depletion.dy

    @staticmethod
    def __generate_standard_grid(name, paths, kernels, sigmas, fill=0):
        """Generates a reservoir grid based on inputs.

        Args:
            name (str): Name of grid to generate
            paths (dict[str, str]): Dictionary with paths to grids/fault data
            kernels (dict[str, int, None], None): Kernel size for median filtering reservoir grids
            sigmas (dict[str, float, None], None): Smoothing length-scale for Gaussian filtering of reservoir grids
            fill (float): Number used to fill in missing grid nodes

        Returns:
            seismod.geomodels.ReservoirGrid: Reservoir grid
        """
        if not paths.get(name):
            return None
        standard_grid = ReservoirGrid.from_csv(paths.get(name), fill)
        if kernels.get(name, None):
            standard_grid.median_filter(kernels.get(name))
        if sigmas.get(name, None):
            standard_grid.snap()
            standard_grid.gaussian_filter(sigmas.get(name))
        return standard_grid

    @staticmethod
    def __generate_depletion_from_pressure(paths, kernels, sigmas, fill=0):
        """Create the depletion grid with respect to the first epoch in pressure.

        Pressure recovery (negative depletion) will be set to zero.

        Args:
            paths (dict[str, str]): Dictionary with paths to grids/fault data
            kernels (dict[str, int, None], None): Kernel size for median filtering reservoir grids
            sigmas (dict[str, float, None], None): Smoothing length-scale for Gaussian filtering of reservoir grids
            fill (float): Number used to fill in missing grid nodes in pressure grid

        Returns:
            seismod.geomodels.ReservoirGridEpochs: Depletion grid (in MPa)
        """
        if not paths.get('pressure'):
            return None
        pressure = ReservoirGridEpochs.from_csv(paths.get('pressure'), fill)
        if kernels.get('pressure', None):
            pressure.median_filter(kernels.get('pressure'))
        if sigmas.get('pressure', None):
            pressure.snap()
            pressure.adapted_gaussian_filter(sigmas.get('pressure'))

        logger.debug(' Generating depletion grid')
        depletion = (pressure.attributes[0] - pressure.attributes) / 10.  # Convert units from bar to MPa
        depletion[depletion < 0] = 0  # Model does not handle pressure recovery, set to zero in those cases
        return ReservoirGridEpochs(pressure.xnodes, pressure.ynodes, pressure.names, depletion, fill=0)

    @staticmethod
    def _compute_stress(depletion, compressibility, topogradients, beta4):
        r"""Calculates incremental Coulomb stress.

        Args:
            depletion (np.ndarray[float], theano.tensor.TensorVariable): Depletion data (in MPa)
            compressibility (np.ndarray[float], theano.tensor.TensorVariable): Compressibility (in MPa\ :sup:`-1`\ )
            topogradients (np.ndarray[float], theano.tensor.TensorVariable): Topographic gradient data
            beta4 (float, theano.tensor.TensorVariable): Skeleton module parameter; :math:`H_s = 10^{\beta_4}`

        Returns:
            np.ndarray[float], theano.tensor.TensorVariable: Incremental Coulomb stress
        """
        hsinv = 10. ** -beta4  # reciprocal of Hs [MPa^-1]
        hrinv = compressibility
        h = 1. / (hrinv + hsinv)
        return h * topogradients * depletion * compressibility

    @lru_cache(maxsize=4)
    def _topographic_gradients(self, beta2, beta3):
        """Generate a topographic gradient grid, based on fault model data and parameters

        Args:
            beta2 (float, theano.tensor.TensorVariable): Smoothing length-scale (in meters)
            beta3 (float, theano.tensor.TensorVariable): Maximum throw-to-thickness ratio

        Returns:
            seismod.geomodels.ReservoirGrid: Topographic gradients grid
        """
        return self.__fault_model.generate_topographic_gradients(rmax=beta3, sigma=beta2)

    def stress_grid(self, start, end, beta2, beta3, beta4, full=True):
        """Calculates incremental Coulomb stress grid within given period.

        Args:
            start (datetime.datetime): Period start date
            end (datetime.datetime): Period end date
            beta2 (float, theano.tensor.TensorVariable): Smoothing length-scale (in meters)
            beta3 (float, theano.tensor.TensorVariable): Maximum throw-to-thickness ratio
            beta4 (float, theano.tensor.TensorVariable): Skeleton module parameter; :math:`H_s = 10^{\beta_4}`
            full (bool): Return a grid for every epochs, otherwise return only first and last epoch

        Returns:
            np.ndarray[float], theano.tensor.TensorVariable: Stress grid [MPa]. Shape ``(nepochs, nx, ny)``
        """
        topogradients = self._topographic_gradients(beta2, beta3).attributes
        compressibility = self.__compressibility.attributes
        depletion = self.__depletion.period_attributes(start, end)
        if full:
            return self._compute_stress(depletion, compressibility, topogradients, beta4)
        return self._compute_stress(depletion[[0, -1]], compressibility, topogradients, beta4)

    def uniform_grid(self, start, end, full=True):
        """Calculates a (spatially) uniform grid within given period.

        The grid values along time linearly increase.

        Args:
            start (datetime.datetime): Period start date
            end (datetime.datetime): Period end date
            full (bool): Return a grid for every epochs, otherwise return only first and last epoch

        Returns:
            np.ndarray[float], theano.tensor.TensorVariable: Uniform grid [day]. Shape ``(nepochs, nx, ny)``
        """
        # Linear increase is required to have constant event density
        epochs = self.__depletion.period_epochs(start, end)
        increment = (epochs - epochs[0]).astype('f8') / (24 * 3600)
        if full:
            return (self.__depletion.period_attributes(start, end) > 0) * increment[:, None, None]
        return (self.__depletion.period_attributes(start, end)[[0, -1]] > 0) * increment[[0, -1], None, None]

    def thickness_grid(self):
        """Compute the full thickness grid.

        Returns:
            np.ndarray[float]: Thickness grid (in meters). Shape ``(nx, ny)``
        """
        return self.__thickness.attributes

    def stress_event(self, time, location, beta2, beta3, beta4, nu=0):
        r"""Computes incremental Coulomb stress for each event.

        It is possible to compute the incremental Coulomb stress partial derivative in time as well.

        Args:
            time (np.ndarray[np.datetime64]): Event times. Shape ``(nevents,)``
            location (np.ndarray[np.float]): Event locations in RD coordinates (in meters). Shape ``(nevents, 3)``
            beta2 (float, theano.tensor.TensorVariable): Smoothing length-scale (in meters)
            beta3 (float, theano.tensor.TensorVariable): Maximum throw-to-thickness ratio
            beta4 (float, theano.tensor.TensorVariable): Skeleton module parameter; :math:`H_s = 10^{\beta_4}`
            nu (int): Temporal partial derivative order to compute

        Returns:
            np.ndarray[float], theano.tensor.TensorVariable: Event Coulomb stress (or rat). If :math:`\nu = 0` then
            units are in [MPa], while for :math:`\nu = 1` units are [MPa day\ :sup:`-1`\ ] and so forth. Shape
            ``(nevents,)``
        """
        topogradients = self._topographic_gradients(beta2, beta3).value_at_x(location)
        compressibility = self.__compressibility.value_at_x(location)
        depletion = self.__depletion.value_at_tx(time, location, nu=nu)
        return self._compute_stress(depletion, compressibility, topogradients, beta4)

    def thickness_event(self, location):
        """Thickness value for each event.

        Args:
            location (np.ndarray[np.float]): Event locations in RD coordinates (in meters). Shape ``(nevents, 3)``

        Returns:
            np.ndarray[float]: Event thickness (in meters). Shape ``(nevents,)``
        """
        return self.__thickness.value_at_x(location)

    def relcoords2abscoords(self, start, end, itime, ix, iy, z):
        """Converts relative coordinates (in terms of grid indices) into absolute coordinates (times and locations).

        Note that the indices are floats. The integer part determines the corresponding index in the grid, while
        the decimal part is the an offset for the precise location: The factor of the corresponding grid spacing
        needed to get the exact location.

        Args:
            start (datetime.datetime): Period start date
            end (datetime.datetime): Period end date
            itime (np.ndarray[float]): Event epoch indices with respect to period. Shape ``(nevents,)``
            ix (np.ndarray[float]): Event x-nodes indices. Shape ``(nevents,)``
            iy (np.ndarray[float]): Event y-nodes indices. Shape ``(nevents,)``
            z (float, np.ndarray[float]): Event depth (in meters). Shape ``(nevents,)``

        Returns:
            tuple[np.ndarray]: Times (`np.datetime64`, shape ``(nevents,)``) and locations (`float`,
            shape ``(nevents, 3)``) of the generated events. Locations are in meters in RD coordinates
        """
        epochs = self.__depletion.period_epochs(start, end)
        iepoch = np.floor(itime).astype('i8')
        time = epochs[iepoch] + np.diff(epochs)[iepoch] * (itime - iepoch)

        if np.shape(ix) != np.shape(iy):
            raise TypeError('Shape for x and y locations do not match')
        location = np.empty(np.shape(ix) + (3,))
        location[..., 0] = self.__depletion.minx + ix * self.__depletion.dx
        location[..., 1] = self.__depletion.miny + iy * self.__depletion.dy
        location[..., 2] = z

        return time, location
