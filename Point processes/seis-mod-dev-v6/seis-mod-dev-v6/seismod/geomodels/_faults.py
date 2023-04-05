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
import numpy.lib.recfunctions as rfn

from ._reservoirs import ReservoirGrid
from .. import utils

logger = logging.getLogger('geomodels')


class FaultModel:
    def __init__(self, locations, throws, thickness, xnodes, ynodes, dmin=85.):
        """Initialize geological fault model from CSV file.

        Args:
            locations (array_like[float]): Fault node locations in RD coordinates (in meters). Shape ``(nfaults, 2)``
            throws (array_like[float]): Fault node throw values (in meters). Shape ``(nfaults,)``
            thickness (array_like[float]): Average fault node thickness (in meters). Shape ``(nfaults,)``
            xnodes (array_like[float]): Easting RD coordinates of grid nodes (in meters). Shape ``(nx,)``
            ynodes (array_like[float]): Northing RD coordinates of grid nodes (in meters). Shape ``(ny,)``
            dmin (float): Only fault nodes separated by more than this distance (in meters) are included in the model
        """
        self._xnodes, self._ynodes = np.array(xnodes), np.array(ynodes)
        nx, ny = self._xnodes.size, self._ynodes.size
        self._dx, *_ = np.diff(self._xnodes)
        self._dy, *_ = np.diff(self._ynodes)
        if not np.allclose(np.diff(self._xnodes), self._dx) or not np.allclose(np.diff(self._ynodes), self._dx):
            raise ValueError('Cannot initialize {}: Grid spacing is not regular'.format(self.__class__.__name__))

        locations, throws, thickness = np.array(locations), np.array(throws), np.array(thickness)
        keep, prev_node = np.ones(thickness.size, dtype='?'), None
        for i, (node, throw, thick) in enumerate(zip(locations, throws, thickness)):
            if i == 0:
                prev_node = node
                continue
            elif np.sqrt(np.sum((node - prev_node) ** 2)) <= dmin:
                keep[i] = False
            else:
                prev_node = node
        locations, throws, thickness = locations[keep], throws[keep], thickness[keep]

        ii = np.clip(np.searchsorted(self._xnodes - self._dx / 2., locations[..., 0]) - 1, 0, nx - 1)
        jj = np.clip(np.searchsorted(self._ynodes - self._dy / 2., locations[..., 1]) - 1, 0, ny - 1)
        unique, counts = np.unique(np.array([ii, jj]).T, axis=0, return_counts=True)
        ncounts = np.max(counts)
        self.__ratio, self.__thick = np.zeros((ncounts, nx, ny)), np.zeros((ncounts, nx, ny))
        for (i, j), count in zip(unique, counts):
            mask = ii + jj * ny == i + j * ny
            self.__ratio[:count, i, j] = throws[mask] / thickness[mask]
            self.__thick[:count, i, j] = thickness[mask]

    @property
    def ratio(self):
        """np.ndarray[float]: Throw-to-thickness ratio values in the grid. Shape ``(ncounts, nx, ny)``"""
        return self.__ratio

    @property
    def thickness(self):
        """np.ndarray[float]: Fault thickness values in the grid (in meters). Shape ``(ncounts, nx, ny)``"""
        return self.__thick

    @classmethod
    def from_csv(cls, path, xnodes, ynodes, dmin=85.):
        """Imports fault model from a CSV file.

        Assumes that the CSV file contains at least the following columns:

        =========  =========  =========  =========  ==============  ==============  ===============  ===============
          Easting   Northing   ZTopLeft  ZTopRight  ThicknessLeft1  ThicknessLeft2  ThicknessRight1  ThicknessRight2
        =========  =========  =========  =========  ==============  ==============  ===============  ===============
        coord [m]  coord [m]  depth [m]  depth [m]   thickness [m]   thickness [m]    thickness [m]    thickness [m]
              ...        ...        ...        ...             ...             ...              ...              ...
        coord [m]  coord [m]  depth [m]  depth [m]   thickness [m]   thickness [m]    thickness [m]    thickness [m]
        =========  =========  =========  =========  ==============  ==============  ===============  ===============

        Additional columns and the order of the columns is ignored. All units must be meters.

        Args:
            path (str): Path to fault model CSV file
            xnodes (array_like[float]): Easting RD coordinates of grid nodes (in meters). Shape ``(nx,)``
            ynodes (array_like[float]): Northing RD coordinates of grid nodes (in meters). Shape ``(ny,)``
            dmin (float): Only fault nodes separated by more than this distance (in meters) are included in the model

        Returns:
            seismod.geomodels.FaultModel: Generated fault model
        """
        logger.debug(' Reading fault data from: {}'.format(path))
        data = utils.import_csv(path)

        locations = rfn.structured_to_unstructured(data[['Easting', 'Northing']])
        throws = np.abs(data['ZTopRight'] - data['ZTopLeft'])
        thickness_columns = ['ThicknessLeft1', 'ThicknessLeft2', 'ThicknessRight1', 'ThicknessRight2']
        thickness = np.mean(rfn.structured_to_unstructured(data[thickness_columns]), axis=-1)

        return cls(locations, throws, thickness, xnodes, ynodes, dmin)

    def generate_topographic_gradients(self, rmax=np.inf, sigma=0):
        r"""Creates a topographic gradient grid based on the fault model.

        The fault nodes with the highest throw (up to a throw-to-thickness ratio of :math:`r_\mathrm{max}`)
        inside each grid node are selected to create an "unsmoothed" gradient grid. This is then passed through
        a Gaussian filter with length-scale given by :math:`\sigma`.

        Args:
            rmax (float, theano.tensor.TensorVariable): Maximum throw-to-thickness ratio
            sigma (float, theano.tensor.TensorVariable): Smoothing length-scale (in meters)

        Returns:
            seismod.geomodels.ReservoirGrid: Topographic gradients grid. Attribute name is ``topo_gradients``
        """
        logger.info(' Generating topographic gradient grid')
        logger.debug('  Filtering using maximum throw-to-thickness ratio: {}'.format(rmax))

        throw = (self.ratio * self.thickness * (self.ratio < rmax)).max(axis=0, keepdims=True) / (self._dx + self._dy)
        grid = ReservoirGrid(xnodes=self._xnodes, ynodes=self._ynodes, attributes=throw, names=('topo_gradients',))
        if sigma:
            grid.gaussian_filter(sigma)
        return grid
