# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
import re
import logging
from datetime import datetime

import numpy as np
import numpy.lib.recfunctions as rfn
import theano.tensor as tt
import theano.tensor.signal.conv as ttsignal
from dateutil.relativedelta import relativedelta
from scipy import signal, interpolate

from .. import utils

logger = logging.getLogger('geomodels')


class ReservoirGrid:
    """Manages 2D gridded reservoir data."""

    def __init__(self, xnodes, ynodes, names, attributes, fill=0):
        """Grid initialization can be performed by loading the grid data based on its attributes.

        Args:
            xnodes (array_like[float]): Easting RD coordinates of grid nodes (in meters). Shape ``(nx,)``
            ynodes (array_like[float]): Northing RD coordinates of grid nodes (in meters). Shape ``(ny,)``
            names (array_like[str]): Attribute names. Length ``nattrs``
            attributes (array_like[float], theano.tensor.TensorVariable): Grid attributes. Shape ``(nattrs, nx, ny)``
            fill (float): Value used for undefined grid data
        """
        self.__fill = fill

        self._names = tuple(names)
        self.__xnodes, self.__ynodes = np.array(xnodes), np.array(ynodes)

        eshape = len(self._names), self.__xnodes.size, self.__ynodes.size
        try:  # Input is a PyMC3Variable
            ashape = attributes.tag.test_value.shape
            self._attributes = attributes
        except AttributeError:
            try:  # Input is a Theano Tensor
                ashape = attributes.shape.eval()
                self._attributes = attributes
            except AttributeError:  # Input is array like
                ashape = np.shape(attributes)
                self._attributes = np.array(attributes)
        if ashape != eshape:
            raise ValueError('Cannot initialize {}: Mismatched dimensions'.format(self.__class__.__name__))

        xdiff, ydiff = np.diff(self.__xnodes), np.diff(self.__ynodes)
        self.__dx, self.__dy = xdiff[0], ydiff[0]
        if not np.allclose(xdiff, self.dx) or not np.allclose(ydiff, self.dx):
            raise ValueError('Cannot initialize {}: Grid spacing is not regular'.format(self.__class__.__name__))

        fmt_shape = '   {:5s}: {:9d} {:>2s}, {:5s}: {:9d} {:>2s}, {:5s}: {:9d} {:>2s}'
        fmt_float = '   {:5s}: {:9.1f} {:>2s}, {:5s}: {:9.1f} {:>2s}'
        logger.debug(fmt_shape.format('nx',   self.nx,   '',  'ny',   self.ny,   '', 'nattrs',   self.nattrs,   ''))
        logger.debug(fmt_float.format('dx',   self.dx,   'm', 'dy',   self.dy,   'm'))
        logger.debug(fmt_float.format('xmin', self.minx, 'm', 'xmax', self.maxx, 'm'))
        logger.debug(fmt_float.format('ymin', self.miny, 'm', 'ymax', self.maxy, 'm'))

    @property
    def minx(self):
        """float: Minimum Easting grid node coordinates (in meters)"""
        return np.min(self.__xnodes)

    @property
    def maxx(self):
        """float: Maximum Easting grid node coordinates (in meters)"""
        return np.max(self.__xnodes)

    @property
    def miny(self):
        """float: Minimum Northing grid node coordinates (in meters)"""
        return np.min(self.__ynodes)

    @property
    def maxy(self):
        """float: Maximum Northing grid node coordinates (in meters)"""
        return np.max(self.__ynodes)

    @property
    def dx(self):
        """float: Grid spacing in the Easting direction (in meters)"""
        return self.__dx

    @property
    def dy(self):
        """float: Grid spacing in the Northing direction (in meters)"""
        return self.__dy

    @property
    def nx(self):
        """int: Number of nodes in the Easting direction"""
        return self.__xnodes.size

    @property
    def ny(self):
        """int: Number of nodes in the Northing direction"""
        return self.__ynodes.size

    @property
    def nattrs(self):
        """int: Number of attributes present in grid"""
        return len(self._names)

    @property
    def names(self):
        """tuple[str]: Attribute names. Length ``nattrs``"""
        return self._names

    @property
    def xnodes(self):
        """np.ndarray: Easting nodes RD coordinates (in meters). Shape ``(nx,)``"""
        return self.__xnodes

    @property
    def ynodes(self):
        """np.ndarray: Northing nodes RD coordinates (in meters). Shape ``(ny,)``"""
        return self.__ynodes

    @property
    def attributes(self):
        """np.ndarray[float], theano.tensor.TensorVariable: Grid attributes. Shape ``(nattrs, nx, ny)``/``(nx, ny)``"""
        return self._attributes.squeeze()  # Remove unused dimension for single attribute case

    @property
    def mask(self):
        """np.ndarray[bool]: Mask of attributes that are considered undefined (see ``fill`` in initialization)"""
        return np.isnan(self._attributes) if np.isnan(self.__fill) else self._attributes == self.__fill

    def __getitem__(self, item):
        """Access grid values of an attribute.

        Args:
            item (str): Attribute name

        Returns:
            np.ndarray[float], pymc3.model.PyMC3Variable: Requested reservoir grid attribute. Shape ``(nx, ny)``
        """
        return self._attributes[self._names.index(item)]

    @classmethod
    def from_csv(cls, path, fill=0):
        """Import reservoir grid data coordinates and attributes from a CSV file.

        Assumes that the CSV file is organized as follows:

        =======  =======  ============  ============  ===  ============
              X        Y   attribute_1   attribute_2  ...   attribute_n
        =======  =======  ============  ============  ===  ============
        x coord  y coord         value         value  ...         value
            ...      ...           ...           ...  ...           ...
        x coord  y coord         value         value  ...         value
        =======  =======  ============  ============  ===  ============

        Args:
            path (str): Path of grid CSV file
            fill (float): Value used for undefined grid data

        Returns:
            seismod.geomodels.ReservoirGrid: Reservoir grid, as read from the file
        """
        logger.debug(' Reading reservoir grid data from: {}'.format(path))
        data = utils.import_csv(path, dtype='f8')
        xnodes, indx = np.unique(data['X'], return_inverse=True)
        ynodes, indy = np.unique(data['Y'], return_inverse=True)

        data = rfn.drop_fields(data, ('X', 'Y'))
        names = data.dtype.names

        attributes = np.full((len(names), xnodes.size, ynodes.size), fill, dtype='f8')
        attributes[:, indx, indy] = rfn.structured_to_unstructured(data).T

        return cls(xnodes, ynodes, names, attributes, fill)

    def snap(self):
        """Fills undefined nodes using spatial nearest neighbours."""
        if not np.any(self.mask):  # Nothing to snap
            return
        mask = self.mask
        nodes = np.vstack((np.repeat(self.__xnodes, self.ny), np.tile(self.__ynodes, self.nx))).T
        for i in range(self.nattrs):
            vn = nodes[~mask[i].ravel()]  # Valid (defined) nodes
            va = self._attributes[i, ~mask[i]]  # Valid (defined) attributes
            self._attributes[i] = interpolate.griddata(vn, va, nodes, method='nearest').reshape((self.nx, self.ny))

    def export_csv(self, path):
        """Write reservoir grid to a CSV file.

        Args:
            path (str): Path to output file
        """
        nodes = np.vstack((np.repeat(self.__xnodes, self.ny), np.tile(self.__ynodes, self.nx))).T
        output = np.hstack((nodes, self._attributes.reshape(self.nattrs, -1).T))
        utils.export_csv(path, output, header=','.join(('X', 'Y') + self._names))

    def _location_index(self, location):
        """Find grid node indices nearest to the given locations.

        Args:
            location (np.ndarray[float]): Requested RD coordinates (in meters). Shape ``(n, 2)``

        Returns:
            tuple[np.ndarray[int]]: Corresponding indices for Easting and Northing nodes. Both have shape ``(n,)``
        """
        i = np.clip(np.searchsorted(self.__xnodes - self.dx / 2., location[..., 0]) - 1, 0, self.nx - 1)
        j = np.clip(np.searchsorted(self.__ynodes - self.dy / 2., location[..., 1]) - 1, 0, self.ny - 1)
        return i, j

    def value_at_x(self, location):
        """Read attribute values of nodes nearest to the given locations.

        Args:
            location (np.ndarray[float]): RD coordinates (in meters). Shape ``(n, 2)``

        Returns:
            np.ndarray[float]: Attribute values at locations. Shape ``(nattrs, n)`` or ``(n,)`` if single attribute
        """
        i, j = self._location_index(location)
        return self._attributes[:, i, j].squeeze()

    def median_filter(self, size=3):
        """Apply a spatial median filter to all reservoir grid attributes.

        Undefined nodes will be ignored.

        Args:
            size (int): Size of kernel (number of nodes) used for filtering (must be positive and odd)
        """
        logger.debug('  Applying median filter with kernel size: {} px'.format(size))
        if size != int(size):
            raise ValueError('Kernel size for median filter must be an integer')
        if size < 0:
            raise ValueError('Kernel size for median filter must be positive')
        elif not size % 2:
            raise ValueError('Kernel size for median filter must be odd')
        dx = (size - 1) // 2
        mask = self.mask

        self._attributes[mask] = np.nan
        filtered = np.full((self.nattrs, self.nx, self.ny), self.__fill)
        for (k, i, j), val in np.ndenumerate(self._attributes):
            ir = np.clip(np.arange(i - dx, i + dx + 1), 0, self.nx - 1)
            jr = np.clip(np.arange(j - dx, j + dx + 1), 0, self.ny - 1)
            if not np.isnan(val):
                filtered[k, i, j] = np.nanmedian(self._attributes[k, ir, jr])
        self._attributes = filtered
        self._attributes[mask] = self.__fill

    def gaussian_filter(self, sigma):
        """Apply a standard Gaussian filter to all reservoir grid attributes.

        Args:
            sigma (float, pymc3.model.PyMC3Variable): Smoothing length-scale (in meters)
        """
        if isinstance(sigma, tt.TensorVariable):
            logger.debug('  Applying standard Gaussian filter with sigma: {}'.format(sigma.name))
        else:
            logger.debug('  Applying standard Gaussian filter with sigma: {} m'.format(sigma))
        nx, ny = self.nx if self.nx % 2 else self.nx + 1, self.ny if self.ny % 2 else self.ny + 1
        xcoord, ycoord = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
        kernel = np.exp(-((xcoord - nx // 2) ** 2 + (ycoord - ny // 2) ** 2) / (2 * (sigma / self.dx) ** 2))
        kernel /= kernel.sum()

        px, py = int(nx // 2), int(ny // 2)
        if isinstance(sigma, tt.TensorVariable) or isinstance(self._attributes, tt.TensorVariable):
            self._attributes = ttsignal.conv2d(self._attributes, kernel, border_mode='full')[:, px:-px, py:-py]
        else:
            padded = np.pad(self._attributes, ((0, 0), (px, px), (py, py)), mode='symmetric')
            self._attributes = signal.convolve(padded, kernel[None, ...], mode='full')[:, 2*px:-2*px, 2*py:-2*py]

    def adapted_gaussian_filter(self, sigma):
        r"""Apply a Gaussian filter to all reservoir grid attributes and update the grid attributes.

        Utilizes a custom method so as to be compatible with the geodesic compressibility data grid input file.
        It uses a standard convolution with constant, zero valued, padding for the convolution itself. It then adds
        the value of each original grid node, weighed by the one minus the area of overlap between the kernel and
        the individual grid node, e.g. for grid nodes around the center, this factor will be zero, while
        around corners it will be closer to 0.75.

        Args:
            sigma (float, pymc3.model.PyMC3Variable): Smoothing length-scale (in meters)
        """
        if isinstance(sigma, tt.TensorVariable):
            logger.debug('  Applying standard Gaussian filter with sigma: {}'.format(sigma.name))
        else:
            logger.debug('  Applying standard Gaussian filter with sigma: {} m'.format(sigma))
        nx, ny = self.nx if self.nx % 2 else self.nx + 1, self.ny if self.ny % 2 else self.ny + 1
        xcoord, ycoord = np.meshgrid(np.arange(nx), np.arange(ny), indexing='ij')
        kernel = np.exp(-((xcoord - nx // 2) ** 2 + (ycoord - ny // 2) ** 2) / (2 * (sigma / self.dx) ** 2))
        kernel /= kernel.sum()

        if isinstance(sigma, tt.TensorVariable) or isinstance(self._attributes, tt.TensorVariable):
            px, py = int(nx // 2), int(ny // 2)
            norm = ttsignal.conv2d(tt.ones_like(self._attributes), kernel, border_mode='full')[:, px:-px, py:-py]
            convolved = ttsignal.conv2d(self._attributes, kernel, border_mode='full')[:, px:-px, py:-py]
        else:
            norm = signal.convolve(np.ones_like(self._attributes), kernel[None, ...], mode='same')
            convolved = signal.convolve(self._attributes, kernel[None, ...], mode='same')
        self._attributes = convolved + self._attributes * (1 - norm)

    @staticmethod
    def match_nodes(*grids):
        """For all given grids, keep only nodes that are common to all grids.

        Args:
            grids (ReservoirGrid, None): Grid to be compared against. It will not be modified by this method
        """
        xnodes, ynodes, dx, dy = None, None, None, None
        for grid in grids:  # Select common nodes
            if grid is None:
                continue
            if xnodes is None and ynodes is None:
                xnodes, ynodes, dx, dy = grid.xnodes, grid.ynodes, grid.dx, grid.dy
            elif np.isclose(dx, grid.dx) and np.isclose(dy, grid.dy):
                xnodes, ynodes = np.intersect1d(xnodes, grid.xnodes), np.intersect1d(ynodes, grid.ynodes)
            else:
                raise ValueError('Cannot match nodes for grids with different grid spacing')

        if xnodes is None and ynodes is None:
            return
        logger.debug('Removing non-common nodes in grids')
        for grid in grids:  # Remove nodes not present in current grid
            if grid is None:
                continue
            xmask, = np.nonzero(np.in1d(grid.xnodes, xnodes))
            ymask, = np.nonzero(np.in1d(grid.ynodes, ynodes))
            if np.all(xmask) and np.all(ymask):
                continue  # Nothing to remove
            grid.__xnodes, grid.__ynodes = grid.xnodes[xmask], grid.ynodes[ymask]
            grid._attributes = grid._attributes[:, xmask][..., ymask]


class ReservoirGridEpochs(ReservoirGrid):
    """Manages 2D reservoir grid data with epochs as attributes.

    Supported formats for epochs attributes are: ``[P]YYYY_MM_DD`` or ``[P]YYYY``, where ``P`` is an optional,
    arbitrary prefix, which must be constant for all epochs. For the second accepted format, the day and month are
    assumed to be 1st of January.

    **Important:** Single epoch grids are not supported.
    """
    def __init__(self, xnodes, ynodes, names, attributes, fill=0):
        super().__init__(xnodes, ynodes, names, attributes, fill)

        # Set additional properties of grid with epochs
        self.__fmt = None
        self.__epochs = np.full(self.nattrs, 'NaT', dtype='M8[s]')
        for i, name in enumerate(self._names):
            self.__get_format(name)
            self.__epochs[i] = datetime.strptime(name, self.__fmt)

        initial, final = self.epochs[:-1], self.epochs[1:]
        dt, *deltas = [relativedelta(t2.astype(datetime), t1.astype(datetime)) for t1, t2 in zip(initial, final)]
        if deltas and any(delta != dt for delta in deltas):  # Check that epochs are equispaced
            raise ValueError('Cannot initialize {}: Non-constant epoch spacing'.format(self.__class__.__name__))
        if dt != relativedelta(years=1):
            raise NotImplementedError('{}: Only annual epochs are currently supported'.format(self.__class__.__name__))

        tmin, *_, tmax = self.epochs
        tmin, tmax = '{:%d/%m/%Y}'.format(tmin.astype(datetime)), '{:%d/%m/%Y}'.format(tmax.astype(datetime))
        logger.debug('   {:5s}: {:>12s}, {:5s}: {:>12s}'.format('tmin', tmin, 'tmax', tmax))

    @property
    def epochs(self):
        """np.ndarray[np.datetime64]: List of epochs contained in the reservoir grid"""
        return self.__epochs

    def __getitem__(self, item):
        """Access grid values of an attribute.

        Args:
            item (str, datetime.datetime): Attribute name or epoch

        Returns:
            np.ndarray: Requested reservoir grid attribute. Shape ``(nx, ny)``
        """
        try:  # Assumes item is a datetime object
            return super().__getitem__(self.__epoch2str(item))
        except AttributeError:
            try:  # Check whether item has np.datetime64 scalar (or other numpy scalar convertible to datetime)
                return super().__getitem__(self.__epoch2str(item.astype(datetime)))
            except AttributeError:  # Final try, defaulting to the ReservoirGrid (should be string)
                return super().__getitem__(item)

    def __epoch2str(self, epoch):
        """Returns given epoch as string formatted according to the grid attribute names.

        Args:
            epoch (datetime.datetime): Epoch to convert to string

        Returns:
            str: Grid formatted epoch string
        """
        return epoch.strftime(self.__fmt)

    def __get_format(self, name):
        """Sets the date formatting parameters based on given name.

        If formatting is already set, checks that the name conforms to existing format.

        Args:
            name (str): Name to convert to epoch
        """
        prefix, date = re.match(r'(\D*)(.+)', name).group(1, 2)
        try:
            datetime.strptime(date, '%Y_%m_%d'.format(date))
            fmt = '{}%Y_%m_%d'.format(prefix)
        except ValueError:
            datetime.strptime(date, '%Y'.format(date))
            fmt = '{}%Y'.format(prefix)
        if self.__fmt is None:
            self.__fmt = fmt
        elif self.__fmt != fmt:
            raise ValueError('{}: Multiple formats for epoch labels'.format(self.__class__.__name__))

    def period_attributes(self, start, end):
        """Grid attributes contained within period.

        Args:
            start (datetime.datetime): Start date for period (must be contained in grid epochs)
            end (datetime.datetime): End date for period (must be contained in grid epochs)

        Returns:
            np.ndarray[float]: Attributes within period. Shape ``(nepochs, nx, ny)``
        """
        return self._attributes[self._mask_period(start, end)].copy()

    def period_epochs(self, start, end):
        """Grid epochs contained within period.

        Args:
            start (datetime.datetime): Start date for period (must be contained in grid epochs)
            end (datetime.datetime): End date for period (must be contained in grid epochs)

        Returns:
            np.ndarray[np.datetime64]: Epochs within period. Shape ``(nepochs,)``
        """
        return self.epochs[self._mask_period(start, end)].copy()

    def _mask_period(self, start, end):
        """Generate mask of all epochs contained within period.

        Args:
            start (datetime.datetime): Start date for period (must be contained in grid epochs)
            end (datetime.datetime): End date for period (must be contained in grid epochs)

        Returns:
            np.ndarray[bool]: Mask of epochs within period. Shape ``(nepochs,)``
        """
        mask = (self.epochs >= start) & (self.epochs <= end)
        ti, *_, tf = self.epochs[mask]
        if ti != start or tf != end:
            raise ValueError('{}: Period start or end not present in grid data'.format(self.__class__.__name__))
        return mask

    def value_at_tx(self, time, location, nu=0):
        """Return grid values at given times and locations.

        Values are taken from the grid node closest to the location and interpolated in time with a cubic spline.

        It is also possible to compute partial derivatives in time, in which case the resulting rate is always
        on a per day basis.

        Args:
            location (np.ndarray[float]): RD coordinates (in meters). Shape ``(n, 2)``
            time (np.datetime64, np.ndarray[np.datetime64]): Times to interpolate into. Shape ``(n,)``
            nu (int): Temporal partial derivative order to compute

        Returns:
            float, np.ndarray[float]: Reservoir attributes at given times and locations. Shape ``(n,)``
        """
        *shape, ncoords = np.shape(location)
        if list(np.shape(time)) != shape:
            raise IndexError('Mismatch shape for location and time')
        grid_dt = (self.epochs - self.epochs[0]).astype('f8') / (24 * 3600)
        spline = interpolate.CubicSpline(grid_dt, self._attributes, extrapolate=False)
        try:
            i, j = self._location_index(np.reshape(location, (-1, ncoords)))
            eval_dt = (np.ravel(time) - self.epochs[0]).astype('f8') / (24 * 3600)
            return spline(eval_dt, nu=nu)[..., i, j].diagonal().reshape(shape)
        except MemoryError:
            if len(shape) <= 1:  # No obvious way of breaking events down
                raise
        # Too many total events for memory, try breaking them down using largest dimension
        out = np.empty(shape)
        laxis = shape.index(max(shape))
        chunk = tuple(size for dim, size in enumerate(shape) if dim != laxis)
        for k in range(max(shape)):
            # noinspection PyTypeChecker
            kk = (slice(None),) * laxis + (k,)
            i, j = self._location_index(np.reshape(location[kk], (-1, ncoords)))
            eval_dt = (np.ravel(time[kk]) - self.epochs[0]).astype('f8') / (24 * 3600)
            out[kk] = spline(eval_dt, nu=nu)[..., i, j].diagonal().reshape(chunk)
        return out
