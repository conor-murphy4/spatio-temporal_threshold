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
import datetime

import numpy as np

from .. import utils

logger = logging.getLogger('geomodels')


def _datetime2time(dates, hours):
    """Converts arrays of date and time strings into ``np.datetime64`` data used by catalogues.

    Args:
        dates (np.ndarray[str]): Date of event as string with format ``dd/mm/yyyy``
        hours (np.ndarray[str]): Time of event as string with format ``HH:MM:SS``

    Returns:
        np.ndarray[np.datetime64]: Dates and times of event (using second precision)
    """
    times = []
    for time in [' '.join(time_pair) for time_pair in zip(dates.astype(str), hours.astype(str))]:
        try:
            times.append(datetime.datetime.strptime(time, '%d/%m/%Y %H:%M:%S'))
        except ValueError:
            times.append('NaT')
    return np.array(times, dtype='<M8[s]')


def _time2datetime(time):
    """Converts arrays of ``np.datetime64`` data used by catalogues into date and time string data.

    Args:
        time (np.ndarray[np.datetime64]): Dates and times of event (using second precision)

    Returns:
        tuple[np.ndarray[str]]: Separated date and time string data
    """
    dates, hours = np.char.partition(np.char.replace(np.datetime_as_string(time), 'NaT', ''), 'T')[:, [0, 2]].T
    # np.datetime_as_string has format %Y-%m-%d. The next step changes it to NAM style date (%d/%m/%Y)
    dates = np.array(['/'.join(date[::-1]) for date in np.char.split(dates, '-')])
    return dates, hours


class EarthquakeCatalogue:
    """Manages a collection of earthquakes."""
    _HEADER = ['Num', 'Place', 'Date', 'Time_UT', 'Magnitude', 'Depth', 'Easting', 'Northing', 'Field', 'Province']
    _IMPORT_INFO = dict(
        num=dict(
            type='<U15',
            default='',
            column='Num'
        ),
        place=dict(
            type='<U35',
            default='',
            column='Place'
        ),
        time=dict(
            type='<M8[s]',
            default='NaT',
            column=['Date', 'Time_UT'],
            input_conversion=_datetime2time,
            output_conversion=_time2datetime
        ),
        magnitude=dict(
            type='f8',
            default=np.nan,
            column='Magnitude'
        ),
        location=dict(
            type=('f8', 3),
            default=np.nan,
            column=['Depth', 'Easting', 'Northing'],
            input_conversion=lambda depth, easting, northing: np.vstack((easting, northing, depth * 1e3)).T,
            output_conversion=lambda location: (location[:, 2] / 1e3, location[:, 0], location[:, 1])
        ),
        field=dict(
            type='<U15',
            default='',
            column='Field'
        ),
        province=dict(
            type='<U15',
            default='',
            column='Province'
        ),
        trial=dict(
            type='i4',
            default=0
        )
    )
    _DTYPE = [(name, form['type']) for name, form in _IMPORT_INFO.items()]

    def __init__(self, **kwargs):
        """Initialize earthquake catalogue with the given attributes attributes directly.

        Keyword Args:
            magnitude (np.ndarray[float]): Magnitude of the events (optional). Shape ``(nevents,)``
            location (np.ndarray[float]): RD coordinates of events (in meters, optional). Shape ``(nevents, 3)``
            time (np.ndarray[np.datetime64]): Event dates (optional). Shape ``(nevents,)``
            place (np.ndarray[str]): Name of the event location (optional). Shape ``(nevents,)``
            field (np.ndarray[str]): Name of the gas-field associated (optional). Shape ``(nevents,)``
        """
        shape = 0
        for keys, values in kwargs.items():
            cshape = np.shape(values)
            if keys == 'location':  # Deal with size 3 dim for East, North, Depth
                cshape = cshape[:-1]
            if shape == 0:
                shape = cshape
            elif shape != cshape:
                raise ValueError('Cannot initialize {}: Mismatched dimensions'.format(self.__class__.__name__))
        if shape == 0:  # Empty catalogue
            self.__catalogue = np.zeros(shape, dtype=self._DTYPE)
            return

        if len(shape) > 2:
            raise ValueError('Cannot initialize {}: Too many dimensions'.format(self.__class__.__name__))
        self.__catalogue = np.zeros(shape, dtype=self._DTYPE)
        for name, form in self._IMPORT_INFO.items():
            self.__catalogue[name] = kwargs.pop(name, form['default'])

        if len(shape) == 2 and np.all(self.trial == self._IMPORT_INFO['trial']['default']):
            self.trial[...] = np.arange(0, shape[1])
        elif len(shape) == 2:
            raise ValueError('Cannot initialize {}: 2D-catalogue with trial'.format(self.__class__.__name__))

        if kwargs:  # Check that there are no extra keywords
            raise KeyError('Cannot initialize {}: Unrecognized arguments'.format(self.__class__.__name__))

    @property
    def magnitude(self):
        """np.ndarray[float]: Event magnitudes. Shape ``(nevents,)``"""
        return self.__catalogue['magnitude']

    @property
    def location(self):
        """np.ndarray[float]: RD coordinates and depth of each event (in meters). Shape ``(nevents, 3)``"""
        return self.__catalogue['location']

    @property
    def time(self):
        """np.ndarray[np.datetime64]: Date and time of each event. Shape ``(nevents,)``"""
        return self.__catalogue['time']

    @property
    def place(self):
        """np.ndarray[str]: Name of locality where each event occurred. Shape ``(nevents,)``"""
        return self.__catalogue['place']

    @property
    def field(self):
        """np.ndarray[str]: Name of gas-field associated with each event. Shape ``(nevents,)``"""
        return self.__catalogue['field']

    @property
    def trial(self):
        """np.ndarray[int]: Catalogue index associated with each event. Shape ``(nevents,)``"""
        return self.__catalogue['trial']

    @property
    def dr(self):
        """np.ndarray[float]: Inter-event distance (in meters), as a matrix. Shape ``(nevents, nevents)``"""
        return np.linalg.norm(self.location[:, None] - self.location, axis=-1)

    @property
    def dt(self):
        """np.ndarray[float]: Inter-event times (in days), as a matrix. Shape ``(nevents, nevents)``"""
        return (self.time[:, None] - self.time).astype('f8') / (24 * 3600)

    @property
    def nevents(self):
        """tuple[int]: Shape of events in the catalogue"""
        return self.__catalogue.shape

    @property
    def Mmax(self):
        """float: Maximum magnitude of events in catalogue, ``NaN`` if there are no events"""
        if np.prod(self.nevents) == 0:
            return np.nan
        if len(self.nevents) > 1:
            return np.max(self.magnitude, axis=0)
        return np.max(self.magnitude)

    @classmethod
    def from_csv(cls, path):
        """Imports earthquake catalogue from CSV file.

        Assumes that the CSV file is organized as follows:

        ===  =====  ==========  ========  =========  ==========  =========  =========  =====  ========
        Num  Place        Date   Time_UT  Magnitude       Depth    Easting   Northing  Field  Province
        ===  =====  ==========  ========  =========  ==========  =========  =========  =====  ========
        tag   name  DD/MM/YYYY  HH:MM:SS      value  depth [km]  coord [m]  coord [m]   name      name
        ...    ...         ...       ...        ...         ...        ...        ...    ...       ...
        tag   name  DD/MM/YYYY  HH:MM:SS      value  depth [km]  coord [m]  coord [m]   name      name
        ===  =====  ==========  ========  =========  ==========  =========  =========  =====  ========

        Args:
            path (str): Path to earthquake catalogue

        Returns:
            seismod.geomodels.EarthquakeCatalogue: Structured array with catalogue data formatted
        """
        logger.debug(' Reading earthquake catalogue from: {}'.format(path))

        data = utils.import_csv(path, filling_values='')
        attributes = dict()
        for name, form in cls._IMPORT_INFO.items():
            if 'input_conversion' in form:
                attributes[name] = form['input_conversion'](*[data[col] for col in form['column']])
            elif 'column' in form:
                attributes[name] = data[form['column']]
        return cls(**attributes)

    def export_csv(self, path):
        """Writes catalogue to CSV file using NAM formatting.

        Args:
            path (str): Path to CSV file where data will be saved
        """
        outdata = np.zeros((self.__catalogue.size, len(self._HEADER)), dtype='O')
        for name, form in self._IMPORT_INFO.items():
            if 'output_conversion' in form:
                for sname, svalue in zip(form['column'], form['output_conversion'](self.__catalogue[name])):
                    outdata[:, self._HEADER.index(sname)] = svalue
            elif 'column' in form:
                outdata[:, self._HEADER.index(form['column'])] = self.__catalogue[name]
        utils.export_csv(path, outdata, header=','.join(self._HEADER))

    def sort(self, order='time'):
        """Sorts catalogues.

        Args:
            order (str): Key to sort catalogues by
        """
        if len(self.nevents) > 1:
            self.__catalogue.sort(order=order, axis=0)
        self.__catalogue.sort(order=order)

    def _selection_mask(self, **kwargs):
        """Generates mask of catalogue which fit the selection criteria (see ``select``)."""
        mask = np.ones(self.__catalogue.shape, dtype='?')
        for key, value in kwargs.items():
            criteria, minmax = key[0], key[1:]
            if criteria == 'M':
                selector = self.magnitude
            elif criteria == 't':
                selector = self.time
            elif criteria == 'x':
                selector = self.location[:, 0]
            elif criteria == 'y':
                selector = self.location[:, 1]
            elif criteria == 'z':
                selector = self.location[:, 2]
            else:
                raise KeyError('Unrecognized keyword argument: {}'.format(key))

            if minmax == 'min':
                mask &= selector >= value
            elif minmax == 'max':
                mask &= selector < value
            else:
                raise KeyError('Unrecognized keyword argument: {}'.format(key))

        return mask

    def select(self, **kwargs):
        """Select events from catalogue that fit all the given parameters.

        Keyword Args:
            Mmin (float): Minimum magnitude
            Mmax (float): Maximum magnitude
            tmin (datetime.datetime): Minimum date and time
            tmax (datetime.datetime): Maximum date and time
            xmin (float): Minimum Easting (in RD coordinates in meters)
            xmax (float): Maximum Easting (in RD coordinates in meters)
            ymin (float): Minimum Northing (in RD coordinates in meters)
            ymax (float): Maximum Northing (in RD coordinates in meters)
            zmin (float): Minimum depth (in meters)
            zmax (float): Maximum depth (in meters)
            trial (np.array[int]):
        """
        if 'trial' in kwargs:
            trial = kwargs.pop('trial')
            if np.any(np.isin(trial, self.trial, invert=True)):  # Check that there are no non-existent trial
                logger.warning('{}: Skipping selection of non-existing trials'.format(self.__class__.__name__))
            selection = np.isin(self.trial, trial)
            nevents = self.nevents
            self.__catalogue = self.__catalogue[selection]
            if len(nevents) > 1:
                self.__catalogue = self.__catalogue.reshape((nevents[0], -1))

        mask = self._selection_mask(**kwargs)
        if np.all(mask):
            return
        if len(self.nevents) > 1:
            logger.warning('{}: Masked selection in 2D catalogue flattens it'.format(self.__class__.__name__))
        self.__catalogue = self.__catalogue[mask]
