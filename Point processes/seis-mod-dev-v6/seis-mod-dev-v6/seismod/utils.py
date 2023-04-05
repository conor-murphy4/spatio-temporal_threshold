# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
"""General, simple utilities."""
import os
import logging

import numpy as np
import numpy.lib.recfunctions as rfn

_c_default = 9.1
_d_default = 1.5


def mag2mom(mag, c=_c_default, d=_d_default):
    r"""Convert magnitudes to seismic moment using,

    .. math::
       M_0 = 10^{d M_w + c}

    Args:
        mag (float, np.ndarray[float]): Magnitudes to convert
        c (float): Parameter :math:`c` in equation
        d (float): Parameter :math:`d` in equation

    Returns:
        float, np.ndarray[float]: Equivalent seismic moment
    """
    return 10 ** (d * mag + c)


def mom2mag(mom, c=_c_default, d=_d_default):
    r"""Convert seismic moment to magnitude using,

    .. math::
       M_w = (\log_{10}M_0 - c) / d

    Args:
        mom (float, np.ndarray[float]): Seismic moment to convert
        c (float): Parameter :math:`c` in equation
        d (float): Parameter :math:`d` in equation

    Returns:
        float, np.ndarray[float]: Equivalent magnitude
    """
    return (np.log10(mom) - c) / d


def import_csv(path, **kwargs):
    """Reads a CSV file into a numpy array.

    By default, it reads the first row as field names for a structured array with the rest of the file contents as
    its contents. Also by default, the types for each field are inferred from the contents of the file.

    This function is a wrapper for ``numpy.genfromtxt`` with the options ``delimiter`` set to ``','``, ``names`` to
    ``True``, ``dtype`` to ``None`` and ``encoding`` to ``None``. With the exception of ``delimiter``, it is
    possible to override all the other options and/or provide additional options.

    Args:
        path (str): Path to CSV file to be read
        **kwargs: Additional keyword arguments. See documentation for ``numpy.genfromtxt``

    Returns:
        np.ndarray: By default, a structured array. See documentation for ``numpy.genfromtxt`` for more details
    """
    kwargs.setdefault('names', True)
    kwargs.setdefault('dtype', None)
    kwargs.setdefault('encoding', None)

    data = np.genfromtxt(path, delimiter=',', **kwargs)
    return data if len(data.shape) else np.array([data])


def export_csv(path, data, **kwargs):
    """Writes given data as CSV file.

    If the given path does not contain an extension, ``.csv`` will be added.

    If file directory does not exist, it will be created.

    This function is a wrapper for ``numpy.savetxt`` with the options ``delimiter`` set to ``','``, ``comments`` to
    ``''`` and ``fmt`` to ``'%s'``. With the exception of ``delimiter``, it is possible to override all the other
    options and/or provide additional options.

    Args:
        path (str): Path of CSV file to write
        data (np.ndarray): Data to be saved in file
        **kwargs: Additional keyword arguments. See documentation for ``numpy.savetxt``
    """
    kwargs.setdefault('comments', '')
    kwargs.setdefault('fmt', '%s')

    _, ext = os.path.splitext(path)
    if not ext:
        path += '.csv'
    if os.path.dirname(path):
        os.makedirs(os.path.dirname(path), exist_ok=True)
    if data.dtype.names is not None and 'header' not in kwargs:
        np.savetxt(path, rfn.structured_to_unstructured(data), delimiter=',', header=','.join(data.dtype.names), **kwargs)
    else:
        np.savetxt(path, data, delimiter=',', **kwargs)


def set_logger(level, loggers=('pymc3', 'geomodels', 'models', 'runners'), clear_handlers=True):
    """Sets up terminal logging for selected loggers.

    A single handler is used for all selected loggers. Optionally, it can remove all existing handlers beforehand.

    Args:
        level (int, string): Logging level in use (see documentation for ``logging``)
        loggers (list[str], tuple[str]): Sequence of logger names to set up
        clear_handlers (bool): Whether to remove all pre-existing handlers before setting the new one
    """
    fmt = logging.Formatter('%(asctime)-15s %(levelname)-7s %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    handler = logging.StreamHandler()
    handler.setFormatter(fmt)
    handler.setLevel(level)

    for name in loggers:
        logger = logging.getLogger(name)
        if clear_handlers:
            logger.handlers = []  # Needs to clean handlers to avoid double logging
        logger.addHandler(handler)
        logger.setLevel(level)
