# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
"""This module contains utility classes for handling and storing geophysical information needed in seismicity models.

The classes present are:

* :class:`~seismod.geomodels.EarthquakeCatalogue`: Reads, stores and performs event selections of catalogues.
* :class:`~seismod.geomodels.ReservoirGrid`: Reads and performs basic processing on grids of reservoir attributes.
* :class:`~seismod.geomodels.ReservoirGridEpochs`: As above, but for attributes whose values are time dependent.
* :class:`~seismod.geomodels.FaultModel`: Reads and processes fault data for generating topographic gradient grids.
* :class:`~seismod.geomodels.ThinSheetModel`: Reads and processes grid and fault data for full reservoir model.
"""
from ._reservoirs import ReservoirGrid, ReservoirGridEpochs
from ._earthquakes import EarthquakeCatalogue
from ._faults import FaultModel
from ._thinsheet import ThinSheetModel


__all__ = [
    'EarthquakeCatalogue',
    'ReservoirGrid',
    'ReservoirGridEpochs',
    'FaultModel',
    'ThinSheetModel',
]
