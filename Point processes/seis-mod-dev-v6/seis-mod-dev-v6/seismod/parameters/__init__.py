# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
"""This module contains utility classes for reading parameter for models/runs from JSON files.

Available classes:

* class::`~seismod.parameters.ModelParameters`: Read model specific parameters (independently of run parameters).
* class::`~seismod.parameters.SampleParameters`: Read sampling parameters (including thin-sheet parameters).
* class::`~seismod.parameters.CompareParameters`: Read comparison parameters (including thin-sheet parameters).
"""
from ._models import ModelParameters
from ._runners import SampleParameters, CompareParameters

__all__ = [
    'ModelParameters',
    'SampleParameters',
    'CompareParameters'
]
