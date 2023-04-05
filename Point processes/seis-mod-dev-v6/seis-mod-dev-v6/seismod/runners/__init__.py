# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
"""This module contains utility classes for handling common runs involving seismological models.

The classes present are:

* :class:`~seismod.runners.ModelSample`: Class for running standard variable sampling of models.
* :class:`~seismod.runners.ModelCompare`: Class for running comparison of model performance over synthetic catalogues.
"""
from ._sample import ModelSample
from ._compare import ModelCompare


__all__ = [
    'ModelSample',
    'ModelCompare'
]


def _main_compare():
    """Runs a comparison of the out-of-sample log-likelihood between two models, based on synthetic catalogues.

    The synthetic catalogues are generated based on the base model and are divided between training catalogues (for
    which the base and target model parameters are sampled) and evaluation catalogues (for which the log-likelihood
    of the generated traces is evaluated). The resulting differences in log-likelihoods are taken as those of the
    base model, minus those of the target model.

    Only stress invariant magnitude models are accepted as target and base models.
    """
    import argparse

    from .. import utils

    parser = argparse.ArgumentParser(description=_main_compare.__doc__)

    parser.add_argument('input', type=str, help='Input JSON with comparison parameters')
    parser.add_argument('--base', type=str, required=True, help='JSON containing the base model parameters')
    parser.add_argument('--target', type=str, required=True, help='JSON containing the target model parameters')
    parser.add_argument('--outdir', type=str, help='Folder to store output files (will override input file)')
    parser.add_argument('--figformat', type=str, help='Extension used for figure outputs (will override input file)')
    parser.add_argument('--dm-lims', metavar='LIM', nargs=2, type=float, default=[-3, 3],
                        help='Range of maximum magnitude difference between catalogues to select and plot')
    parser.add_argument('--dm-diffs', metavar='DIFF', nargs='+', type=float, default=[-1, 0, 1],
                        help='Maximum magnitude differences selected for trace and ologp difference plot')
    parser.add_argument('--debug', action='store_true', help='Print debug information during run')

    options = parser.parse_args()

    if options.debug:
        utils.set_logger('DEBUG')
    else:
        utils.set_logger('INFO')

    runner = ModelCompare(options.input, options.base, options.target)
    runner(options.dm_diffs, options.dm_lims, options.outdir, options.figformat)
    return 0


def _main_sample():
    """Generate a sample of model parameters."""
    import argparse

    from .. import utils

    parser = argparse.ArgumentParser(description=_main_sample.__doc__)

    parser.add_argument('input', type=str, help='Input JSON with comparison parameters')
    parser.add_argument('--model', type=str, required=True, help='JSON containing the base model parameters')
    parser.add_argument('--outdir', type=str, help='Folder to store output files (will override input file)')
    parser.add_argument('--figformat', type=str, help='Extension used for figure outputs (will override input file)')
    parser.add_argument('--debug', action='store_true', help='Print debug information during run')

    options = parser.parse_args()

    if options.debug:
        utils.set_logger('DEBUG')
    else:
        utils.set_logger('INFO')

    runner = ModelSample(options.input, options.model)
    runner(options.outdir, options.figformat)
    return 0
