# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
import os
import logging

import pymc3 as pm

from .. import plots
from ..geomodels import ThinSheetModel
from ..parameters import SampleParameters, ModelParameters

logger = logging.getLogger('runners')


class ModelSample:
    def __init__(self, path_sample, path_model):
        """Initialize required seismological model and read in parameters for sampling.

        Args:
            path_sample (str): Path to JSON file with sampling parameters
            path_model (str): Path to JSON file with model parameters
        """
        self._params = SampleParameters(path_sample)
        thinsheet = ThinSheetModel(self._params.paths, self._params.kernels, self._params.sigmas)

        mparams = ModelParameters(path_model, Mmin=self._params.Mmin)
        self._dof = mparams.dof
        self._generator = mparams.model_generator(thinsheet)

    def __call__(self, outdir=None, figformat=None):
        """Run sampling based on parameters read at initialization.

        This generates a series of output files:

        * ``SeismologicalModel_{type}.{label}_{start}-{end}.csv``: Start and end dates correspond to training period
        * ``posterior_{label}.{format}``: Plot with posterior distribution of sampled variables
        * ``joint_{label}.{format}``: Plot with joint distribution of sampled variable pairs

        Args:
            outdir (str, None): Path to output directory. If not provided, uses the one defined in parameters
            figformat (str, None): Format for figure files. If not provided, uses the one defined in parameters
        """
        outdir = outdir or self._params.outdir
        figformat = figformat or self._params.format

        catalogue, start, end = self._params.training_catalogue, self._params.start, self._params.end
        model = self._generator.create_model(self._dof, catalogue, start, end)
        trace = pm.sample(model=model, return_inferencedata=False, **self._params.mcmc)

        os.makedirs(outdir, exist_ok=True)

        outcsv = self._generator.export_trace(outdir, trace, start, end)
        logger.info('Generated CSV table with traces: {}'.format(outcsv))

        figname = os.path.join(outdir, 'posterior_{}.{}'.format(self._generator.label, figformat))
        plots.posterior(figname, model, trace)
        logger.info('Generated plot: {}'.format(figname))

        figname = os.path.join(outdir, 'joint_{}.{}'.format(self._generator.label, figformat))
        try:
            plots.pair(figname, model, trace)
        except Exception:  # Unfortunately this is the exception raised by ArviZ for multiple parameters
            logger.info('Skipping joint distribution plot: Single variable model')
        else:
            logger.info('Generated plot: {}'.format(figname))
