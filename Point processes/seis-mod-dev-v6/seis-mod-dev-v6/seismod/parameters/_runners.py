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
from datetime import datetime

from .._bases import parameters
from ..geomodels import EarthquakeCatalogue


class SampleParameters(parameters.Parameters):
    """Class to read sampling and evaluation parameters from a JSON file."""
    _def_params = dict(
        Mmin=1.5,
        training=dict(
            start=[1995, 1, 1],
            end=[2019, 1, 1]
        ),
        evaluation=dict(
            start=[2013, 1, 1],
            end=[2019, 1, 1]
        ),
        output_dir='outputs',
        figure_format='png',
        mcmc=dict(
            cores=1,
            chains=2,
            draws=10000,
            tune=10000,
            init='advi+adapt_diag',
            target_accept=0.95
        ),
        smoothing_sigmas=dict(
            pressure=None,
            compressibility=None,
            thickness=None
        ),
        median_kernel_sizes=dict(
            pressure=None,
            compressibility=None,
            thickness=None
        ),
        paths=dict(
            catalogue='',
            faults='',
            pressure='',
            compressibility='',
            thickness=''
        )
    )

    def __init__(self, path, **kwargs):
        r"""Accepted keys in JSON file are:

        * ``Mmin``: Minimum magnitude of accepted events
        * ``training``: Time range for selecting training events:

          * ``start``: List with year, month and day of starting date
          * ``end``:  List with year, month and day of final date

        * ``evaluation``: Time range for selecting evaluation events (analogous to ``training``)

          * ``start``: List with year, month and day of starting date
          * ``end``:  List with year, month and day of final date

        * ``output_dir``: Directory where output files are stored relative to running directory (will be created if
          it does not exist)
        * ``figure_format``: Format of output figures (usually either ``png`` or ``pdf``)
        * ``mcmc``: MCMC parameters used when sampling:

          * ``cores``: Number of cores used when sampling. There are known issues in Windows when using more than one
          * ``chains``: Number of chains sampled
          * ``draws``: Number of samples to draw for each chain
          * ``tune``: Number of tuning steps for each chain. These are samples at the beginning, which are discarded
          * ``init``: Initialization method for NUTS sampler (see `documentation
            <https://docs.pymc.io/api/inference.html#pymc3.sampling.sample>`_)
          * ``target_accept``: Step size is tuned so that it approximates this value (see `documentation
            <https://docs.pymc.io/api/inference.html#pymc3.sampling.sample>`_)

        * ``smoothing_sigmas``: Length-scale used for Gaussian smoothing of reservoir grids (optional):

          * ``pressure``: Pressure grid :math:`\sigma`
          * ``compressibility``: Compressibility grid :math:`\sigma`
          * ``thickness``: Thickness grid :math:`\sigma`

        * ``median_kernel_sizes``: Kernel size used for median filtering reservoir grids (optional):

          * ``pressure``: Kernel size for pressure median filtering
          * ``compressibility``: Kernel size for compressibility median filtering
          * ``thickness``: Kernel size for thickness median filtering

        * ``paths``: Paths to observed catalogue, reservoir grids and fault model (**important:** These paths must be
          relative to the location of the input JSON):

          * ``catalogue``: Path to observed catalogue (see :meth:`seismod.geomodels.EarthquakeCatalogue.from_csv`
            for details on the expected structure)
          * ``faults``: Path to fault model (see :meth:`seismod.geomodels.FaultModel.from_csv` for details on the
            expected structure)
          * ``pressure``: Path to pressure grid (see :meth:`seismod.geomodels.ReservoirGrid.from_csv` for details on
            the expected structure). Pressure must be given in bar
          * ``compressibility``: Path to compressibility grid (see :meth:`seismod.geomodels.ReservoirGrid.from_csv`
            for details on the expected structure). Compressibility must be given in MPa\ :sup:`-1`
          * ``thickness``: Path to thickness grid (see :meth:`seismod.geomodels.ReservoirGrid.from_csv` for details
            on the expected structure). Thickness must be given in meters

        Args:
            path (str): Path to JSON with input parameters
            **kwargs: Keyword arguments can be used to override parameters from input file
        """
        super().__init__(path, **kwargs)
        for key in self.paths:
            if self.paths[key]:
                self.paths[key] = os.path.join(self._dirname, self.paths[key])

    @property
    def Mmin(self):
        """float: Minimum magnitude accepted for models"""
        return self._access('Mmin')

    @property
    def format(self):
        """str: Format used for output figures"""
        return self._access('figure_format')

    @property
    def outdir(self):
        """str: Directory where output files will be created"""
        return self._access('output_dir')

    @property
    def training_catalogue(self):
        """seismod.geomodels.EarthquakeCatalogue: Catalogue of training events"""
        catalogue = EarthquakeCatalogue.from_csv(self.paths['catalogue'])
        catalogue.select(Mmin=self.Mmin, tmin=self.start, tmax=self.end)
        return catalogue

    @property
    def evaluation_catalogue(self):
        """seismod.geomodels.EarthquakeCatalogue: Catalogue of evaluation events"""
        catalogue = EarthquakeCatalogue.from_csv(self.paths['catalogue'])
        catalogue.select(Mmin=self.Mmin, tmin=self.eval_start, tmax=self.eval_end)
        return catalogue

    @property
    def start(self):
        """datetime.datetime: Starting date of training period"""
        return datetime(*self._access('training/start'))

    @property
    def end(self):
        """datetime.datetime: End date of training period"""
        return datetime(*self._access('training/end'))

    @property
    def eval_start(self):
        """datetime.datetime: Starting date of evaluation period"""
        return datetime(*self._access('evaluation/start'))

    @property
    def eval_end(self):
        """datetime.datetime: End date of evaluation period"""
        return datetime(*self._access('evaluation/end'))

    @property
    def mcmc(self):
        """dict[str, int, float]: Dictionary with parameters for MCMC sampling"""
        return self._access('mcmc')

    @property
    def sigmas(self):
        """dict[str, float, None]: Dictionary with length-scales for Gaussian smoothing of reservoir grids"""
        return self._access('smoothing_sigmas')

    @property
    def kernels(self):
        """dict[str, float, None]: Dictionary with kernel sizes for median filtering of reservoir grids"""
        return self._access('median_kernel_sizes')

    @property
    def paths(self):
        """dict[str]: Dictionary with paths to observed catalogue, reservoir grids and fault model"""
        return self._access('paths')


class CompareParameters(SampleParameters):
    """Class to read running comparison runs parameters.

    Additional keys accepted in JSON:

    * ``Mmax``: Maximum magnitude for simulated events
    * ``simulation``: Parameters used for generating synthetic training and evaluation catalogues

      * ``seed``: Seed used to generate catalogues for reproducibility (optional)
      * ``catalogues_simulate``: Number of catalogues to simulate
      * ``catalogue_select``: Number of simulated catalogues selected for the comparison
      * ``events_training``: Number of events to simulate in each training catalogues
      * ``events_evaluation``: Number of events to simulate in each evaluation catalogues
    """
    _def_params = dict(
        Mmax=7,
        simulation=dict(
            seed=None,
            catalogues_simulate=100000,
            catalogues_select=501,
            events_training=280,
            events_evaluation=280,
        )
    )

    @property
    def Mmax(self):
        """float: Maximum magnitude accepted in model"""
        return self._access('Mmax')

    @property
    def simulation(self):
        """dict[str, int, None]: Parameters used for simulation of training and evaluation catalogues"""
        return self._access('simulation')
