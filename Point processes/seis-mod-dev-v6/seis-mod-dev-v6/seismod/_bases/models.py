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
import abc
import logging

import numpy as np
import pymc3 as pm
import theano.tensor

from .. import distributions, utils
from ..geomodels import EarthquakeCatalogue
from .distributions import BaseDistribution

logger = logging.getLogger('models')


class Priors:
    UNI = 'uniform'
    EXP = 'exponential'
    UNILOG = 'uniform-log'
    CONST = 'constant'

    WITH_LOWER = UNI, EXP, UNILOG
    WITH_UPPER = UNI, UNILOG
    WITH_VALUE = CONST,
    WITH_MEAN = EXP,


class Options:
    PRIOR = 'distr'
    VALUE = 'value'
    LOWER = 'lower'
    UPPER = 'upper'
    MEAN = 'mean'


class _BaseModel(abc.ABC):
    """Base class for model building.

    The distribution used, expected variables and used covariates must be specified in sub-classes.
    """
    ObsDistr = BaseDistribution
    _type = ''
    _vars = ('Mmin',)
    _covs = ()

    def __init__(self, thinsheet, label=''):
        """Details on the distribution used for observed parameters can be found in :py:attr:`ObsDistr`.

        Args:
            thinsheet (seismod.geomodels.ThinSheetModel): Reservoir properties
            label (str): Label used to identify the model
        """
        self._thinsheet = thinsheet
        self._label = label
        logger.info('Initializing {} model (label: {})'.format(self.__class__.__name__, self._label))
        logger.debug('Using PyMC3 version {}'.format(pm.__version__))

    @property
    def label(self):
        """str: Label identifying the model"""
        return self._label

    def ologp(self, evals, trace, cutoff=None):
        """Calculate the log-likelihood of the given earthquake catalogue for each parameter sample in the trace.

        Args:
            evals (seismod.geomodels.EarthquakeCatalogue): Catalogue for which to calculate the log-likelihood
            trace (pymc3.backends.base.MultiTrace): Trace values on which to condition the log-likelihood
            cutoff (float, None): Earthquakes with magnitudes below this value will not be included in computation

        Returns:
            float, np.ndarray[float]: Log-likelihood for each sample. Shape ``(nsamples,)``
        """
        cutoff = cutoff if cutoff else -np.inf
        logger.info('Computing out-of-sample log-likelihood ({})...'.format(self._label))
        out = []
        for var in self._iterate_trace(trace):
            covariates = self._generate_all_covariates(var, evals.time, evals.location)
            parameters = self._generate_all_parameters(var, covariates)
            values = np.ma.masked_less_equal(self._observations(evals, var), cutoff)
            out.append(self.ObsDistr.loglike(values, **parameters))
        return np.ma.filled(np.ma.array(out), np.nan).T

    def export_trace(self, outdir, trace, start, end):
        """Exports a CSV with the trace for each variable in model (``Mmin`` is always skipped).

        The output file name is generated based on the model parameters and arguments, it is always set as
        ``SeismologicalModel_{type}.{label}_{start}-{end}.csv``. The type is normally either ``mom`` or ``ar``
        for moment or activity rate models, the label is the model label provided at initialization and the start
        and end date provided are assumed to correspond to the training period, although this is not verified.

        Args:
            outdir (str): Base directory where output file will be stored
            trace (pymc3.backends.base.MultiTrace): Trace with variables sampled values
            start (datetime.datetime): Start date for file labeling (usually training period)
            end (datetime.datetime): End date for file labeling (usually training period)
        """
        path = 'SeismologicalModel_{}.{}_{:%Y%m%d}-{:%Y%m%d}.csv'.format(self._type, self.label, start, end)
        path = os.path.join(outdir, path)

        nsamples = len(trace) * trace.nchains
        variables = [name for name in sorted(self._vars) if name != 'Mmin']
        output = np.empty((nsamples, len(variables)))
        for i, name in enumerate(variables):
            output[:, i] = trace[name]
        utils.export_csv(path, output, header=','.join(variables))
        return path

    @classmethod
    def _check_dof(cls, dof):
        """Checks whether the dictionary with the parametrization for the model variables is consistent.

        The parametrization is expected to be a nested dictionary, with first level keys containing all model
        variables. The second level dictionary must contain the distribution used for the variable and the parameters
        needed for such distribution (e.g., if the distribution is constant, it requires that lower and upper bounds
        are defined).

        Args:
            dof (dict[str, dict[str, float]]): Parametrization of the model's degrees-of-freedom
        """
        if set(dof) - set(cls._vars):
            raise ValueError('Invalid DOF(s) {} for {}'.format(set(dof) - set(cls._vars), cls.__name__))
        if set(cls._vars) - set(dof):
            raise ValueError('Missing DOF(s) {} for {}'.format(set(cls._vars) - set(dof), cls.__name__))
        for var, params in dof.items():
            if params[Options.PRIOR] in Priors.WITH_MEAN and Options.MEAN not in params:
                raise ValueError('Missing parameter \'mean\' for DOF {}'.format(var))
            if params[Options.PRIOR] in Priors.WITH_LOWER and Options.LOWER not in params:
                raise ValueError('Missing parameter \'lower\' for DOF {}'.format(var))
            if params[Options.PRIOR] in Priors.WITH_VALUE and Options.VALUE not in params:
                raise ValueError('Missing parameter \'value\' for DOF {}'.format(var))
            if params[Options.PRIOR] in Priors.WITH_MEAN and params[Options.PRIOR] in Priors.WITH_LOWER:
                if params[Options.MEAN] < params[Options.LOWER]:
                    raise ValueError('Mean value for {} is below lower bound'.format(var))
            if params[Options.PRIOR] in Priors.WITH_MEAN and params[Options.PRIOR] in Priors.WITH_UPPER:
                if params[Options.MEAN] > params[Options.UPPER]:
                    raise ValueError('Mean value for {} is above upper bound'.format(var))

    def create_model(self, dof, catalogue, start=None, end=None):
        """Creates a PyMC3 model.

        The parametrization for the degrees-of-freedom is expected to be a nested dictionary, with first level keys
        containing all model variables. The second level dictionary must contain the distribution used for the variable
        and the parameters needed for such distribution (e.g., if the distribution is constant, it requires that lower
        and upper bounds are defined).

        Args:
            dof (dict[str, dict[str, float]]): Parametrization of the model's degrees-of-freedom
            catalogue (seismod.geomodels.EarthquakeCatalogue): Catalogue used for observed values
            start (datetime.datetime, None): Model starting date
            end (datetime.datetime, None): Model ending date

        Returns:
            pymc3.model.Model: PyMC3 model generated based on the input parameters
        """
        self._check_dof(dof)
        try:
            _, shape = catalogue.magnitude.shape
        except ValueError:
            shape = ()
        if start is not None and np.any(catalogue.time < start):
            raise ValueError('Cannot handle event times before period start')
        if end is not None and np.any(catalogue.time > end):
            raise ValueError('Cannot handle event times before period end')
        with pm.Model() as model:
            variables = self.__generate_all_variables(dof, shape)
            covariates = self._generate_all_covariates(variables, catalogue, start, end)
            parameters = self._generate_all_parameters(variables, covariates)
            self.ObsDistr(name='obs', observed=self._observations(catalogue, variables), shape=shape, **parameters)
        return model

    @abc.abstractmethod
    def _observations(self, catalogue, var):
        """Generate the observed value for the underlying model distribution.

        If generating them requires the use of variables, the output will be a PyMC3 distribution.

        Args:
            catalogue (seismod.geomodels.EarthquakeCatalogue): Catalogue used for observed values
            var (dict[str, np.ndarray[float], pymc3.distributions.Distribution]): Model variables

        Returns:
            np.ndarray[float], pymc3.distributions.Distribution: Observed parameter
        """
        pass

    def __generate_all_variables(self, dof, shape):
        """Create the distribution of all model variables.

        Args:
            dof (dict[str, dict[str, float]]): Parametrization of the model's degrees-of-freedom
            shape (int): Shape of the prior for the variable (typically amounts to the number of catalogues)

        Returns:
            dict[str, np.ndarray[float], pymc3.distributions.Distribution]: Distributions for all model variables
        """
        return {name: self.__generate_variable(dof[name], name, shape) for name in self._vars}

    @staticmethod
    def __generate_variable(par, name, shape):
        """Creates the distribution for the variable with the given name.

        Args:
            par (dict[str, float]): Parametrization of the variable distribution
            name (str): Name of the variable for which to generate the prior
            shape (int): Shape of the prior for the variable (typically amounts to the number of catalogues)

        Returns:
            np.ndarray[float], pymc3.distributions.Distribution: Distribution for the variable
        """
        prior, value = par.get(Options.PRIOR), par.get(Options.VALUE)
        lower, upper, mean = par.get(Options.LOWER), par.get(Options.UPPER), par.get(Options.MEAN)

        msg = ' Parameter {:8s}: {:12s}'.format(name, prior)
        if par[Options.PRIOR] == Priors.UNI:
            logger.info('{} (lower={}; upper={})'.format(msg, lower, upper))
            return pm.Uniform(name=name, lower=lower, upper=upper, shape=shape)
        if par[Options.PRIOR] == Priors.UNILOG:
            logger.info('{} (lower={}; upper={})'.format(msg, lower, upper))
            uniform = pm.Uniform(name='{}_log__'.format(name), lower=np.log(lower), upper=np.log(upper), shape=shape)
            return pm.Deterministic(name=name, var=theano.tensor.exp(uniform))
        if par[Options.PRIOR] == Priors.EXP:
            logger.info('{} (lower={}; mean={})'.format(msg, lower, mean))
            if par[Options.LOWER] == 0:
                return pm.Exponential(name=name, lam=1 / mean, shape=shape)
            return pm.Bound(pm.Exponential, lower=lower)(name=name, lam=1 / mean, shape=shape)
        if par[Options.PRIOR] == Priors.CONST:
            logger.info('{} (value={})'.format(msg, value))
            pm.Deterministic(name=name, var=theano.tensor.constant(par[Options.VALUE]))
            return par[Options.VALUE]
        raise NotImplementedError('Distribution {} not found'.format(par[Options.PRIOR]))

    def _generate_all_covariates(self, var, catalogue=None, start=None, end=None):
        """Creates all covariates used in model.

        Args:
            var (dict[str, np.ndarray[float], pymc3.distributions.Distribution]): Model variables
            catalogue (seismod.geomodels.EarthquakeCatalogue): Catalogue used for observed values
            start (datetime.datetime, None): Start date for model
            end (datetime.datetime, None): End date for model

        Returns:
            dict[str, np.ndarray[float], pymc3.distributions.Distribution]: Model covariates
        """
        return {name: self._generate_covariate(name, var, catalogue, start, end) for name in self._covs}

    @abc.abstractmethod
    def _generate_covariate(self, name, var, catalogue, start, end):
        """Creates the covariate with the given name.

        Args:
            name (str): Name of the covariate to generate
            var (dict[str, np.ndarray[float], pymc3.distributions.Distribution]): Model variables
            catalogue (seismod.geomodels.EarthquakeCatalogue): Catalogue used for observed values
            start (datetime.datetime, None): Start date for model
            end (datetime.datetime, None): Start date for model

        Returns:
            float, np.ndarray[float], pymc3.distributions.Distribution: Required covariate
        """
        raise ValueError('Model {} cannot handle covariate {}'.format(self.__class__.__name__, name))

    def _generate_all_parameters(self, var, covs):
        """Creates all parameters used in the observed distribution.

        Args:
            var (dict[str, np.ndarray[float], pymc3.distributions.Distribution]): Model variables
            covs (dict[str, np.ndarray[float], pymc3.distributions.Distribution]): Model covariates

        Returns:
            dict[str, np.ndarray[float], pymc3.distributions.Distribution]: All observed distribution parameters
        """
        return {name: self._generate_parameter(name, var, covs) for name in self.ObsDistr.parameters}

    @abc.abstractmethod
    def _generate_parameter(self, name, var, covs):
        """Creates the observed distribution parameter with the given name.

        Args:
            name (str): Name of observed distribution parameter to generate
            var (dict[str, np.ndarray[float], pymc3.distributions.Distribution]): Model variables
            covs (dict[str, np.ndarray[float], pymc3.distributions.Distribution]): Model covariates

        Returns:
            float, np.ndarray[float], pymc3.distributions.Distribution: Observed distribution parameter
        """
        raise ValueError('Model {} cannot handle parameter {}'.format(self.__class__.__name__, name))

    def _iterate_trace(self, trace):
        """Iterates over the given trace.

        Assumes that the trace is compatible with the model (i.e., the variables match)

        Args:
            trace (pymc3.backends.base.MultiTrace): Trace values to iterate over

        Yields:
            dict[str, float, np.ndarray[float]]: Current variable set in the trace
        """
        ndraws = len(trace) * trace.nchains
        for i in range(ndraws):
            yield {name: trace[name][i] for name in self._vars}

    @abc.abstractmethod
    def simulate(self, var, rng=None, **kwargs):
        pass


class StressInvariantMagnitudeModel(_BaseModel, abc.ABC):
    ObsDistr = distributions.PowerLawExponentialTaper
    _type = 'mom'

    def _observations(self, catalogue, var):
        return catalogue.magnitude

    def _generate_covariate(self, name, var, catalogue, start, end):
        return super()._generate_covariate(name, var, catalogue, start, end)

    def simulate(self, var, rng=None, **kwargs):
        """Simulate magnitudes for the given catalogue. Catalogue is modified in-place.

        Args:
            var (dict[str, np.ndarray[float]]): Model variables to use for simulation
            rng (np.random.Generator, None): Random number generator used to generate the samples

        Keyword Args:
            catalogue (seismod.geomodels.EarthquakeCatalogue): Catalogue with times and locations
        """
        catalogue = kwargs.pop('catalogue')
        covariates = self._generate_all_covariates(var, catalogue=catalogue)
        parameters = self._generate_all_parameters(var, covariates)
        catalogue.magnitude[...] = self.ObsDistr.rvs(size=catalogue.nevents, rng=rng, **parameters)


class StressInvariantActivityRateModel(_BaseModel, abc.ABC):
    ObsDistr = distributions.SpaceTime
    _covs = ('uniform_grid',)
    _type = 'ar'

    def _generate_all_covariates(self, var, catalogue=None, start=None, end=None, full=False):
        return {name: self._generate_covariate(name, var, catalogue, start, end, full) for name in self._covs}

    def _generate_covariate(self, name, var, catalogue, start, end, full=False):
        if name == 'uniform_grid':
            return self._thinsheet.uniform_grid(start, end, full)
        return super()._generate_covariate(name, var, catalogue, start, end)

    def simulate(self, var, rng=None, **kwargs):
        """Simulate the time and location of a given number of events (alternatively, randomly sample it as well).

        Args:
            var (dict[str, np.ndarray[float]]): Model variables to use for simulation
            rng (np.random.Generator, None): Random number generator used to generate the samples

        Keyword Args:
            nevents (int, tuple[int]): Number of events. If not given, it will be randomly sampled
            depth (float): Depth for all events. If not given, it defaults to 3000 m
            start (datetime.datetime): Starting date for simulation (required)
            end (datetime.datetime): End date for simulation (required)

        Returns:
            tuple[np.ndarray]: Times (`np.datetime64`, shape ``(nevents,)``) and locations (`float`,
            shape ``(nevents, 3)``) of the generated events. Locations are in meters in RD coordinates
        """
        nevents, depth = kwargs.pop('nevents', None), kwargs.pop('depth', 3000.)
        start, end = kwargs.pop('start'), kwargs.pop('end')
        covariates = self._generate_all_covariates(var, start=start, end=end, full=True)
        parameters = self._generate_all_parameters(var, covariates)

        nevents = nevents if nevents else self.ObsDistr.nevents(rng, **parameters)
        it, ix, iy = self.ObsDistr.rvs(size=nevents, rng=rng, **parameters)
        time, location = self._thinsheet.relcoords2abscoords(start, end, it, ix, iy, z=depth)
        return EarthquakeCatalogue(time=time, location=location)


class StressDependentMagnitudeModel(StressInvariantMagnitudeModel, abc.ABC):
    _vars = StressInvariantMagnitudeModel._vars + ('beta2', 'beta3', 'beta4')
    _covs = ('stress_event',)

    def _generate_covariate(self, name, var, catalogue, start, end):
        if name == 'stress_event':
            b2, b3, b4 = var['beta2'], var['beta3'], var['beta4']
            return self._thinsheet.stress_event(catalogue.time, catalogue.location, b2, b3, b4, nu=0)
        return super()._generate_covariate(name, var, catalogue, start, end)


class StressDependentActivityRateModel(StressInvariantActivityRateModel, abc.ABC):
    _vars = StressInvariantActivityRateModel._vars + ('beta2', 'beta3', 'beta4')
    _covs = ('stress_grid', 'thickness_grid', 'magnitude_event')

    def _generate_covariate(self, name, var, catalogue, start, end, full=False):
        if name == 'stress_grid':
            b2, b3, b4 = var['beta2'], var['beta3'], var['beta4']
            return self._thinsheet.stress_grid(start, end, b2, b3, b4, full)
        if name == 'thickness_grid':
            return self._thinsheet.thickness_grid()
        if name == 'magnitude_event':
            return catalogue.magnitude
        return super()._generate_covariate(name, var, catalogue, start, end)
