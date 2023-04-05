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

import numpy as np
import pymc3 as pm

from .. import plots
from ..geomodels import ThinSheetModel, EarthquakeCatalogue
from ..parameters import CompareParameters, ModelParameters

logger = logging.getLogger('runners')


class ModelCompare:
    def __init__(self, path_compare, path_base, path_target):
        """Initialize required seismological models and read in parameters for comparison.

        Args:
            path_compare (str): Path to JSON file with comparison parameters
            path_base (str): Path to JSON file with parameters for base model
            path_target (str): Path to JSON file with parameters for target model
        """
        self._params = CompareParameters(path_compare)
        self._bparams = ModelParameters(path_base, Mmin=self._params.Mmin)
        if not self._bparams.is_magnitude or not self._bparams.is_invariant:
            raise ValueError('Base model must be an invariant magnitude model')
        self._tparams = ModelParameters(path_target, Mmin=self._params.Mmin)
        if not self._tparams.is_magnitude or not self._tparams.is_invariant:
            raise ValueError('Target model must be an invariant magnitude model')
        self._thinsheet = ThinSheetModel(self._params.paths, self._params.kernels, self._params.sigmas)

    def _simulate_events(self, dmmin, dmmax):
        base = self._bparams.model_generator(self._thinsheet)

        nc, ns = self._params.simulation['catalogues_simulate'], self._params.simulation['catalogues_select']
        nt, ne = self._params.simulation['events_training'], self._params.simulation['events_evaluation']

        rng = np.random.default_rng(seed=self._params.simulation['seed'])
        tcatalogue = EarthquakeCatalogue(time=np.empty((nt, nc)), location=np.empty((nt, nc, 3)))
        ecatalogue = EarthquakeCatalogue(time=np.empty((ne, nc)), location=np.empty((ne, nc, 3)))

        base.simulate(var=self._bparams.values, catalogue=tcatalogue, rng=rng)
        base.simulate(var=self._bparams.values, catalogue=ecatalogue, rng=rng)

        dm = np.max(tcatalogue.magnitude, axis=0) - np.max(ecatalogue.magnitude, axis=0)
        order = np.argsort(dm)
        if nc > ns:
            logger.info('Selecting {} catalogues'.format(ns))
            dmmin, dmmax = dmmin if dmmin else np.min(dm), dmmax if dmmax else np.max(dm)
            selection = np.clip(np.searchsorted(dm[order], np.linspace(dmmin, dmmax, ns)), 0, nc - 1)
            order = order[selection]
        else:
            logger.info('Using all catalogues')
        tcatalogue.select(trial=order)
        ecatalogue.select(trial=order)
        return tcatalogue, ecatalogue

    def __call__(self, selected, lims=(None, None), outdir=None, figformat=None):
        r"""Generates synthetic catalogues from base model and run sampling and comparison of out-of-sample
        log-likelihood from both the base and target model.

        The comparisons per catalogue are organized based on the maximum magnitude difference between the training
        events and the evaluation events within the catalogue, :math:`\Delta M_\mathrm{max}`.

        This generates a series of output files:

        * ``trace_{label}.{format}``: Variable traces and chains at selected :math:`\Delta M_\mathrm{max}` values
        * ``dlogp_{base}_{target}.{format}``: Distribution of log-likelihood difference between target and base
          models at selected :math:`\Delta M_\mathrm{max}` values
        * ``logp_diff_{base}_{target}.{format}``: Mean, median and 95% spread around the median of log-likelihood
          difference between target and base models as a function of :math:`\Delta M_\mathrm{max}`
        * ``logp_ratio_{base}_{target}.{format}``: Mean log-likelihood ratio between target and base models as a
          function of :math:`\Delta M_\mathrm{max}`
        * ``performance_weights_{base}_{target}.{format}``: Ratio of samples that lead to higher likelihoods in the
          base model with respect to the total number of samples as a function of :math:`\Delta M_\mathrm{max}`

        Args:
            selected (list[float], tuple[float]): Selected values of :math:`\Delta M_\mathrm{max}` used in plots
            lims (tuple[float, None], list[float, None]): :math:`\Delta M_\mathrm{max}` limits to select and plot
            outdir (str, None): Path to output directory. If not provided, uses the one defined in parameters
            figformat (str, None): Format for figure files. If not provided, uses the one defined in parameters
        """
        outdir = outdir or self._params.outdir
        figformat = figformat or self._params.format

        tcatalogue, ecatalogue = self._simulate_events(*lims)
        dm = tcatalogue.Mmax - ecatalogue.Mmax

        base, target = self._bparams.model_generator(self._thinsheet), self._tparams.model_generator(self._thinsheet)
        bmodel = base.create_model(self._bparams.dof, tcatalogue)
        tmodel = target.create_model(self._tparams.dof, tcatalogue)

        btrace = pm.sample(model=bmodel, return_inferencedata=False, **self._params.mcmc)
        ttrace = pm.sample(model=tmodel, return_inferencedata=False, **self._params.mcmc)

        dlogp = base.ologp(ecatalogue, btrace) - target.ologp(ecatalogue, ttrace)
        indices, labels = self._indices_and_labels(dm, selected)

        os.makedirs(outdir, exist_ok=True)
        figname = os.path.join(outdir, 'trace_{{}}.{}'.format(figformat))
        plots.traces(figname.format(base.label), bmodel, btrace, coords=indices, labels=labels)
        logger.info('Generated plot: {}'.format(figname.format(base.label)))
        plots.traces(figname.format(target.label), tmodel, ttrace, coords=indices, labels=labels)
        logger.info('Generated plot: {}'.format(figname.format(target.label)))

        figname = os.path.join(outdir, 'dlogp_{}_{}.{}'.format(base.label, target.label, figformat))
        plots.log_likelihood_difference_distribution(figname, dlogp[indices], labels)
        logger.info('Generated plot: {}'.format(figname))

        figname = os.path.join(outdir, 'logp_diff_{}_{}.{}'.format(base.label, target.label, figformat))
        plots.log_likelihood_difference(figname, dm, dlogp, lims)
        logger.info('Generated plot: {}'.format(figname))

        figname = os.path.join(outdir, 'logp_ratio_{}_{}.{}'.format(base.label, target.label, figformat))
        plots.log_likelihood_ratio(figname, dm, dlogp, lims)
        logger.info('Generated plot: {}'.format(figname))

        figname = os.path.join(outdir, 'performance_weights_{}_{}.{}'.format(base.label, target.label, figformat))
        plots.performance_weights(figname, dm, dlogp, lims)
        logger.info('Generated plot: {}'.format(figname))

    @staticmethod
    def _indices_and_labels(dm, selected):
        """Picks the indices and generate labels for the magnitude deltas closest to the selected difference(s).

        Args:
            dm (np.ndarray): Existing maximum magnitude differences. Shape ``(nselect,)``
            selected (list[float]): List of maximum magnitude differences to be selected. Length ``ndiffs``

        Returns:
            tuple: Length 2 tuple with:

            * *list[int]*: List of indices closest to requested magnitude deltas. Length ``ndiffs``
            * *list[str]*: List of labels associated to the selection (used within ``matplotlib``). Length ``ndiffs``
        """
        ind = [np.argmin(np.abs(dm - diff)) for diff in selected]
        labels = [r'$\Delta M_\mathrm{{max}} = {:.1f}$'.format(dm) for dm in dm[ind]]
        return ind, labels
