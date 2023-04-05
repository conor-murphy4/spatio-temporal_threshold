# This is a CONFIDENTIAL document.
#
# Copyright Royal Dutch/Shell Group of Companies
#
# Neither the whole or any part of this document may be
# reproduced, stored in any retrieval system or trans-
# mitted in any form or by any means (electronic,
# mechanical, reprographic, recording or otherwise)
# without the prior consent of the copyright owner.
"""Functions to generate different figures (using `matplotlib <https://matplotlib.org/3.3.1/index.html>`_)"""
import numpy as np
import pymc3 as pm
import scipy.stats
import matplotlib.pyplot as plt


def _select_variables(trace):
    return [var for var in trace.varnames if not var.endswith('__') and np.min(trace[var]) != np.max(trace[var])]


def traces(path, model, trace, coords=None, labels=None):
    """Plots the variable traces and chains.

    Args:
        path (str): Path to output file where figure will be saved
        model (pymc3.model.Model): Model from which the traces were generated
        trace (pymc3.backends.base.MultiTrace): Trace containing the samples
        coords (list[int], None): Catalogue indices in trace to select. Default is using all
        labels (list[str], None): Labels for the selected coordinates. Assumes same length as ``coords``
    """
    variables = _select_variables(trace)
    with model:
        kwargs_plot = dict()
        if coords:
            kwargs_plot.update(coords={var + '_dim_0': coords for var in variables})
        pm.traceplot(trace, var_names=variables, **kwargs_plot)
        if labels:
            handlers = plt.gcf().axes[-2].lines[::trace.nchains]  # second to last has the traces proper
            plt.gcf().set_constrained_layout(False)
            plt.gcf().legend(handlers, labels, loc='lower center', ncol=5, bbox_to_anchor=(.5, -.01))
            plt.tight_layout(rect=[0, 0.1, 1, 1])
    plt.savefig(path)


def posterior(path, model, trace, coords=None):
    """Plots the variable traces and chains.

    Args:
        path (str): Path to output file where figure will be saved
        model (pymc3.model.Model): Model from which the traces were generated
        trace (pymc3.backends.base.MultiTrace): Trace containing the samples
        coords (list[int], None): Catalogue indices in trace to select. Default is using all
    """
    variables = _select_variables(trace)
    with model:
        kwargs_plot = dict(hdi_prob=0.95)
        if coords:
            kwargs_plot.update(coords={var + '_dim_0': coords for var in variables})
        pm.plot_posterior(trace, var_names=variables, **kwargs_plot)
    plt.savefig(path)


def pair(path, model, trace, coords=None):
    variables = _select_variables(trace)
    with model:
        kwargs_plot = dict()
        if coords:
            kwargs_plot.update(coords={var + '_dim_0': coords for var in variables})
        pm.pairplot(trace, var_names=variables, kind='kde', **kwargs_plot)
    plt.savefig(path)


def magnitude_exceedance(path, mags, sfs, model_labels, labels, trains=None, evals=None, interval=95):
    """Plot magnitude exceedance probability.

    Plots the survival functions for each model provided in ``sfs`` in terms of the median and the prediction
    interval requested. The plotting of training and evaluation sets is optional

    Args:
        path (str): Path to output file where figure will be saved. Figure will have ``ndiffs`` subplots
        mags (np.ndarray[float]): Magnitudes for which survival functions have been calculated. Assumes all of them
            match. Shape ``(nmags,)``
        sfs (tuple[np.ndarray[float]]): List of survival functions to plot. The individual survival functions must
            have shape ``(nmags, nsamples, ndiffs)``. The median and prediction intervals are derived along samples
        model_labels (tuple[str]): Model labels used in plot legend. Assumed to have same length as ``sfs``
        labels (list[str]): Titles for each subplot. Length is assumed to be ``ndiffs``
        trains (np.ndarray[float], None): Training magnitude set. Shape ``(ntrain, ndiffs)``
        evals (np.ndarray[float], None): Evaluation magnitude set. Shape ``(neval, ndiffs)``
        interval (float): Prediction interval to plot. Checks that it is within 0 and a 100
    """
    colors_error = ['orange', 'lightblue']
    colors_median = ['darkorange', 'steelblue']
    if interval <= 0 or interval > 100:
        raise ValueError('Confidence interval must be between 0 and 100. Set to: {}'.format(interval))
    plt.figure(figsize=(5, 3 * len(labels)))
    for i, label in enumerate(labels):
        plt.subplot(len(labels), 1, i + 1)
        for j, (sf, model_label) in enumerate(zip(sfs, model_labels)):
            median = np.median(sf[..., i], axis=1)
            errors = np.percentile(sf[..., i], (50 - interval / 2, 50 + interval / 2), axis=1)
            try:
                plt.fill_between(mags, *errors, facecolor=colors_error[j], alpha=0.5)
                plt.plot(mags, median, color=colors_median[j], lw=2, label='Model: {}'.format(model_label))
            except IndexError:
                color = plt.fill_between(mags, *errors, alpha=0.5).get_facecolor()  # Catch default color
                color[-1][-1] = 1  # Set opacity to one
                plt.plot(mags, median, lw=2, label='Model: {}'.format(model_label))
        if evals is not None:
            plt.plot(*_observed_sf_stairstep(evals[:, i]), color='grey', lw=1, label='Data: Evaluation')
        if trains is not None:
            plt.plot(*_observed_sf_stairstep(trains[:, i]), color='black', lw=1, label='Data: Training')
        plt.legend()
        plt.yscale('log')
        plt.ylim(1e-4, 1)
        plt.xlim(1.5, 6)
        if i == len(labels) - 1:
            plt.xlabel('Magnitude')
        if i == len(labels) // 2:
            plt.ylabel('Exceedance probability')
        plt.title(label)
    plt.tight_layout()
    plt.savefig(path)


def log_likelihood_difference(path, dmags, dlogp, lims, interval=95):
    """Plot log-likelihood difference between models.

    Plots the difference as a function of maximum magnitude difference between training and evaluation sets. The
    median, mean and confidence interval are plotted.

    Args:
        path (str): Path to output file where figure will be saved
        dmags (np.ndarray[float]): Selected maximum magnitude differences. Shape ``(nselect,)``
        dlogp (np.ndarray[float]): Samples log-likelihood differences. Shape ``(nselect, nsamples)``
        lims (list[float or None]): Length 2 list with limits for maximum magnitude difference to plot
        interval (float): Prediction interval to plot. Checks that it is within 0 and a 100
    """
    if interval <= 0 or interval >= 100:
        raise ValueError('Confidence interval must be between 0 and 100. Set to: {}'.format(interval))
    plt.figure(figsize=(5, 2.5))
    mask = ~np.isnan(dlogp[:, 0])  # Catalogues with NaN give a full NaN slice
    mean = np.mean(dlogp[mask], axis=-1)
    median = np.median(dlogp[mask], axis=-1)
    err = np.percentile(dlogp[mask], (50 - interval / 2, 50 + interval / 2), axis=-1)
    plt.plot(dmags[mask], mean, 'o', ms=1, mfc='darkorange', mec='darkorange', label='Mean')
    plt.errorbar(dmags[mask], median, np.abs(err - median), ls='', marker='o', ms=1, mfc='steelblue',
                 mec='steelblue', ecolor='lightblue', elinewidth=1, label='Median (with {}% CI)'.format(interval))
    plt.legend(frameon=False, loc='best', fontsize=8)
    plt.xlim(lims)
    plt.ylim(-20, 50)
    plt.xlabel(r'Magnitude difference, $M_\mathrm{max, t} - M_\mathrm{max, e}$')
    plt.ylabel('Log likelihood difference,\n' + r'$\ell_1 - \ell_2$')
    plt.tight_layout()
    plt.savefig(path)


def log_likelihood_ratio(path, dmags, dlogp, lims):
    """Plot log-likelihood ratio between models.

    Plots the ratio as a function of maximum magnitude difference between training and evaluation sets. Only the
    mean is plotted.

    Args:
        path (str): Path to output file where figure will be saved
        dmags (np.ndarray[float]): Selected maximum magnitude differences. Shape ``(nselect,)``
        dlogp (np.ndarray[float]): Samples log-likelihood differences. Shape ``(nselect, nsamples)``
        lims (list[float or None]): Length 2 list with limits for maximum magnitude difference to plot
    """
    plt.figure(figsize=(5, 3.5))
    mask = ~np.isnan(dlogp[:, 0])  # Catalogues with NaN give a full NaN slice
    plt.plot(dmags[mask], np.exp(np.mean(dlogp[mask], axis=-1)), 'o', ms=1, mfc='steelblue', mec='steelblue')
    plt.vlines(0, 7.5e-2, 1e4, colors='k', lw=.5)
    plt.xlim(lims)
    plt.ylim(7.5e-2, 1e4)
    plt.yscale('log')
    plt.xlabel(r'$M_\mathrm{max}\mathrm{[training]} - M_\mathrm{max}\mathrm{[testing]}$')
    plt.ylabel('likelihood ratio no-taper/taper')
    plt.tight_layout()
    plt.savefig(path)


def performance_weights(path, dmags, dlogp, lims):
    """Plot performance weights.

    Performance weights are calculated as the fraction of samples where the base model outperforms the
    target model (i.e., the log-likelihood difference is positive). Plots the weights as a function of
    maximum magnitude difference between training and evaluation sets.

    Args:
        path (str): Path to output file where figure will be saved
        dmags (np.ndarray[float]): Selected maximum magnitude differences. Shape ``(nselect,)``
        dlogp (np.ndarray[float]): Samples log-likelihood differences. Shape ``(nselect, nsamples)``
        lims (list[float or None]): Length 2 list with limits for maximum magnitude difference to plot
    """
    plt.figure(figsize=(5, 2.5))
    mask = ~np.isnan(dlogp[:, 0])  # Catalogues with NaN give a full NaN slice
    plt.plot(dmags[mask], np.mean(dlogp[mask] >= 0, axis=-1), 'o', ms=1, mfc='steelblue', mec='steelblue')
    lims = plt.xlim(lims)  # Needed to handle None limits
    plt.hlines(.5, *lims, colors='k', lw=.5)
    plt.ylim(0, 1)
    plt.xlabel(r'Magnitude difference, $M_\mathrm{max,t} - M_\mathrm{max,e}$')
    plt.ylabel(r'Performance weight, $P_{12}$')
    plt.tight_layout()
    plt.savefig(path)


def log_likelihood_difference_distribution(path, dlogp, labels):
    """Plot the distribution in log-likelihood difference.

    Plots distribution for a selection of maximum magnitude differences. The density is estimated using Gaussian
    kernels.

    Args:
        path (str): Path to output file where figure will be saved. Figure will have ``ndiffs`` subplots
        dlogp (np.ndarray[float]): Samples log-likelihood differences. Shape ``(ndiffs, nsamples)``
        labels (list[str]): Titles for each subplot. Length is assumed to be ``ndiffs``
    """
    with plt.style.context('seaborn'):
        plt.figure(figsize=(5, 2.5 * len(labels)))
        for i, label in enumerate(labels):
            if np.sum(~np.isnan(dlogp[i])) <= 1:
                # logger.warning('Not enough valid data to produce dlogp distribution plot for {}'.format(label))
                continue
            kde = scipy.stats.gaussian_kde(dlogp[i, ~np.isnan(dlogp[i])])
            logp_range = np.linspace(np.min(dlogp[i]), np.max(dlogp[i]), 1000)
            plt.subplot(len(labels), 1, i + 1)
            plt.plot(logp_range, kde(logp_range))
            plt.tick_params(labelleft=False, left=False, right=False, bottom=False, top=False)
            plt.grid(axis='y')
            if i == len(labels) - 1:
                plt.xlabel(r'Log likelihood difference, $\ell_1 - \ell_2$')
            if i == len(labels) // 2:
                plt.ylabel(r'PDF')
            plt.title(label)
        plt.tight_layout()
        plt.savefig(path)


def _observed_sf_stairstep(catalogue):
    """Generates the observed survival function as a stepped series, ready for plotting.

    Args:
        catalogue (np.ndarray[float]): Observed event magnitudes. Shape ``(nevents,)``

    Returns:
        tuple[np.ndarray[float]]: Length 2 tuple with:

        * *np.ndarray[float]*: Stair-step array of magnitudes. Shape ``(2 * nevents - 1,)``
        * *np.ndarray[float]*: Stair-step array of survival function values. Shape ``(2 * nevents - 1,)``
    """
    catalogue = np.sort(catalogue)
    sf = (np.arange(catalogue.size)[::-1] + 1) / catalogue.size
    return np.repeat(catalogue, 2)[:-1], np.repeat(sf, 2)[1:]
