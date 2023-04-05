Comparison script
=================

The script ``ologp-compare`` is used to compare the out-of-sample log-likelihood between two models.
This uses a set of synthetic catalogues, dividing them in training and evaluation events. The results of
the comparison are displayed in terms of the difference in the event maximum magnitude between the training
and evaluation catalogues.

The base and target models used for comparison *must* be stress invariant magnitude models.

Command line options
--------------------

A basic run uses:

.. code-block:: bash

   $ ologp-compare <compare.json> --base <base-model.json> --target <target-model.json>

The structure of the file ``<compare.json>`` is shown below, while ``<base-model.json>`` and
``<target-model.json>`` have the same structure, which is described in :ref:`models`.

The following command line options are optional:

* ``--outdir``: Directory where output files will be stored (overrides JSON if provided).
* ``--figformat``: Format for output figures (overrides JSON if provided). Typically one of ``png`` or ``pdf``.
* ``--dm-lims``: Lower and higher limit for maximum magnitude differences to select and plot.
* ``--dm-diffs``: Selection of maximum magnitude differences used in a few plots. Defaults to ``-1, 0, 1``.
* ``--debug``: Print debugging information during the run.

.. _comparison_json:

Input JSON
----------

The main input defining the sampling parameters is a JSON file with the following structure:

.. code-block:: json

   {
       "Mmin": 1.5,
       "Mmax": 7,
       "output_dir": "outputs",
       "figure_format": "png",
       "simulation": {
           "random_seed": null,
           "catalogues_simulate": 10000,
           "catalogues_select": 501,
           "events_training": 140,
           "events_evaluation": 140,
       },
       "mcmc": {
           "cores": 1,
           "chains": 2,
           "draws": 10000,
           "tune": 10000,
           "init": "advi+adapt_diag",
           "target_accept": 0.95
       }
   }

Some fields can be left out and default values will be provided (defaults are shown
above). An example file is included as ``inputs/compare.json``.

A brief summary of the options:

* ``Mmin``: Minimum event magnitude accepted in the model.
* ``Mmax``: Maximum event magnitude possible in simulated catalogues.
* ``output_dir``: Directory where output files will be stored.
* ``figure_format``: Format for output figures. Typically one of ``png`` or ``pdf``.
* ``simulation``: Parameters for simulating synthetic catalogues.

  * ``random_seed``: Seed used for producing catalogues (for reproducibility).
  * ``catalogues_simulate``: Number of catalogues to simulate.
  * ``catalogues_select``: Number of simulated catalogues to select for comparison.
  * ``events_training``: Number of event in training catalogues.
  * ``events_evaluation``: Number of event in evaluation catalogues catalogues.

* ``mcmc``: MCMC sampling parameters.

  * ``cores``: Number of cores used when doing the MCMC sampling.
  * ``chains``: Number of chains to sample.
  * ``draws``: Number of samples generated for each chain.
  * ``tune``: Number of MCMC tuning steps.
  * ``init``: Option passed directly to PyMC3 sampling. This is an advanced option and better left unchanged.
  * ``target_accept``: Option passed directly to PyMC3 sampling. This is an advanced option and better left unchanged.

Outputs
-------

The following output files are generated (examples below each description):

* ``trace_<base-label>.<format>``: Base model variable traces and chains at selected :math:`\Delta M_\mathrm{max}`
  values.

  .. image:: figures/trace_uni.b1.png
     :align: center
     :alt: Base model traces

* ``trace_<target-label>.<format>``: Target model variable traces and chains at selected :math:`\Delta M_\mathrm{max}`
  values.

  .. image:: figures/trace_uni.b1z1.png
     :align: center
     :alt: Target model traces

* ``dlogp_<base-label>_<target-label>.<format>``: Distribution of log-likelihood difference between target and base
  models at selected :math:`\Delta M_\mathrm{max}` values.

  .. image:: figures/dlogp_uni.b1_uni.b1z1.png
     :align: center
     :alt: Distribution of out-of-sample log-likelihood difference

* ``logp_diff_<base-label>_<target-label>.<format>``: Mean, median and 95% spread around the median of log-likelihood
  difference between target and base models as a function of :math:`\Delta M_\mathrm{max}`.

  .. image:: figures/logp_diff_uni.b1_uni.b1z1.png
     :align: center
     :alt: Log-likelihood difference as a function of maximum magnitude difference

* ``logp_ratio_<base-label>_<target-label>.<format>``: Mean log-likelihood ratio between target and base models as a
  function of :math:`\Delta M_\mathrm{max}`.

  .. image:: figures/logp_ratio_uni.b1_uni.b1z1.png
     :align: center
     :alt: Log-likelihood ratio as a function of maximum magnitude difference

* ``performance_weights_<base-label>_<target-label>.<format>``: Ratio of samples that lead to higher likelihoods in the
  base model with respect to the total number of samples as a function of :math:`\Delta M_\mathrm{max}`.

  .. image:: figures/performance_weights_uni.b1_uni.b1z1.png
     :align: center
     :alt: Performance weights as a function of maximum magnitude difference
