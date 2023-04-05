Sampling script
===============

The script ``seismod-sample`` is used to generate sets of parameters, the result of which are stored in
a CSV file. Figures showing the output distributions are generated as well. It is intended to be run
from the command line and directly uses two input files.

Command line options
--------------------

A basic run uses:

.. code-block:: bash

   $ seismod-sample <sample.json> --model <model.json>

The structure of the file ``<sample.json>`` is shown below, while ``<model.json>`` is described in
:ref:`models`.

The following command line options are optional:

* ``--outdir``: Directory where output files will be stored (overrides JSON if provided).
* ``--figformat``: Format for output figures (overrides JSON if provided). Typically one of ``png`` or ``pdf``.
* ``--debug``: Print debugging information during the run.

.. _sampling_json:

Input JSON
----------

The main input defining the sampling parameters is a JSON file with the following structure:

.. code-block:: json

   {
       "Mmin": 1.5,
       "output_dir": "outputs",
       "figure_format": "png",
       "training": {
           "start": [1995, 1, 1],
           "end": [2019, 1, 1]
       },
       "mcmc": {
           "cores": 1,
           "chains": 2,
           "draws": 10000,
           "tune": 10000,
           "init": "advi+adapt_diag",
           "target_accept": 0.95
       },
       "smoothing_sigmas": {
           "pressure": null,
           "compressibility": null,
           "thickness": null
       },
       "median_kernel_sizes": {
           "pressure": null,
           "compressibility": null,
           "thickness": null
       },
       "paths": {
           "catalogue": "",
           "faults": "",
           "pressure": "",
           "compressibility": "",
           "thickness": ""
       }
   }

Some fields can be left out and default values will be provided (defaults are shown
above). An example file is included as ``inputs/sample.json``.

A brief summary of the options:

* ``Mmin``: Minimum event magnitude accepted in the model.
* ``output_dir``: Directory where output files will be stored.
* ``figure_format``: Format for output figures. Typically one of ``png`` or ``pdf``.
* ``training``: Definition of the training period.

  * ``start``: Starting date as a list with year, month and day.
  * ``end``: Ending date as a list with year, month and day.

* ``mcmc``: MCMC sampling parameters.

  * ``cores``: Number of cores used when doing the MCMC sampling.
  * ``chains``: Number of chains to sample.
  * ``draws``: Number of samples generated for each chain.
  * ``tune``: Number of MCMC tuning steps.
  * ``init``: Option passed directly to PyMC3 sampling. This is an advanced option and better left unchanged.
  * ``target_accept``: Option passed directly to PyMC3 sampling. This is an advanced option and better left unchanged.

* ``paths``: Paths to additional input CSVs. **Paths must be absolute or relative to the JSON file.**

  * ``catalogue``: Earthquake catalogue. See :meth:`seismod.geomodels.EarthquakeCatalogue.from_csv` for the structure.
  * ``faults``: Fault model. See :meth:`seismod.geomodels.FaultModel.from_csv` for the structure.
  * ``pressure``: Reservoir pressure in bar for all epochs. See :meth:`seismod.geomodels.ReservoirGrid.from_csv` for
    the structure.
  * ``compressibility``: Reservoir compressibility in MPa\ :sup:`-1`\ . See
    :meth:`seismod.geomodels.ReservoirGrid.from_csv` for the structure.
  * ``thickness``: Reservoir thickness in meters. See :meth:`seismod.geomodels.ReservoirGrid.from_csv` for the
    structure.

* ``smoothing_sigmas``: Length-scale in meters used for optional Gaussian smoothing of reservoir grids.

  * ``pressure``: Length-scale for pressure grid smoothing.
  * ``compressibility``: Length-scale for compressibility grid smoothing.
  * ``thickness``: Length-scale for thickness grid smoothing.

* ``median_kernel_sizes``: Kernel size in pixels used for optional median filtering of reservoir grids.

  * ``pressure``: Kernel size for pressure median filtering.
  * ``compressibility``: Kernel size for compressibility median filtering.
  * ``thickness``: Kernel size for thickness median filtering.

Outputs
-------

The following output files are generated (examples for the figures are below each description):

* ``SeismologicalModel_<type>.<label>_<start>-<end>.csv``: Sampled variable values.
* ``posterior_<label>.<format>``: Plot with posterior distribution of sampled variables.

  .. image:: figures/posterior_ets0.htb3.png
     :align: center
     :alt: Posterior distribution for all variables

* ``joint_<label>.<format>``: Plot with joint distribution of sampled variable pairs.

  .. image:: figures/joint_ets0.htb3.png
     :align: center
     :alt: Joint distribution of variable pairs
