.. _models:

Model input JSONs
=================

All scripts here require a seismological model, which is defined using a JSON file. Most common models
are already distributed here, inside the folder ``inputs/models``.

The structure of these JSON files is as follows:

.. code-block:: json

   {
       "label": "label",
       "type": "ModelName",
       "dof": {
           "variable_name": {
               "distr": "distribution_name",
               "value": 0,
               "lower": 0,
               "upper": 1,
               "mean": 0
           }
       }
   }

A brief summary of the parameters:

* ``label``: A label for the model, used for setting output file names.
* ``type``: Name of the model type used. Must be one of the class names in :mod:`~seismod.seismodels`.
* ``dof``: Parameters for all the variables required by the model (with the exception of ``Mmin``). Each element
  must have a variable name with the following parameters:

  * ``distr``: Distribution used. One of ``constant`` (fixed value), ``uniform`` (uniform distribution),
    ``uniform-log`` (distribution uniform in log-space) or ``exponential`` (bounded exponential distribution).
  * ``value``: Only required for simulations (as for the comparison script) or ``constant`` distributions.
  * ``lower``: Lower limit of the distribution. Required by ``uniform``, ``uniform-log`` and ``exponential``.
  * ``upper``: Upper limit of the distribution. Required by ``uniform`` and ``uniform-log``.
  * ``mean``: Mean of the distribution. Only required by ``exponential``.
