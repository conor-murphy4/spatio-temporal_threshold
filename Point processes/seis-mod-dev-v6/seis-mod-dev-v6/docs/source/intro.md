# Introduction

This package contains the probability distributions and models necessary to sample variables for the seismological model 
and evaluate the likelihood out-of-sample. All distributions and models are built on top of the [PyMC3](https://docs.pymc.io)
library. 

A console script will be added at installation (`ologp-compare`), which is used to run a comparison 
between different models, whose settings must be specified with a JSON file. 

The code has been implemented for running under Python 3.7 or newer (it has been tested under Python 3.7.2).

## Setting environment

Due to the specific versions required, it is recommended to install this package within a dedicated environment.

### Using a virtual environment (Linux)

To create a virtual environment using the Python module `venv`, run:
```bash
$ python3 -m venv <venv>
```
where `<venv>` is the name of the directory where the environment will reside. This command will only create the 
virtual environment, it still needs to be activated using,
```bash
$ source <venv>/bin/activate
```
Note that the environment needs to be activated independently for every terminal where the code is intended to run. In
case there is a need to deactivate the environment, simply issue the command `deactivate`.

Once the environment is activated, the requirements can be installed using,
```bash
$ pip install -r requirements.txt
```

### Using conda (Linux/Windows)

> **Note:** Installing and running the scripts in Windows should be done using the Anaconda Power Shell.

To create a dedicated environment using `conda` with the requirements installed, run:
```bash
$ conda env create -f environment.yml
```
This only creates a new environment (named `seismod-env`), which need to be activated using,
```bash
$ conda activate seismod-env
```
Finally, to deactivate the environment, run:
```bash
$ conda deactivate
```

## Install

To install the Python library `seismod` and its main scripts, make sure that the environment is active and run:
```bash
$ pip install .
```
from the repository folder (the same that contains `requirements.txt` and `environment.yml` among others).

Alternatively, to install the library in development mode, run:
```bash
$ pip install -e .
```
This allows for modifications made to the source code to affect the installed version without the need to reinstall.

## Testing

To make sure that all the requirements and the library itself are working as expected, it is recommended to run the 
full test suite after installation. To do so, run:
```bash
$ python -m unittest
```
from the repository folder. If no tests have failed, the code is ready for running.

## Building documentation

In order to build the documentation, go to the `docs` folder,
```bash
$ cd docs
```

For Linux, run:
```bash
$ make html
```

For Windows, run:
```bash
> make.bat html
```

The above commands will generate the html documentation inside `build/html` (to start browsing, open `index.html`).
