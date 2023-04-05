from setuptools import setup, find_packages

import seismod

setup(
    name='seismod',
    version=seismod.__version__,
    packages=find_packages(exclude=['tests']),
    entry_points={
        "console_scripts": [
            "ologp-compare = seismod.runners:_main_compare",
            "seismod-sample = seismod.runners:_main_sample"
        ]
    },
    url='',
    license='',
    author='NAM PSRA team',
    author_email='',
    description='Tool for sampling seismological model parameters and running comparisons between models'
)
