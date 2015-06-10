VISPR - A visualization server for CRISPR data.
===============================================

VISPR is a visualization framework for CRISPR/Cas9 knockout screen experiments.
To install VISPR, we recommend the Miniconda Python distribution (see below).
Currently, VISPR is under heavy development and not yet recommended for public
use.


Installing VISPR with the Miniconda Python distribution
-------------------------------------------------------

Install Miniconda for Python 3 from here:

http://conda.pydata.org/miniconda.html

This will install a minimal Python 3, together with the conda
package manager (if preferred, you can also use Python 2).

Then, open a terminal and execute

    conda install -c johanneskoester vispr

to install VISPR with all dependencies.
See below for running a test instance of VISPR.


Running VISPR
-------------

After successful installation, you can run VISPR by executing the command

    vispr server path/to/config.yaml

See below for the config file format.
If you only want to test VISPR, you can run a test instance with example
data by executing

    vispr test


Installing VISPR with another Python distribution (for experts)
---------------------------------------------------------------

Make sure that you have numpy, scipy, scikit-learn and pandas installed.
Else, their automatic compilation with the command below would take very long.
Then, you can issue

    pip install vispr

or

    pip install vispr --user

if you want to install VISPR without admin rights.
All remaining dependencies will be installed automatically.
