VISPR - A visualization server for CRISPR data.
===============================================

VISPR is a visualization framework for CRISPR/Cas9 knockout screen experiments.
To install VISPR, we recommend the Miniconda Python distribution (see below).
Currently, VISPR is under heavy development and not yet recommended for public
use.

Installation
------------

The easiest way to install VISPR is to use the Miniconda Python distribution (http://conda.pydata.org/miniconda.html). 
We recommend to use the Python 3.4 variant.
With miniconda installed, you can issue

    conda install -c johanneskoester vispr

to install VISPR with all dependencies.
To update vispr, you can issue

    conda update -c johanneskoester vispr

Usage
-----

After successful installation, you can run VISPR by executing the command

    vispr server path/to/config.yaml

This starts a server process. The VISPR user interface will be rendered in a webbrowser.
See below for the config file format.
If you only want to test VISPR, you can run a test instance with example
data by executing

    vispr test


Configuring VISPR
-----------------

VISPR takes MAGeCK and FastQC results as input. These are provided along with additional parameters as config files. One config file defines one set of results (i.e. one experiment).
VISPR can be invoked with multiple config files (i.e. multiple experiments), allowing to select and compare experiments via the user interface.
The command

    vispr config

prints a config file template to the terminal.


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

Author
------

Johannes Köster (koester@jimmy.harvard.edu),
Liu lab,
Department of Biostatistics and Computational Biology,
Dana-Farber Cancer Institute

License
-------

Licensed under the MIT license (http://opensource.org/licenses/MIT). This project may not be copied, modified, or distributed except according to those terms.