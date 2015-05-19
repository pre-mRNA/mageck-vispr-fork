VISPR - A visualization server for CRISPR data.
===============================================

There are two ways to test VISPR.

Installing VISPR with the Miniconda Python distribution
-------------------------------------------------------

Install Miniconda for Python 3 from here:

http://conda.pydata.org/miniconda.html

This will install a minimal Python 3, together with the conda
package manager.

Then, open a terminal in this directory and execute

    conda install --file requirements.txt

to install all required packages.
Install VISPR with

    python setup.py install

See below for running a test instance of VISPR.

Installing VISPR with another Python distribution (for experts)
---------------------------------------------------------------

Make sure that you have numpy and pandas installed. Else, their
automatic compilation with the command below would take very long.
Then, you can issue

    python setup.py install

or

    python setup.py install --user

if you want to install VISPR to your home folder.
All remaining dependencies will be installed automatically.

Running VISPR
-------------

After successful installation, you can run VISPR by executing the command

    vispr server path/to/config.yaml

See below for the config file format.
If you only want to test VISPR, you can run a test instance with example
data by executing

    vispr test


Config file format
------------------

The config file has to be in YAML format. It has to contain a section for
each experiment. Each experiment contains two entries, one pointing to the
target (or gene) result file provided by "MAGeCK test", one pointing to a
file with normalized read counts as provided by "MAGeCK count":

    myexperiment:
        target_results: gene_summary.txt
        rna_counts: normalized_counts.txt

