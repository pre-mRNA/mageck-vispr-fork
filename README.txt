VISPR - A visualization server for CRISPR data.
===============================================

There are two ways to test VISPR.

Testing VISPR with the Miniconda Python distribution
----------------------------------------------------

Install Miniconda for Python 3 from here:

http://conda.pydata.org/miniconda.html

This will install a minimal Python 3, together with the conda
package manager.

Then, open a terminal in this directory and execute

    conda install --file conda.txt

to install all required packages.
Then, issue

    python setup.py nosetests

to run a test instance of VISPR.


Testing VISPR with another Python distribution (for experts)
------------------------------------------------------------

Make sure that you have numpy and nose installed. Then, you can
issue

    python setup.py nosetests

which will install remaining dependencies locally.
