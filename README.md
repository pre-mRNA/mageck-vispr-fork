# Testing

There are two ways to test VISPR.

## I have a working Python installation

Open a terminal in this directory and issue

    python setup.py nosetests

If you can choose, use Python 3 for VISPR, because it will be faster.

## I don't have a working Python installation

Install Miniconda for Python 3 from here:
http://conda.pydata.org/miniconda.html
This will install a minimal Python 3, together with the conda
package manager.

Then, open a terminal in this directory and execute

    conda install --file conda.txt

to install all required packages.
Finally, follow the procedure described under
"I have a working Python installation".
