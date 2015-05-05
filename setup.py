# coding: utf-8
from __future__ import absolute_import, division, print_function

import sys

try:
    from setuptools import setup
except ImportError:
    print("Could not load setuptools. Please install the setuptools package.", file=sys.stderr)


# load version info
exec(open("vispr/version.py").read())


setup(
    name="vispr",
    version=__version__,
    author="Johannes Köster",
    author_email="johannes.koester@tu-dortmund.de",
    description="Interactive HTML5 visualization for GeCKO screens.",
    license="MIT",
    url="",
    packages=["vispr", "vispr.tests"],
    include_package_data=True,
    zip_safe=False,
    install_requires=["pandas", "flask", "vincent", "pyyaml", "numpy", "nose"],
    entry_points={"console_scripts": ["vispr = vispr.cli:main"]},
    classifiers=[
        "Development Status :: 4 - Beta",
        #"Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 2",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
