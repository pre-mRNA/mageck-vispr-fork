# coding: utf-8

import sys

from setuptools import setup


if sys.version_info < (3,3):
    sys.stdout.write("At least Python 3.3 is required.\n")
    sys.exit(1)


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
    packages=["vispr"],
    zip_safe=True,
    install_requires=["pandas", "flask", "vincent"],
    entry_points={"console_scripts": ["vispr = vispr.cli:main"]},
    package_data={'': ['*.html']},
    classifiers=[
        "Development Status :: 4 - Beta",
        #"Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ]
)
