# coding: utf-8
__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


from __future__ import absolute_import, division, print_function

from vispr.version import __version__

from vispr.results import Screens, Screen


class VisprError(Exception):
    pass
