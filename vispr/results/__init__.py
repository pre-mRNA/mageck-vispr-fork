# coding: utf-8
from __future__ import absolute_import, division, print_function

import os

from vispr.results.target import TargetResults
from vispr.results.rna import RNAResults


class Screens:
    def __init__(self):
        self.screens = {}

    def add(self, config, parentdir="."):
        screen = config["experiment"]
        self.screens[screen] = Screen(config, parentdir=parentdir)

    def __iter__(self):
        return iter(self.screens.keys())

    def __getitem__(self, screen):
        return self.screens[screen]


class Screen:
    def __init__(self, config, parentdir="."):
        def get_path(key):
            return os.path.join(parentdir, config[key])

        self.name = config["experiment"]
        self.targets = TargetResults(get_path("target_results"))
        self.rnas = RNAResults(get_path("rna_counts"))
        is_genes = config.get("genes", False)
        self.is_genes = is_genes
        if is_genes:
            self.species = config["species"]
