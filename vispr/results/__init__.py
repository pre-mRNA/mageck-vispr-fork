# coding: utf-8
from __future__ import absolute_import, division, print_function

import os

from vispr.results import target
from vispr.results import rna
from vispr.results import fastqc
from vispr.results import mapstats


class Screens:
    def __init__(self):
        self.screens = {}

    def add(self, config, parentdir="."):
        screen = config["experiment"]
        self.screens[screen] = Screen(config, parentdir=parentdir)

    def __iter__(self):
        return iter(sorted(self.screens.keys()))

    def __getitem__(self, screen):
        return self.screens[screen]

    def _overlap_targets(self, fdr=0.05, items=None):
        if items is None:
            items = [(screen, positive)
                     for screen in self.screens for positive in (True, False)]
        selection = ["-", "+"]
        return {
            "{} {}".format(screen, selection[positive]):
            self.screens[screen].targets.targets(fdr,
                                                 positive=positive)
            for screen, positive in items
        }

    def plot_overlap_chord(self, fdr=0.05, items=None):
        return target.plot_overlap_chord(**self._overlap_targets(fdr=fdr,
                                                                 items=items))

    def plot_overlap_venn(self, fdr=0.05, items=None):
        return target.plot_overlap_venn(**self._overlap_targets(fdr=fdr,
                                                                items=items))


class Screen:
    def __init__(self, config, parentdir="."):
        def get_path(relpath):
            if relpath is None:
                return None
            return os.path.join(parentdir, relpath)

        self.name = config["experiment"]

        self.targets = target.Results(
            get_path(config["targets"]["results"]),
            controls=get_path(config["targets"].get("controls", None)))
        self.is_genes = config["targets"].get("genes", False)
        if self.is_genes:
            self.species = config["targets"]["species"]

        self.rnas = rna.Results(get_path(config["sgrnas"]["counts"]), info=get_path(config["sgrnas"].get("info", None)))
        self.mapstats = None
        if "mapstats" in config["sgrnas"]:
            self.mapstats = mapstats.Results(
                get_path(config["sgrnas"]["mapstats"]))

        self.fastqc = None
        if "fastqc" in config:
            self.fastqc = fastqc.Results(**{
                sample: get_path(path)
                for sample, path in config["fastqc"].items()
            })
