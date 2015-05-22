# coding: utf-8
from __future__ import absolute_import, division, print_function

import os

import pandas as pd

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
            self.screens[screen].targets(positive).ids(fdr)
            for screen, positive in items
        }

    def plot_overlap_chord(self, fdr=0.05, items=None):
        return target.plot_overlap_chord(**self._overlap_targets(fdr=fdr,
                                                                 items=items))

    def plot_overlap_venn(self, fdr=0.05, items=None):
        plt = target.plot_overlap_venn(**self._overlap_targets(fdr=fdr,
                                                               items=items))
        return plt

    def overlap(self, fdr=0.05, items=None):
        return target.overlap(*self._overlap_targets(fdr=fdr,
                                                     items=items).values())


class Screen:
    def __init__(self, config, parentdir="."):
        def get_path(relpath):
            if relpath is None:
                return None
            return os.path.join(parentdir, relpath)

        self.name = config["experiment"]

        targets = get_path(config["targets"]["results"])
        self.pos_targets = target.Results(targets, positive=True)
        self.neg_targets = target.Results(targets, positive=False)
        self.is_genes = config["targets"].get("genes", False)
        if self.is_genes:
            self.species = config["targets"]["species"]

        self.rnas = rna.Results(
            get_path(config["sgrnas"]["counts"]),
            info=get_path(config["sgrnas"].get("info", None)))
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

        self.control_targets = set()
        if "controls" in config["targets"]:
            self.control_targets = set(
                pd.read_table(config["targets"]["controls"],
                              header=None,
                              squeeze=True,
                              na_filter=False))

    def targets(self, positive=True):
        return self.pos_targets if positive else self.neg_targets
