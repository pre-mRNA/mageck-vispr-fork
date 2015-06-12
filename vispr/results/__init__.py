# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os
from itertools import filterfalse

import pandas as pd

from vispr.results import target
from vispr.results import rna
from vispr.results import fastqc
from vispr.results import mapstats


class Screens(object):
    def __init__(self):
        self.screens = {}

    def add(self, config, parentdir="."):
        screen = config["experiment"]
        self.screens[screen] = Screen(config, parentdir=parentdir)

    def __iter__(self):
        return map(self.__getitem__, sorted(self.screens.keys()))

    def __getitem__(self, screen):
        return self.screens[screen]

    def _overlap_targets(self, fdr=0.05, items=None):
        assert items is not None
        def label(screen, condition, selection):
            if condition != "default":
                return screen, condition, selection
            else:
                return screen, selection
        return {
            " | ".join(label(screen, condition, selection)):
            self.screens[screen].targets[condition][selection].ids(fdr)
            for screen, condition, selection in items
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


class Screen(object):
    def __init__(self, config, parentdir="."):
        def get_path(relpath):
            if relpath is None:
                return None
            return os.path.join(parentdir, relpath)

        self.name = config["experiment"]

        self.targets = parse_target_results(
            get_path(config["targets"]["results"]))
        self.is_genes = config["targets"].get("genes", False)
        self.species = config["species"].upper()
        self.assembly = config["assembly"]

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
                sample: map(get_path, paths)
                for sample, paths in config["fastqc"].items()
            })

        self.control_targets = set()
        if "controls" in config["targets"]:
            self.control_targets = set(
                pd.read_table(get_path(config["targets"]["controls"]),
                              header=None,
                              squeeze=True,
                              na_filter=False))

    def targets(self, positive=True):
        return self.pos_targets if positive else self.neg_targets


def parse_target_results(path,
                         selections=["negative selection",
                                     "positive selection"]):
    results = pd.read_table(path, na_filter=False)
    paths = [col.split("|") for col in results.columns]

    max_depth = max(map(len, paths))
    if max_depth == 3:
        # MLE format
        def get_results(selection):
            res = results["id {cond}|beta {cond}|{sel}|p-value {cond}|{sel}|fdr".format(
                cond=condition,
                sel=selection[:3]).split()]
            res.columns = ["target", "beta", "p-value", "fdr"]
            return target.Results(res.copy())

        conditions = [path[0] for path in paths if len(path) > 1]
        targets = {
            condition:
            {selection: get_results(selection)
             for selection in selections}
        }

    elif max_depth == 2:
        # RRA format
        def get_results(selection):
            res = results["id {sel}|score {sel}|p-value {sel}|fdr".format(
                sel=selection[:3]).split()]
            res.columns = ["target", "score", "p-value", "fdr"]
            return target.Results(res.copy())

        targets = {
            "default": {
                selection: get_results(selection)
                for selection in selections
            }
        }
    else:
        raise IOError("Invalid target results format.")
    return targets
