# coding: utf-8
from __future__ import absolute_import, division, print_function

import json
try:
    from functools import lru_cache
except ImportError:
    # dummy cache, i.e. Python 2 version will be a bit slower.
    def lru_cache():
        def dummy(func):
            return func
        return dummy

import numpy as np
import pandas as pd
import vincent
from vincent.marks import MarkProperties, MarkRef, Mark
from vincent.transforms import Transform
from vincent.properties import PropertySet
from vincent.values import ValueRef

from vispr.version import __version__


class VisprError(Exception):
    pass


class Results:
    def __init__(self, config):
        self.targets = {}
        self.rnas = {}
        try:
            for screen, files in config.items():
                self.targets[screen] = TargetResults(files["target_results"])
                self.rnas[screen] = RNAResults(files["rna_counts"])
        except KeyError:
            raise VisprError("No results for screen {}".format(screen))

    @property
    def screens(self):
        return self.targets.keys()


class AbstractResults:
    def __init__(self, dataframe):
        """
        Arguments:
        dataframe -- a pandas data frame or its path consisting of per gene results as produced by MAGeCK
        """
        if isinstance(dataframe, str):
            dataframe = pd.read_table(dataframe, na_filter=False)
        self.df = dataframe

    def __getitem__(self, slice):
        return self.df.__getitem__(slice)


class TargetResults(AbstractResults):
    """Keep and display feature results."""

    @lru_cache()
    def get_pvals(self, positive=True):
        # select column and sort
        col = "p.{}".format("pos" if positive else "neg")
        pvals = -np.log10(self.df[[col]])
        data = pd.concat([self.df[["id"]], pvals], axis=1).sort(col, ascending=False).reset_index(drop=True)
        return data, col

    def plot_pvals(self, positive=True):
        """
        Plot the gene ranking in form of their p-values as line plot.

        Arguments
        positive -- if true, plot positive selection scores, else negative selection
        """
        data, col = self.get_pvals(positive=positive)

        # create plot
        plt = vincent.Line(data, columns=[col], width=300, height=200)

        # Circles for each data point
        from_ = MarkRef(
            data='table',
            transform=[Transform(type='facet', keys=['data.idx'])])
        enter_props = PropertySet(
            x=ValueRef(scale='x', field="data.idx"),
            y=ValueRef(scale='y', field="data.val"),
            size=ValueRef(value=100),
            opacity=ValueRef(value=0),
            fill=ValueRef(scale="color", field='data.idx'))
        update_props = PropertySet(
            opacity=ValueRef(value=0)
        )
        hover_props = PropertySet(
            opacity=ValueRef(value=1)
        )
        marks = [Mark(type='symbol',
                      properties=MarkProperties(enter=enter_props, update=update_props, hover=hover_props))]
        mark_group = Mark(type='group', from_=from_, marks=marks)
        plt.marks.append(mark_group)

        # format plot
        plt.axes[0].ticks = 1
        plt.axis_titles(x="Targets", y="-log10 p-value")

        return plt

    def get_pvals_idx(self, target, positive=True):
        data, _ = self.get_pvals(positive=positive)
        idx = data.index.values[(data["id"] == target).values]
        assert len(idx) == 1
        return int(idx[0])


class RNAResults(AbstractResults):
    def by_target(self, target):
        print(target)
        return self.df.loc[self.df["Gene"] == target].ix[:, self.df.columns != "Gene"]
