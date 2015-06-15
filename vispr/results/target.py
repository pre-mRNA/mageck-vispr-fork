# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import json
from itertools import combinations
from operator import itemgetter

from flask import render_template
import pandas as pd
import numpy as np

from vispr.results.common import lru_cache, AbstractResults, templates


class Results(AbstractResults):
    """Keep and display target results."""

    def __init__(self, dataframe):
        """
        Arguments

        dataframe -- path to file containing MAGeCK target (gene) summary. Alternatively, a dataframe.
        controls  -- path to file containing control genes. Alternatively, a dataframe.
        """
        self.df = dataframe
        self.df.sort("p-value", inplace=True)
        self.df.reset_index(drop=True, inplace=True)
        self.df["log10-p-value"] = -np.log10(self.df[["p-value"]])
        self.df["idx"] = self.df.index
        self.df.index = self.df["target"]

    def plot_pvals(self, control_targets=None, mode="hide"):
        """
        Plot the gene ranking in form of their p-values as line plot.

        Arguments
        positive -- if true, plot positive selection scores, else negative selection
        """
        data = self.df[["idx", "log10-p-value"]]

        i5 = self.df["fdr"].searchsorted(0.05)[0]
        i25 = self.df["fdr"].searchsorted(0.25)[0]

        fdr5 = data.iloc[i5]["log10-p-value"]
        fdr25 = data.iloc[i25]["log10-p-value"]
        fdr5label = "{:.0%} FDR".format(self.df.iloc[i5]["fdr"])
        fdr25label = "{:.0%} FDR".format(self.df.iloc[i25]["fdr"])

        top5 = self.df[["idx", "log10-p-value", "target"]]
        if control_targets:
            if mode == "hide":
                valid = self.df["target"].apply(lambda target: target not in
                                                control_targets)
                top5 = top5[valid]
            elif mode == "show-only":
                valid = self.df["target"].apply(lambda target: target in
                                                control_targets)
                top5 = top5[valid]
        top5 = top5.ix[:5]

        plt = templates.get_template("plots/pvals.json").render(
            pvals=data.to_json(orient="records"),
            highlight=top5.to_json(orient="records"),
            fdr5=fdr5,
            fdr25=fdr25,
            fdr5label=fdr5label,
            fdr25label=fdr25label)
        return plt

    def get_pvals_highlight_targets(self, highlight_targets):
        data = self.df[["idx", "log10-p-value", "target"]]
        return data.ix[highlight_targets]

    def plot_pval_hist(self):
        edges = np.arange(0, 1.1, 0.1)
        counts, _ = np.histogram(self.df["p-value"], bins=edges)
        bins = edges[1:]

        hist = pd.DataFrame({"bin": bins, "count": counts})
        return templates.get_template("plots/pval_hist.json").render(
            hist=hist.to_json(orient="records"))

    def ids(self, fdr):
        valid = self.df["fdr"] <= fdr
        return set(self.df.ix[valid, "target"])


def overlap(*targets):
    isect = set(targets[0])
    for other in targets[1:]:
        isect &= other
    return isect


def overlaps(order, **targets):
    """
    Arguments
    order   -- 1: single condition, 2: overlap of 3 conditions, 3: overlap of 3 conditions...
    targets -- labels and targets to compare
    """
    for c in combinations(targets.items(), order):
        isect = overlap(*map(itemgetter(1), c))
        labels = list(map(itemgetter(0), c))
        yield labels, len(isect)


def plot_overlap_chord(**targets):
    ids = {label: i for i, label in enumerate(targets)}
    data = []
    for s in range(2, len(targets) + 1):
        for labels, isect in overlaps(s, **targets):
            data.append([{"group": ids[label],
                          "value": isect} for label in labels])
    for label, t in targets.items():
        excl = set(t)
        for l, t in targets.items():
            if l != label:
                excl -= t
        data.append([{"group": ids[label], "value": len(excl)}])
    return json.dumps({
        "connections": data,
        "labels": {i: label
                   for label, i in ids.items()}
    })


def plot_overlap_venn(**targets):
    data = []
    for s in range(1, len(targets) + 1):
        for labels, isect in overlaps(s, **targets):
            data.append({"sets": labels, "size": isect})
    return json.dumps(data)
