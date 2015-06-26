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
        self.df["log10-p-value"] = -np.log10(self.df["p-value"])
        self.df["idx"] = self.df.index
        self.df.index = self.df["target"]

        pval_cdf = self.df.replace([np.inf, -np.inf], np.nan).dropna()["log10-p-value"].value_counts(normalize=True, sort=False, bins=1000).cumsum()
        pval_cdf.index = np.maximum(0, pval_cdf.index)
        self.pval_cdf = pd.DataFrame({"p-value": pval_cdf.index, "cdf": pval_cdf})

    def get_pval_cdf_points(self, pvals):
        idx = self.pval_cdf["p-value"].searchsorted(pvals, side="right") - 1
        d = self.pval_cdf.iloc[idx]
        d.reset_index(drop=True, inplace=True)
        return d

    def plot_pvals(self, control_targets=None, mode="hide"):
        """
        Plot the gene ranking in form of their p-values as CDF plot.
        """

        fdr_idx = self.df["fdr"].searchsorted([0.05, 0.25], side="right")

        #import pdb; pdb.set_trace()
        fdrs = pd.DataFrame(self.get_pval_cdf_points(self.df.iloc[fdr_idx]["log10-p-value"]))
        fdrs.loc[:, "label"] = self.df.iloc[fdr_idx]["fdr"].apply("{:.0%} FDR".format).values

        top5 = self.df.index
        if control_targets:
            if mode == "hide":
                valid = self.df["target"].apply(lambda target: target not in
                                                control_targets)
                top5 = top5[valid]
            elif mode == "show-only":
                valid = self.df["target"].apply(lambda target: target in
                                                control_targets)
                top5 = top5[valid]
        top5targets = top5[:5]
        top5 = pd.DataFrame(self.get_pval_cdf_points(self.df.ix[top5targets, "log10-p-value"]))
        top5.loc[:, "target"] = top5targets.values

        plt = templates.get_template("plots/pvals.json").render(
            pvals=self.pval_cdf.to_json(orient="records"),
            highlight=top5.to_json(orient="records"),
            fdrs=fdrs.to_json(orient="records"))
        return plt

    def get_pvals_highlight_targets(self, highlight_targets):
        data = pd.DataFrame(self.get_pval_cdf_points(self.df.ix[highlight_targets, "log10-p-value"]))
        data.loc[:, "target"] = highlight_targets
        return data

    def plot_pval_hist(self):
        edges = np.arange(0, 1.01, 0.05)
        counts, _ = np.histogram(self.df["p-value"], bins=edges)
        bins = edges[1:]

        hist = pd.DataFrame({"bin": bins, "count": counts})
        return templates.get_template("plots/pval_hist.json").render(
            hist=hist.to_json(orient="records"))

    def plot_pval_qq(self):
        data = self.df[["log10-p-value"]]
        data["nulldist"] = -np.log10(np.arange(1, len(self) + 1) / len(self))
        return templates.get_template("plots/pval_qq.json").render(
            data=data.to_json(orient="records"))

    def ids(self, fdr):
        valid = self.df["fdr"] <= fdr
        return set(self.df.ix[valid, "target"])

    def __len__(self):
        return self.df.shape[0]


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
