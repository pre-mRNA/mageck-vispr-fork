import json
from itertools import combinations
from operator import itemgetter

from flask import render_template
import pandas as pd
import numpy as np

from vispr.results.common import lru_cache, AbstractResults


class Results(AbstractResults):
    """Keep and display feature results."""

    @lru_cache()
    def get_pvals(self, positive=True):
        # select column and sort
        col = "pos" if positive else "neg"
        pvals = -np.log10(self.df[["p." + col]])
        fdr = self.df[["fdr." + col]]
        data = pd.concat([self.df[["id"]], pvals, fdr],
                         axis=1).sort("p." + col,
                                      ascending=False).reset_index(drop=True)
        return data

    def plot_pvals(self, positive=True):
        """
        Plot the gene ranking in form of their p-values as line plot.

        Arguments
        positive -- if true, plot positive selection scores, else negative selection
        """
        data = self.get_pvals(positive=positive)

        pvals = pd.DataFrame(
            {"idx": data.index,
             "pval": data.ix[:, 1],
             "fdr": data.ix[:, 2]})

        return render_template("plots/pvals.json",
                               pvals=pvals.to_json(orient="records"))

    @lru_cache()
    def get_pvals_highlight(self, positive=True):
        pvals = self.get_pvals(positive=positive)
        pvals = pd.DataFrame({
            "idx": pvals.index,
            "pval": pvals.ix[:, 1],
            "label": pvals["id"]
        })
        pvals.index = pvals["label"]
        return pvals

    def get_pvals_highlight_targets(self, highlight_targets, positive=True):
        pvals = self.get_pvals_highlight(positive=positive)
        pvals = pvals.ix[highlight_targets]

        return pvals

    def plot_pval_hist(self, positive=True):
        data = self.get_pvals(positive=positive)
        edges = np.arange(0, 1.1, 0.1)
        counts, _ = np.histogram(data.ix[:, 1], bins=edges)
        bins = edges[1:]

        hist = pd.DataFrame({"bin": bins, "count": counts})
        return render_template("plots/pval_hist.json",
                               hist=hist.to_json(orient="records"))

    def get_pvals_idx(self, target, positive=True):
        data, _ = self.get_pvals(positive=positive)
        idx = data.index.values[(data["id"] == target).values]
        assert len(idx) == 1
        return int(idx[0])

    def targets(self, fdr, positive=True):
        col = "pos" if positive else "neg"
        valid = self.df["fdr." + col] <= fdr
        return set(self.df.ix[valid, "id"])


def plot_overlap_chord(**targets):
    ids = {label: i for i, label in enumerate(targets)}
    data = []
    for s in range(2, len(targets) + 1):
        for c in combinations(targets.items(), s):
            isect = set(c[0][1])
            for other in map(itemgetter(1), c[1:]):
                isect &= other
            data.append([{"group": ids[label], "value": len(isect)} for label in map(itemgetter(0), c)])
    for label, t in targets.items():
        excl = set(t)
        for l, t in targets.items():
            if l != label:
                excl -= t
        data.append([{"group": ids[label], "value": len(excl)}])

    return json.dumps({"connections": data, "labels": {i: label for label, i in ids.items()}})
