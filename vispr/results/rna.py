# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import json
from functools import partial

from flask import render_template
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import complete, leaves_list
from scipy.spatial.distance import pdist, squareform

from vispr.results.common import lru_cache, AbstractResults, templates


class Results(AbstractResults):
    def __init__(self, dataframe, info=None):
        super(Results, self).__init__(dataframe)
        # reorder columns lexicographically
        columns = list(self.df.columns[:2]) + sorted(self.df.columns[2:])
        self.df = self.df[columns]
        self.df.columns = ["rna", "target"] + list(self.samples)
        self.df.index = self.df["target"]
        self.info = None
        if info is not None:
            self.info = pd.read_table(
                info,
                na_filter=False,
                index_col=3,
                names=["chrom", "start", "stop", "score"])

    @property
    def samples(self):
        return list(self.df.columns[2:])

    def by_target(self, target):
        first_sample = self.df.columns[2]
        data = self.df.ix[target, self.df.columns != "target"]

        data.sort(first_sample, inplace=True)
        if self.info is not None:
            data.index = data["rna"]
            info = self.info.ix[data["rna"]]
            if not info["score"].hasnans() and not info["start"].hasnans():
                data.insert(1, "efficiency", info["score"])
                data["chrom pos"] = info["start"]
                data.sort("efficiency", inplace=True)
        return data

    def track(self):
        data = self.df.ix[:, self.df.columns != "target"]
        data["desc"] = [
            "na|@chr{}:{}-{}|".format(*row) for _, row in
            self.info.ix[data["rna"], ["chrom", "start", "stop"]].iterrows()
        ]
        return "#1.2\n{rnas}\t{samples}\n{data}".format(
            rnas=data.shape[0],
            samples=data.shape[1] - 2,
            data=data.to_csv(sep="\t",
                             index=False,
                             columns=["rna", "desc"] + self.samples,
                             header=["Name", "Description"] + self.samples))

    def target_locus(self, target):
        loci = self.info.ix[self.df.ix[target, "rna"], ["chrom", "start",
                                                        "stop"]]
        return loci["chrom"][0], loci["start"].min(), loci["stop"].max()

    def plot_normalization(self):
        _, leaves, _ = self.clustering()

        counts = np.log10(self.counts() + 1)
        # sort columns by clustering results
        counts = counts.ix[:, leaves]
        data = pd.DataFrame({
            "label": counts.columns,
            "median": counts.median(),
            "lo": counts.quantile(0.25),
            "hi": counts.quantile(0.75),
            "q10": counts.quantile(0.1),
            "q90": counts.quantile(0.9),
        })
        width = 20 * counts.shape[1]
        return templates.get_template("plots/normalization.json").render(
                               data=data.to_json(orient="records"),
                               width=width)

    def plot_readcount_cdf(self):
        counts = np.log10(self.counts() + 1)

        data = []
        for sample in counts.columns:
            d = counts[sample].value_counts(normalize=True, sort=False, bins=100).cumsum()
            data.append(pd.DataFrame({
                "density": d,
                "count": d.index,
                "sample": sample
            }))
        data = pd.concat(data)
        return templates.get_template("plots/readcounts.json").render(
                               data=data.to_json(orient="records"))


    def counts(self):
        return self.df.ix[:, 2:]

    @lru_cache()
    def pca(self, n_components=3):
        counts = np.log10(self.counts().transpose() + 1)
        pca = PCA(n_components=3)
        data = pd.DataFrame(pca.fit_transform(counts))
        min_coeff = data.min().min()
        max_coeff = data.max().max()
        fields = ["PC{} ({:.0%})".format(i + 1, expl_var)
                  for i, expl_var in enumerate(pca.explained_variance_ratio_)]
        data.columns = fields
        data["sample"] = counts.index
        return data

    def plot_pca(self, comp_x=1, comp_y=2, legend=True):
        pca = self.pca(3)
        data = pd.DataFrame({
            "label": pca["sample"],
            "x": pca[pca.columns[comp_x - 1]],
            "y": pca[pca.columns[comp_y - 1]],
        })

        plt = templates.get_template("plots/pca.json").render(
                              data=data.to_json(orient="records"),
                              xlabel=pca.columns[comp_x - 1],
                              ylabel=pca.columns[comp_y - 1],
                              legend=legend)
        return plt

    @lru_cache()
    def clustering(self):
        counts = np.log10(self.counts().transpose() + 1)
        # calculate correlation
        corr = 1 - pdist(counts, 'correlation')
        # calculate distance from absolute correlation
        dist = 1 - np.abs(corr)
        # cluster
        clustering = complete(dist)
        leaves = leaves_list(clustering)

        # calculate correlation matrix
        corr = np.round(squareform(corr), 2)
        # fix diagonal, that will contain zeros because squareform expects a dist
        np.fill_diagonal(corr, 1)
        labels = counts.index.values
        return corr, leaves, labels

    def plot_correlation(self):
        corr, leaves, labels = self.clustering()

        # convert to json records
        data = [{"a": labels[i],
                 "b": labels[j],
                 "value": corr[i, j]} for i in leaves for j in leaves]

        size = max(min(50 * len(labels), 700), 300)

        plt = templates.get_template("plots/correlation.json").render(
                              data=json.dumps(data),
                              size=size)
        return plt
