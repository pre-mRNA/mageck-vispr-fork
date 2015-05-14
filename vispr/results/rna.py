import json

from flask import render_template
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import complete, leaves_list
from scipy.spatial.distance import pdist, squareform

from vispr.results.common import lru_cache, AbstractResults


class RNAResults(AbstractResults):
    def by_target(self, target):
        first_sample = self.df.columns[2]
        return self.df.loc[self.df["Gene"] == target].ix[:, self.df.columns !=
                                                         "Gene"].sort(
                                                             first_sample)

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
        return render_template("plots/normalization.json",
                               data=data.to_json(orient="records"),
                               width=width)

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

        plt = render_template("plots/pca.json",
                              data=data.to_json(orient="records"),
                              xlabel=pca.columns[comp_x - 1],
                              ylabel=pca.columns[comp_y - 1],
                              legend=legend, )
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

        plt = render_template("plots/correlation.json", data=json.dumps(data), size=size)
        return plt
