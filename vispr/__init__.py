import json

import pandas as pd
import vincent

from vispr.version import __version__


class VisprError(Exception):
    pass


class Results:
    def __init__(self, config):
        self.targets = {}
        self.rnas = {}
        try:
            for screen, files in config.items():
                self.targets[screen] = TargetResults(files["feature_results"])
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
            dataframe = pd.read_table(dataframe)
        self.df = dataframe

    def __getitem__(self, slice):
        return self.df.__getitem__(slice)


class TargetResults(AbstractResults):
    """Keep and display feature results."""

    def plot_rra(self, positive=True):
        """
        Plot the gene ranking in form of their p-values as line plot.

        Arguments
        positive -- if true, plot positive selection scores, else negative selection
        """
        # select column and sort
        col = "lo.{}".format("pos" if positive else "neg")
        rra = self.df[["id", col]].sort(col, ascending=False).reset_index(drop=True)

        # create plot
        plt = vincent.Line(rra, columns=[col], width=300, height=200)

        print(plt.to_json())

        # format plot
        plt.axes[0].ticks = 1
        plt.axis_titles(x="Targets", y="RRA-score")
        plt.colors(brew='Set1')

        return plt


class RNAResults(AbstractResults):
    def by_target(self, target):
        return self.df.loc[self.df["Gene"] == target].ix[:, self.df.columns != "Gene"]
