import json

import pandas as pd
import vincent

from vispr.version import __version__


class UserError(Exception):
    pass


class GeneResults:
    """Keep and display gene results."""

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
        plt.axis_titles(x="Genes", y="RRA score ({} selection)".format("positive" if positive else "negative"))
        plt.colors(brew='Set1')

        return plt
