# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import numpy as np
import pandas as pd
from flask import render_template

from vispr.results.common import AbstractResults, templates

class Results(AbstractResults):
    def __init__(self, dataframe):
        super(Results, self).__init__(dataframe)
        self.df.columns = self.df.columns.str.lower()

    def plot_mapstats(self):
        data = self.df.loc[:, ["label", "reads", "mapped"]]
        data["unmapped_percentage"] = ((self.df["reads"] - self.df["mapped"]) / self.df["reads"]).apply("{:.1%}".format)

        width = 20 * data.shape[0]
        return templates.get_template("plots/mapstats.json").render(
                               data=data.to_json(orient="records"),
                               width=width)

    def plot_zerocounts(self):
        data = self.df.loc[:, ["label"]]
        data["zerocounts"] = np.log10(self.df["zerocounts"])

        width = 20 * data.shape[0]
        return templates.get_template("plots/zerocounts.json").render(
                               data=data.to_json(orient="records"),
                               width=width)

    def plot_gini_index(self):
        data = self.df[["label", "giniindex"]]

        width = 20 * data.shape[0]
        return templates.get_template("plots/gini_index.json").render(
                               data=data.to_json(orient="records"),
                               width=width)
