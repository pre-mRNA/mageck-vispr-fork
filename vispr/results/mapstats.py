import pandas as pd
from flask import render_template

from vispr.results.common import AbstractResults

class Results(AbstractResults):
    def plot_mapstats(self):
        data = self.df[["Label", "Reads", "Mapped"]]
        data.columns = ["label", "reads", "mapped"]
        data.insert(2, "unmapped_percentage", ((data["reads"] - data["mapped"]) / data["reads"]).apply("{:.1%}".format))
        width = 20 * data.shape[0]
        return render_template("plots/mapstats.json",
                               data=data.to_json(orient="records"),
                               width=width)
