import numpy as np
import pandas as pd
from flask import render_template

from vispr.results.common import AbstractResults

class Results(AbstractResults):
    def plot_mapstats(self):
        data = pd.DataFrame({
            "label": self.df["Label"],
            "reads": self.df["Reads"],
            "mapped": self.df["Mapped"],
            "unmapped_percentage": ((self.df["Reads"] - self.df["Mapped"]) / self.df["Reads"]).apply("{:.1%}".format)
        })
        width = 20 * data.shape[0]
        return render_template("plots/mapstats.json",
                               data=data.to_json(orient="records"),
                               width=width)

    def plot_zerocounts(self):
        data = pd.DataFrame({
            "label": self.df["Label"],
            "zerocounts": np.log10(self.df["Zerocounts"])
        })
        width = 20 * data.shape[0]
        return render_template("plots/zerocounts.json",
                               data=data.to_json(orient="records"),
                               width=width)
