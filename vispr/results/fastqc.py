import json
import re
from io import StringIO

import pandas as pd
from flask import render_template


class Results:
    def __init__(self, **fastqc_data):
        self.gc_content = []
        self.base_quality = []
        self.seq_quality = []
        for sample, path in fastqc_data.items():
            for name, data in parse_fastqc_data(path):
                if name == "Per sequence GC content":
                    data.columns = ["gc", "density"]
                    data["sample"] = sample
                    data["density"] /= data["density"].sum()
                    self.gc_content.append(data)
                elif name == "Per base sequence quality":
                    data.columns = ["base", "mean", "median", "lo", "hi", "q10", "q90"]
                    data["sample"] = sample
                    self.base_quality.append(data)
                elif name == "Per sequence quality scores":
                    data.columns = ["qual", "density"]
                    data["density"] /= data["density"].sum()
                    data["sample"] = sample
                    self.seq_quality.append(data)
        self.gc_content = pd.concat(self.gc_content)
        self.base_quality = pd.concat(self.base_quality)
        self.seq_quality = pd.concat(self.seq_quality)

    def plot_gc_content(self):
        plt = render_template("plots/gc_content.json", data=self.gc_content.to_json(orient="records"))
        return plt

    def plot_base_quality(self):
        plt = render_template(
            "plots/base_quality.json",
            data=self.base_quality.to_json(orient="records"))
        return plt

    def plot_seq_quality(self):
        plt = render_template(
            "plots/seq_quality.json",
            data=self.seq_quality.to_json(orient="records"))
        return plt

def parse_fastqc_data(path, pattern=re.compile(r"\n\>\>(?P<name>[\w ]+)\t(pass|fail)\n\#(?P<data>.+?)\n\>\>END_MODULE", flags=re.MULTILINE | re.DOTALL)):
    with open(path) as f:
        f = f.read()
        for match in pattern.finditer(f):
            data = pd.read_table(StringIO(match.group("data")))
            yield match.group("name"), data
