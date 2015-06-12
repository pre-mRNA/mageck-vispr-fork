# coding: utf-8
__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import os

import yaml


def postprocess_config(config):
    for key in "library species assembly targets sgrnas experiments samples".split():
        if key not in config:
            raise Exception("Missing key in config file: {}".format(key))
    for sample in config["samples"]:
        if not isinstance(config["samples"][sample], list):
            config["samples"][sample] = [config["samples"][sample]]


def vispr_config(input, output, wildcards, config):
    relpath = lambda path: os.path.relpath(path, "results")
    vispr_config = {
        "experiment": wildcards.experiment,
        "species": config["species"],
        "assembly": config["assembly"],
        "targets": {
            "results": relpath(input.results),
            "genes": config["targets"]["genes"]
        },
        "sgrnas": {
            "counts": relpath(input.counts),
        }
    }
    if "fastqc" in input.keys():
        vispr_config["fastqc"] = {
            sample: [
                relpath(os.path.join(fastqc, data, "fastqc_data.txt"))
                for data in sorted(os.listdir(fastqc))
            ]
            for sample, fastqc in zip(config["samples"], input.fastqc)
        }
    if "mapstats" in input.keys():
        vispr_config["sgrnas"]["mapstats"] = relpath(input.mapstats)
    if "controls" in config["targets"]:
        vispr_config["targets"]["controls"] = config["targets"]["controls"]
    if config["sgrnas"] is not None and "info" in config["sgrnas"]:
        vispr_config["sgrnas"]["info"] = config["sgrnas"]["info"]
    with open(output[0], "w") as f:
        yaml.dump(
            vispr_config,
            f,
            default_flow_style=False
        )
