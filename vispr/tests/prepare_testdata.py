import sys
import bz2

import yaml


def compress(f, out):
    with open(f, "rb") as f, bz2.open(out, "wb") as out:
        out.write(f.read())


def url(filename):
    return "https://bitbucket.org/liulab/vispr/downloads/{}".format(filename)


if __name__ == "__main__":
    name = sys.argv[1]
    config = sys.argv[2]
    newconfig = sys.argv[3]
    with open(config) as config:
        config = yaml.load(config)

    if "fastqc" in config:
        new = {}
        for sample, fastqs in config["fastqc"].items():
            new[sample] = []
            for i, f in enumerate(fastqs):
                out = "{}.{}.{}.fastqc_data.txt.bz2".format(name, sample, i)
                compress(f, out)
                new[sample].append(url(out))
        config["fastqc"] = new

    out = "{}.count.normalized.txt.bz2".format(name)
    compress(config["sgrnas"]["counts"], out)
    config["sgrnas"]["counts"] = url(out)

    if "mapstats" in config["sgrnas"]:
        out = "{}.countsummary.txt.bz2".format(name)
        compress(config["sgrnas"]["mapstats"], out)
        config["sgrnas"]["mapstats"] = url(out)

    if "annotation" in config["sgrnas"]:
        out = "{}.sgnra_annotation.bed.bz2".format(name)
        compress(config["sgrnas"]["annotation"], out)
        config["sgrnas"]["annotation"] = url(out)

    out = "{}.gene_summary.txt.bz2".format(name)
    compress(config["targets"]["results"], out)
    config["targets"]["results"] = url(out)

    if "controls" in config["targets"]:
        out = "{}.controls.txt.bz2".format(name)
        compress(config["targets"]["controls"], out)
        config["targets"]["controls"] = url(out)

    with open(newconfig, "w") as newconfig:
        yaml.dump(config, newconfig, default_flow_style=False)
