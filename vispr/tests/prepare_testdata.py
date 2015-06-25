import sys


def compress(f, out):
    with open(f, "rb") as f, bz2.open(out, "wb") as out:
        out.write(f.read())


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
                out = "{}.{}.{}.fastq_data.txt".format(name, sample, i)
                compress(f, out)
                new[sample].append(out)
        config["fastqc"] = new

    out = "{}.count.normalized.txt".format(name)
    compress(config["sgrnas"]["counts"], out)
    config["sgrnas"]["counts"] = out

    if "mapstats" in config["sgrnas"]:
        out = "{}.countsummary.txt".format(name)
        compress(config["sgrnas"]["mapstats"], out)
        config["sgrnas"]["mapstats"] = out

    if "annotation" in config["sgrnas"]:
        out = "{}.sgnra_info.bed".format(name)
        compress(config["sgrnas"]["annotation"], out)
        config["sgrnas"]["annotation"] = out

    out = "{}.gene_summary.txt".format(name)
    compress(config["targets"]["results"], out)
    config["targets"]["results"] = out

    if "controls" in config["targets"]:
        out = "{}.controls.txt".format(name)
        compress(config["targets"]["controls"], out)
        config["targets"]["controls"] = out

    with open(newconfig, "w") as newconfig:
        yaml.dump(config, newconfig)
