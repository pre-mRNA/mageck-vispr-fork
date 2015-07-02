import sys
import os
import tarfile
from io import BytesIO

import yaml


def archive(config, out):
    parentdir = os.path.dirname(config)

    def get_path(relpath):
        if relpath.startswith("http"):
            return relpath
        return os.path.join(parentdir, relpath)

    with open(config) as config:
        config = yaml.load(config)

    def get_out(name):
        return os.path.join(config["experiment"], name)

    mode = None
    if out.endswith(".tar"):
        mode = "w"
    elif out.endswith(".tar.gz"):
        mode = "w:gz"
    elif out.endswith(".tar.bz2"):
        mode = "w:bz2"

    with tarfile.open(out, mode) as tar:
        if "fastqc" in config:
            new = {}
            for sample, fastqs in config["fastqc"].items():
                new[sample] = []
                for i, f in enumerate(fastqs):
                    out = "{}.{}.fastqc_data.txt".format(sample, i)
                    tar.add(get_path(f), out)
                    new[sample].append(out)
            config["fastqc"] = new

        out = get_out("all.count.normalized.txt")
        tar.add(get_path(config["sgrnas"]["counts"]), out)
        config["sgrnas"]["counts"] = out

        if "mapstats" in config["sgrnas"]:
            out = get_out("all.countsummary.txt")
            tar.add(get_path(config["sgrnas"]["mapstats"]), out)
            config["sgrnas"]["mapstats"] = out

        if "annotation" in config["sgrnas"]:
            out = get_out("sgnra_annotation.bed")
            tar.add(get_path(config["sgrnas"]["annotation"]), out)
            config["sgrnas"]["annotation"] = out

        out = get_out("all.gene_summary.txt")
        tar.add(get_path(config["targets"]["results"]), out)
        config["targets"]["results"] = out

        if "controls" in config["targets"]:
            out = get_out("all.controls.txt")
            tar.add(get_path(config["targets"]["controls"]), out)
            config["targets"]["controls"] = out

        newconfig = yaml.dump(config, default_flow_style=False)
        tarinfo = tarfile.TarInfo(get_out("vispr.yaml"))
        tarinfo.size = len(newconfig)
        tar.addfile(tarinfo, fileobj=BytesIO(newconfig.encode()))
