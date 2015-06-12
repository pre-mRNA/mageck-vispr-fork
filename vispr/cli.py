# coding: utf-8
from __future__ import absolute_import, division, print_function

__author__ = "Johannes Köster"
__copyright__ = "Copyright 2015, Johannes Köster, Liu lab"
__email__ = "koester@jimmy.harvard.edu"
__license__ = "MIT"


import argparse
import logging
import sys
import string
import random
import os
import shutil

import yaml

from vispr import Screens, VisprError, Screen
from vispr.server import app
from vispr.version import __version__


def init_server(*configs):
    app.screens = Screens()
    for path in configs:
        with open(path) as f:
            config = yaml.load(f)
        try:
            app.screens.add(config, parentdir=os.path.dirname(path))
        except KeyError as e:
            raise VisprError(
                "Syntax error in config file {}. Missing key {}.".format(path,
                                                                         e))
    app.secret_key = ''.join(
        random.choice(string.ascii_uppercase + string.digits)
        for _ in range(30))
    print("Starting server.", "",
          "Open:  go to http://127.0.0.1:5000 in your browser.",
          "Close: hit Ctrl-C in this terminal.",
          file=sys.stderr,
          sep="\n")
    app.run()


def init_workflow(directory):
    try:
        os.makedirs(directory)
    except OSError:
        # either directory exists (then we can ignore) or it will fail in the
        # next step.
        pass
    for f in ["Snakefile", "config.yaml", "README.md", "conda.txt"]:
        source = os.path.join(os.path.dirname(__file__), "workflow", f)
        target = os.path.join(directory, f)
        if f in ["Snakefile", "config.yaml"] and os.path.exists(target):
            shutil.copy(target, target + ".old")
        shutil.copy(source, target)


def test_server():
    os.chdir(os.path.join(os.path.dirname(__file__), "tests"))
    init_server("leukemia.yaml", "melanoma.yaml")


def print_example_config():
    print(open(os.path.join(os.path.dirname(__file__),
                            "example.config.yaml")).read())


def plots(configpath, prefix):
    directory = os.path.dirname(prefix)
    if not os.path.exists(directory):
        os.makedirs(directory)

    with open(configpath) as f:
        screen = Screen(yaml.load(f), parentdir=os.path.dirname(configpath))

    def write(json, name):
        with open(prefix + name + ".vega.json", "w") as out:
            print(json, file=out)

    if screen.fastqc is not None:
        write(screen.fastqc.plot_gc_content(), "gc-content")
        write(screen.fastqc.plot_base_quality(), "base-quality")
        write(screen.fastqc.plot_seq_quality(), "seq-quality")
    if screen.mapstats is not None:
        write(screen.mapstats.plot_mapstats(), "mapped-unmapped")
        write(screen.mapstats.plot_zerocounts(), "zerocounts")
        write(screen.mapstats.plot_gini_index(), "gini-index")
    write(screen.rnas.plot_normalization(), "readcounts")
    write(screen.rnas.plot_readcount_cdf(), "readcount-cdf")
    write(screen.rnas.plot_correlation(), "correlation")
    write(screen.rnas.plot_pca(1, 2), "pca-1-2")
    write(screen.rnas.plot_pca(1, 3), "pca-1-3")
    write(screen.rnas.plot_pca(2, 3), "pca-2-3")
    for condition, results in screen.targets.items():
        for selection, results in results.items():
            pre = ".".join(([] if condition == "default" else [condition]) +
                           [selection.replace(" ", "-")])
            write(results.plot_pvals(), pre + ".p-values")
            write(results.plot_pval_hist(), pre + ".p-value-hist")


def main():
    # create arg parser
    parser = argparse.ArgumentParser(
        "An HTML5-based interactive visualization of CRISPR/Cas9 screen data.")
    parser.add_argument("--version",
                        action="store_true",
                        help="Print version info.")
    parser.add_argument("--debug",
                        action="store_true",
                        help="Print debug info.")
    subparsers = parser.add_subparsers(dest="subcommand")

    config = subparsers.add_parser(
        "config",
        help="Print an example VISPR config file. Pipe the output into a file "
        "and edit it to setup a new experiment to be displayed in VISPR.")

    server = subparsers.add_parser("server", help="Start the VISPR server.")
    server.add_argument(
        "config",
        nargs="+",
        help="YAML config files. Each file points to the results of one "
        "MAGeCK run.")

    plot = subparsers.add_parser(
        "plot",
        help="Plot visualizations in VEGA JSON format.")
    plot.add_argument("config",
                      help="YAML config file pointing to MAGeCK results.")
    plot.add_argument("prefix",
                      help="Prefix of all resulting plots. "
                      "This can be a path to a subdirectory.")

    workflow = subparsers.add_parser(
        "init-workflow",
        help="Initialize the MAGeCK/VISPR workflow "
        "in a given directory. This will "
        "install a Snakefile, a README and a "
        "config file in this directory. "
        "Configure the config file according "
        "to your needs, and run the workflow "
        "with Snakemake "
        "(https://bitbucket.org/johanneskoester/snakemake).")
    workflow.add_argument("directory",
                          help="Path to the directory where the "
                          "workflow shall be initialized.")

    subparsers.add_parser("test",
                          help="Start the VISPR server with some included "
                          "test data.")

    args = parser.parse_args()

    logging.basicConfig(format="%(message)s",
                        level=logging.DEBUG if args.debug else logging.INFO,
                        stream=sys.stderr)
    logging.getLogger('werkzeug').setLevel(logging.ERROR)

    try:
        if args.version:
            print(__version__)
            exit(0)
        if args.subcommand == "server":
            init_server(*args.config)
        elif args.subcommand == "test":
            test_server()
        elif args.subcommand == "init-workflow":
            init_workflow(args.directory)
        elif args.subcommand == "config":
            print_example_config()
        elif args.subcommand == "plot":
            plots(args.config, args.prefix)
        else:
            parser.print_help()
            exit(1)
    except VisprError as e:
        logging.error(e)
        exit(1)
    except ImportError as e:
        print("{}. Please ensure that all dependencies from "
              "requirements.txt are installed.".format(e),
              file=sys.stderr)
        exit(1)
    exit(0)
