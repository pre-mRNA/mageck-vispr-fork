import argparse
import logging
import sys
import string
import random

import yaml

from vispr import TargetResults, Results, VisprError


def write_plot(plt, path):
    with open(path, "w") as f:
        f.write(plt.to_json())


def main():
    # create arg parser
    parser = argparse.ArgumentParser("An HTML5-based interactive visualization of CRISPR/Cas9 screen data.")
    parser.add_argument("--debug", action="store_true", help="Print debug info.")
    subparsers = parser.add_subparsers(dest="subcommand")

    # server subcommand
    server = subparsers.add_parser("server")
    server.add_argument("config", help="YAML config file containing results.")

    # plot subcommand
    plots = subparsers.add_parser("plot")
    plots.add_argument("feature_results", help="Gene results as provided by MAGeCK.")
    plots.add_argument("--feature-rra-pos", help="Plot RRA scores for positive selection.")
    plots.add_argument("--feature-rra-neg", help="Plot RRA scores for negative selection.")

    args = parser.parse_args()

    logging.basicConfig(
        format=" %(asctime)s: %(message)s",
        level=logging.DEBUG if args.debug else logging.INFO,
        stream=sys.stderr,
    )

    try:

        if args.subcommand == "plot":
            gene_results = FeatureResults(args.gene_results)
            if args.feature_rra_pos:
                write_plot(feature_results.plot_rra(positive=True), args.gene_rra_pos)
            if args.feature_rra_neg:
                write_plot(feature_results.plot_rra(positive=False), args.gene_rra_neg)
        elif args.subcommand == "server":
            from vispr.server import app

            with open(args.config) as f:
                config = yaml.load(f)
            app.results = Results(config)
            app.secret_key = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(30))
            app.run()
    except VisprError as e:
        logging.error(e)
        exit(1)
    exit(0)
