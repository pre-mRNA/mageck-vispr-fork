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
    parser.add_argument("config", help="YAML config file containing results.")

    args = parser.parse_args()

    logging.basicConfig(
        format=" %(asctime)s: %(message)s",
        level=logging.DEBUG if args.debug else logging.INFO,
        stream=sys.stderr,
    )

    try:
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
