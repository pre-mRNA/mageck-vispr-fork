# coding: utf-8
from __future__ import absolute_import, division, print_function

import argparse
import logging
import sys
import string
import random
import os

import yaml

from vispr import TargetResults, Results, VisprError
from vispr.server import app


def init_server(config):
    app.results = Results(config)
    app.secret_key = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(30))
    app.run()


def main():
    # create arg parser
    parser = argparse.ArgumentParser("An HTML5-based interactive visualization of CRISPR/Cas9 screen data.")
    parser.add_argument("--debug", action="store_true", help="Print debug info.")
    subparsers = parser.add_subparsers(dest="subcommand")

    server = subparsers.add_parser("server")
    server.add_argument("config", help="YAML config file containing results.")

    test = subparsers.add_parser("test")

    args = parser.parse_args()

    logging.basicConfig(
        format="%(message)s",
        level=logging.DEBUG if args.debug else logging.INFO,
        stream=sys.stderr,
    )

    try:
        if args.subcommand == "server":
            with open(args.config) as f:
                config = yaml.load(f)
                init_server(config)
        elif args.subcommand == "test":
            os.chdir(os.path.join(os.path.dirname(__file__), "tests"))
            with open("config.yaml") as f:
                config = yaml.load(f)
                init_server(config)
        else:
            parser.print_help()
            exit(1)
    except VisprError as e:
        logging.error(e)
        exit(1)
    except ImportError as e:
        print(
            "{}. Please ensure that all depencies from "
            "requirements.txt are installed.".format(e),
            file=sys.stderr
        )
        exit(1)
    exit(0)
