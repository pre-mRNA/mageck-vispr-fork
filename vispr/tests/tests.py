# coding: utf-8
from __future__ import absolute_import, division, print_function

import os, sys

from vispr.cli import init_server

def test_server():
    os.chdir("vispr/tests")
    print("Server starting. Please open http://localhost:5000 in your browser.", file=sys.stderr)
    with open("config.yaml") as f:
        config = yaml.load(f)
        init_server(config)
