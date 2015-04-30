# coding: utf-8
from __future__ import absolute_import, division, print_function

import re, json

from flask import Flask, render_template, request, session

app = Flask(__name__)


@app.route("/")
def index():
    screen = request.args.get("screen")
    if screen:
        session["screen"] = screen
    return render_template("index.html", screens=app.results.screens, screen=get_screen())


@app.route("/targets/<selection>")
def targets(selection):
    return render_template(
        "targets.html",
        screens=app.results.screens,
        screen=get_screen(),
        selection=selection,
    )


@app.route("/plt/pvals/<selection>")
def plt_pvals(selection):
    return app.results.targets[get_screen()].plot_pvals(positive=selection == "positive").to_json()


@app.route("/idx/pvals/<selection>/<target>")
def idx_pvals(selection, target):
    return str(app.results.targets[get_screen()].get_pvals_idx(target, positive=selection == "positive"))


@app.route("/plt/pvalhist/<selection>")
def plt_pval_hist(selection):
    plt = app.results.targets[get_screen()].plot_pval_hist(positive=selection == "positive").to_json()
    return plt


@app.route("/tbl/targets/<selection>", methods=["GET"])
def tbl_targets(selection):
    offset = int(request.args.get("offset", 0))
    perpage = int(request.args.get("perPage", 20))

    # sort and slice records
    records = app.results.targets[get_screen()][:]
    total_count = records.shape[0]
    filter_count = total_count  # TODO add filtering

    columns, ascending = get_sorting()
    if columns:
        records = records.sort(columns, ascending=ascending)
    else:
        records = records.sort("p.pos" if selection == "positive" else "p.neg")
    records = records[offset:offset + perpage]

    return render_template(
        "dyntable.json",
        records=records.to_json(orient="records", double_precision=15),
        filter_count=filter_count,
        total_count=total_count,
    )


@app.route("/tbl/rnas/<target>")
def tbl_rnas(target):
    table = app.results.rnas[get_screen()].by_target(target)
    return table.to_json(orient="records")
    print(    table.columns)
    return render_template(
        "parcoords.json",
        dimensions=json.dumps(list(table.columns)),
        values=table.to_json(orient="values")
    )


def get_sorting(pattern=re.compile("sorts\[(?P<col>.+)\]")):
    cols, ascending = [], []
    for arg, val in request.args.items():
        m = pattern.match(arg)
        if m:
            cols.append(m.group("col"))
            ascending.append(int(val) == 1)
    return cols, ascending


def get_screen():
    return session.get("screen", next(iter(app.results.screens)))
