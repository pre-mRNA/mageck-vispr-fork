# coding: utf-8
from __future__ import absolute_import, division, print_function

import re, json

import numpy as np
from flask import Flask, render_template, request, session

app = Flask(__name__)


@app.route("/")
def index():
    screen = request.args.get("screen")
    if screen:
        session["screen"] = screen
    return render_template("index.html",
                           screens=app.screens,
                           screen=get_screen())


@app.route("/targets/<selection>")
def targets(selection):
    return render_template("targets.html",
                           screens=app.screens,
                           selection=selection,
                           screen=get_screen(),
                           control_targets=get_screen().control_targets,
                           hide_control_targets=session.get("hide_control_targets", False))


@app.route("/qc")
def qc():
    return render_template("qc.html",
                           screens=app.screens,
                           screen=get_screen(),
                           fastqc=get_screen().fastqc is not None,
                           mapstats=get_screen().mapstats is not None)


@app.route("/compare")
def compare():
    overlap_items = ["{} {}".format(screen, sel)
                     for screen in app.screens for sel in "+-"]
    return render_template("compare.html",
                           screens=app.screens,
                           screen=get_screen(),
                           overlap_items=overlap_items)


@app.route("/plt/pvals/<selection>")
def plt_pvals(selection):
    plt = get_targets(selection).plot_pvals()
    return plt


@app.route("/plt/pvalhist/<selection>")
def plt_pval_hist(selection):
    plt = get_targets(selection).plot_pval_hist()
    return plt


@app.route("/tbl/targets/<selection>", methods=["GET", "POST"])
def tbl_targets(selection):
    offset = int(request.args.get("offset", 0))
    perpage = int(request.args.get("perPage", 20))

    # sort and slice records
    records = get_targets(selection)[:]
    total_count = records.shape[0]
    filter_count = total_count

    # restrict to overlap
    if "fdr" in request.form and "overlap-items" in request.form:
        pass

    search = get_search()
    if search or session.get("hide_control_targets", False):
        if search:
            filter = records["target"].str.contains(search)
        else:
            control_targets = get_screen().control_targets
            filter = records["target"].apply(lambda target: target not in control_targets)

        if np.any(filter):
            records = records[filter]
            filter_count = records.shape[0]
        else:
            return render_template("dyntable.json",
                                   records="[]",
                                   filter_count=0,
                                   total_count=total_count)

    columns, ascending = get_sorting()
    if columns:
        records = records.sort(columns, ascending=ascending)
    else:
        records = records.sort("p-value")
    records = records[offset:offset + perpage]

    def fmt_col(col):
        if col.dtype == np.float64:
            return col.apply("{:.2g}".format)
        return col

    records = records.apply(fmt_col)

    return render_template("dyntable.json",
                           records=records.to_json(orient="records",
                                                   double_precision=15),
                           filter_count=filter_count,
                           total_count=total_count)


@app.route("/tbl/pvals_highlight/<selection>/<targets>")
def tbl_pvals_highlight(selection, targets):
    targets = targets.split("|")
    records = get_targets(selection).get_pvals_highlight_targets(
        targets)
    return records.to_json(orient="records")


@app.route("/tbl/rnas/<target>")
def tbl_rnas(target):
    table = get_screen().rnas.by_target(target)
    return table.to_json(orient="records")
    return render_template("parcoords.json",
                           dimensions=json.dumps(list(table.columns)),
                           values=table.to_json(orient="values"))


@app.route("/plt/normalization")
def plt_normalization():
    plt = get_screen().rnas.plot_normalization()
    return plt


@app.route("/plt/pca/<int:x>/<int:y>/<int:legend>")
def plt_pca(x, y, legend):
    plt = get_screen().rnas.plot_pca(comp_x=x, comp_y=y, legend=legend == 1, )
    return plt


@app.route("/plt/correlation")
def plt_correlation():
    plt = get_screen().rnas.plot_correlation()
    return plt


@app.route("/plt/gc_content")
def plt_gc_content():
    plt = get_screen().fastqc.plot_gc_content()
    return plt


@app.route("/plt/base_quality")
def plt_base_quality():
    plt = get_screen().fastqc.plot_base_quality()
    return plt


@app.route("/plt/seq_quality")
def plt_seq_quality():
    plt = get_screen().fastqc.plot_seq_quality()
    return plt


@app.route("/plt/mapstats")
def plt_mapstats():
    plt = get_screen().mapstats.plot_mapstats()
    return plt


@app.route("/plt/zerocounts")
def plt_zerocounts():
    plt = get_screen().mapstats.plot_zerocounts()
    return plt


@app.route("/plt/overlap_chord", methods=["POST"])
def plt_overlap_chord():
    return app.screens.plot_overlap_chord(*get_overlap_args())


@app.route("/plt/overlap_venn", methods=["POST"])
def plt_overlap_venn():
    return app.screens.plot_overlap_venn(*get_overlap_args())


@app.route("/set/hide_control_targets/<int:value>")
def set_hide_control_targets(value):
    session["hide_control_targets"] = value == 1
    return ""


def get_overlap_args():
    def parse_item(item):
        screen, sel = item.split()
        return screen, sel == "+"

    if "fdr" not in request.form and "overlap-items" not in request.form:
        return None

    fdr = float(request.form.get("fdr", 0.25))
    items = list(map(parse_item, request.form.getlist("overlap-items")))
    return fdr, items


def get_sorting(pattern=re.compile("sorts\[(?P<col>.+)\]")):
    cols, ascending = [], []
    for arg, val in request.args.items():
        m = pattern.match(arg)
        if m:
            cols.append(m.group("col"))
            ascending.append(int(val) == 1)
    return cols, ascending


def get_search(pattern=re.compile("search\[(?P<target>.+)\]")):
    return request.args.get("queries[search]", None)


def get_screen():
    return app.screens[session.get("screen", next(iter(app.screens)))]


def get_targets(selection):
    return get_screen().targets(selection == "positive")
