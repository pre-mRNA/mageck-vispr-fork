from flask import Flask, render_template, request

app = Flask(__name__)

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/genes")
def genes():
    return render_template("genes.html")


@app.route("/plt/rra/<type>")
def plt_rra(type):
    json = app.gene_results.plot_rra(positive=type == "pos").to_json()
    return json

@app.route("/tbl/genes", methods=["GET"])
def tbl_genes():
    offset = int(request.args.get("offset", 0))
    perpage = int(request.args.get("perPage", 20))
    return render_template("table.json", records=app.gene_results[offset:offset + perpage].to_json(orient="records"))
