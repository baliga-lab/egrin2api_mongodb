#!/usr/bin/env python3
import logging
import json
import os

from collections import defaultdict

from flask import Flask, Response, url_for, redirect, render_template, request, session, flash, jsonify
import flask

# because we have an API, we need to allow cross-origin here
from flask_cors import CORS, cross_origin

from werkzeug import secure_filename
import requests
import pymongo

app = Flask(__name__)
CORS(app)

app.config.from_envvar('EGRIN2API_SETTINGS')
client = pymongo.MongoClient(host=app.config['MONGODB_HOST'], port=app.config['MONGODB_PORT'])
db = client[app.config['MONGODB_DB']]

def make_sites(cluster_id, motif_num, start, stop):
    sites = db.fimo.find({"cluster_id": cluster_id, "motif_num": motif_num,
                              "start": {"$gte": start}, "stop": {"$lte": stop}},
                             {"_id": 0, "start": 1, "stop": 1, "strand": 1})
    result = []
    for s in sites:
        result.append({"start": s["start"], "stop": s["stop"], "strand": s["strand"]})
    return result

def unique(l):
    added = set()
    result = []
    for e in l:
        key = "%s:%s-%s" % (e["strand"], e["start"], e["stop"])
        if key not in added:
            added.add(key)
            result.append(e)
    return result

def make_counts(docs):
    counts = defaultdict(int)
    for doc in docs:
        for i in range(doc['start'] + 1, doc['stop'] + 1):
            counts[i] += 1
    result = sorted([{'pos': x, 'count': count } for x, count in counts.items()], key=lambda d: d['pos'])
    pos0 = result[0]['pos']
    posn = result[-1]['pos']
    # add artificial 0 GRE counts for a prettier presentation
    result.insert(0, { 'pos': pos0 - 1, 'count': 0})
    result.append({ 'pos': posn + 1, 'count': 0})
    return result


@app.route('/api/v1.0.0/gene_gres/<gene>')
def gene_gre_counts(gene):
    row_ids = [(r["row_id"], r["start"], r["stop"], r["strand"])
                   for r in db.row_info.find({"sysName": gene},
                                                 {"_id": 0, "row_id": 1, "start": 1, "stop": 1, "strand": 1})]
    gene_num, gene_start, gene_stop, gene_strand = row_ids[0]
    biclusters = db.bicluster_info.find({"rows": gene_num}, {"_id": 1})
    cluster_ids = [c['_id'] for c in biclusters]
    motif_infos = db.motif_info.find({"$and": [{"cluster_id": {"$in": cluster_ids}},
                                                   {"gre_id": {"$ne": "NaN"}},
                                                   {"gre_id": {"$ne": None}}]},
                                         {"_id": 0, "gre_id": 1, "motif_num": 1, "cluster_id": 1})
    gres = defaultdict(list)

    window_start = gene_start - 1000
    window_stop = gene_stop + 2000
    for m in motif_infos:
        gres[m["gre_id"]].extend(make_sites(m["cluster_id"], m["motif_num"], window_start, window_stop))
    gres = {'GRE_%d' % gre_id: make_counts(sorted(unique(sites), key=lambda e: e['start']))
            for gre_id, sites in gres.items() if len(sites) > 0}
    return jsonify(gene=gene, gres=gres)


@app.route('/api/v1.0.0/corem_genes/<corem_num>')
def cluster_genes(corem_num):
    corem_num = int(corem_num)
    corems = [corem for corem in db.corem.find({'corem_id': corem_num},{"_id": 0, "rows": 1})]
    rows = []
    if len(corems) > 0:
        rows = corems[0]["rows"]
    genes = db.row_info.find({"row_id": {"$in": rows}},
                             {"_id": 0, "row_id": 1, "sysName": 1, "start": 1, "stop": 1, "strand": 1, "accession": 1})
    gene_names = [gene["sysName"] for gene in genes]
    return jsonify(corem=corem_num, genes=gene_names)

@app.route('/')
def index():
    return render_template('index.html')


if __name__ == '__main__':
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)

    app.debug = True
    app.secret_key = 'supercalifragilistic'
    app.logger.addHandler(handler)
    app.run(host='0.0.0.0', debug=True)
