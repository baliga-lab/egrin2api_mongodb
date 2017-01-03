#!/usr/bin/env python
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
from bson.objectid import ObjectId
import pandas as pd
import numpy as np

import egrin2_query as e2q

BATCH_SIZE = 100

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


@app.route('/api/v1.0.0/corems_with_gene/<gene>')
def corems_with_gene(gene):
    gene_ids = [r["row_id"] for r in db.row_info.find({"sysName": gene}, {"_id": 0, "row_id": 1})]
    gene_id = gene_ids[0]
    corem_infos = [{"corem_id": r["corem_id"], "genes": r["rows"]}
                       for r in db.corem.find({"rows": gene_id}, {"_id": 0, "corem_id": 1, "rows": 1})]
    return jsonify(corem_infos=corem_infos)


@app.route('/api/v1.0.0/gene_gres/<gene>')
def gene_gre_counts(gene):
    """Get the GRE counts for a specific gene
    TODO: make corem contexts available
    """
    gene_infos = [(r["row_id"], r["start"], r["stop"], r["strand"])
                    for r in db.row_info.find({"sysName": gene},
                                                  {"_id": 0, "row_id": 1, "start": 1, "stop": 1, "strand": 1})]
    gene_num, gene_start, gene_stop, gene_strand = gene_infos[0]
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


@app.route('/api/v1.0.0/corem_info/<corem_num>')
def corem_info(corem_num):
    """
    TODO: add the counts for the GREs in the corem associated with the corem's genes
    """
    genes = e2q.corem_genes(db, corem_num)
    corem_gres = e2q.agglom(db, genes, x_type='genes', y_type='gres', logic='or')
    gre_ids = sorted([int(i) for i in corem_gres.index])
    return jsonify(corem=corem_num, genes=genes, gres=gre_ids)


@app.route('/api/v1.0.0/condition_info/<condition_id>')
def condition_info(condition_id):
    cond_docs = [{ "id": c['col_id'], "name": c['egrin2_col_name']}
                  for c in db.col_info.find({'col_id': int(condition_id)},
                                            {'_id': 0, 'col_id': 1, 'egrin2_col_name': 1})]
    cond_blocks_path = app.config["COND_BLOCKS_FILE"]
    df = pd.read_csv(cond_blocks_path)
    cond2block = defaultdict(set)
    block2cond = defaultdict(set)

    for i in range(df.shape[0]):
        egrin2_block = df['EGRIN2.block'][i]
        conds = df["conds"][i].split()
        for cond in conds:
            cond2block[cond].add(egrin2_block)
            block2cond[egrin2_block].add(cond)
    cond_blocks = {block: list(block2cond[block]) for block in cond2block[conds[0]]}
    return jsonify(condition=cond_docs[0], blocks=cond_blocks)


@app.route('/api/v1.0.0/gene_info/<gene>')
def gene_info(gene):
    chroms = db.genome.find({}, {"_id": 0, "scaffoldId": 1, "NCBI_RefSeq": 1 })
    chrom_map = { int(c['scaffoldId']): c['NCBI_RefSeq'] for c in chroms }
    genes = [{ "id": r['row_id'],
                   "gene_name": r['sysName'],
                   "common_name": r['name'],
                   "accession":  str(r['accession']),
                   "description": r['desc'],
                   "start": r['start'], "stop": r['stop'], "strand": r['strand'],
                   "chromosome": chrom_map[r['scaffoldId']]}
                  for r in db.row_info.find({"sysName": gene},
                                            {'_id': 0, 'row_id': 1, 'sysName': 1, 'name': 1,
                                                 'accession': 1, 'desc': 1, 'start': 1, 'stop': 1,
                                            'strand': 1, 'scaffoldId': 1})]
    return jsonify(gene=genes[0])


@app.route('/api/v1.0.0/gene_biclusters/<gene>')
def gene_biclusters(gene):
    row_id = db.row_info.find({'sysName': 'Rv0001'}, {'_id': 0, 'row_id': 1})[0]["row_id"]
    clusters = [{ "id": str(c['_id']), "num_genes": len(c['rows']), "num_conditions": len(c["columns"]), "residual": c["residual"]}
                    for c in db.bicluster_info.find({"rows": row_id}, {'_id': 1, 'rows': 1, 'columns': 1, 'residual': 1})]
    return jsonify(biclusters=clusters)


@app.route('/api/v1.0.0/condition_biclusters/<condition_id>')
def condition_biclusters(condition_id):
    clusters = [{ "id": str(c['_id']), "num_genes": len(c['rows']), "num_conditions": len(c["columns"]), "residual": c["residual"]}
                    for c in db.bicluster_info.find({"columns": int(condition_id)}, {'_id': 1, 'rows': 1, 'columns': 1, 'residual': 1})]
    return jsonify(biclusters=clusters)


@app.route('/api/v1.0.0/bicluster_info/<cluster_id>')
def bicluster_info(cluster_id):
    cluster_id = ObjectId(cluster_id)
    cluster = [{ "id": str(c['_id']), "num_genes": len(c['rows']), "num_conditions": len(c["columns"]), "residual": c["residual"]}
                   for c in db.bicluster_info.find({"_id": cluster_id}, {'_id': 1, 'rows': 1, 'columns': 1, 'residual': 1})][0]
    return jsonify(bicluster=cluster)


@app.route('/api/v1.0.0/bicluster_genes/<cluster_id>')
def bicluster_genes(cluster_id):
    cluster_id = ObjectId(cluster_id)
    gene_ids = db.bicluster_info.find({"_id": cluster_id}, {'_id': 0, 'rows': 1})[0]["rows"]
    chroms = db.genome.find({}, {"_id": 0, "scaffoldId": 1, "NCBI_RefSeq": 1 })
    chrom_map = { int(c['scaffoldId']): c['NCBI_RefSeq'] for c in chroms }
    genes = [{ "id": r['row_id'],
                   "gene_name": r['sysName'],
                   "common_name": r['name'],
                   "accession":  str(r['accession']),
                   "description": r['desc'],
                   "start": r['start'], "stop": r['stop'], "strand": r['strand'],
                   "chromosome": chrom_map[r['scaffoldId']]}
                  for r in db.row_info.find({"row_id": { "$in": gene_ids }},
                                            {'_id': 0, 'row_id': 1, 'sysName': 1, 'name': 1,
                                                 'accession': 1, 'desc': 1, 'start': 1, 'stop': 1,
                                            'strand': 1, 'scaffoldId': 1})]
    return jsonify(genes=genes)


@app.route('/api/v1.0.0/bicluster_conditions/<cluster_id>')
def bicluster_conditions(cluster_id):
    cluster_id = ObjectId(cluster_id)
    condition_ids = db.bicluster_info.find({"_id": cluster_id}, {'_id': 0, 'columns': 1})[0]["columns"]
    cond_docs = [{ "id": c['col_id'], "name": c['egrin2_col_name']}
                    for c in db.col_info.find({'col_id': { "$in": condition_ids }},
                                                  {'_id': 0, 'col_id': 1, 'egrin2_col_name': 1})]
    return jsonify(conditions=cond_docs)


@app.route('/api/v1.0.0/bicluster_pssms/<cluster_id>')
def bicluster_pssms(cluster_id):
    cluster_id = ObjectId(cluster_id)
    motifs = []
    for mi in db.motif_info.find({"cluster_id": cluster_id}, {"_id": 0, "pwm": 1, "motif_num": 1}):
        motifs.append({"motif_num": mi["motif_num"], "alphabet": ["A", "C", "G", "T"],
                           "values": [[r["a"], r["c"], r["g"], r["t"]] for r in mi["pwm"]]})
    return jsonify(motifs=motifs)


@app.route('/api/v1.0.0/corem_genes/<corem_id>')
def corem_genes(corem_id):
    gene_ids = db.corem.find({"corem_id": int(corem_id)}, {'_id': 0, 'rows': 1})[0]["rows"]
    chroms = db.genome.find({}, {"_id": 0, "scaffoldId": 1, "NCBI_RefSeq": 1 })
    chrom_map = { int(c['scaffoldId']): c['NCBI_RefSeq'] for c in chroms }
    genes = [{ "id": r['row_id'],
                   "gene_name": r['sysName'],
                   "common_name": r['name'],
                   "accession":  str(r['accession']),
                   "description": r['desc'],
                   "start": r['start'], "stop": r['stop'], "strand": r['strand'],
                   "chromosome": chrom_map[r['scaffoldId']]}
                  for r in db.row_info.find({"row_id": { "$in": gene_ids }},
                                            {'_id': 0, 'row_id': 1, 'sysName': 1, 'name': 1,
                                                 'accession': 1, 'desc': 1, 'start': 1, 'stop': 1,
                                            'strand': 1, 'scaffoldId': 1})]
    return jsonify(genes=genes)


@app.route('/api/v1.0.0/corem_conditions/<corem_id>')
def corem_conditions(corem_id):
    conds = db.corem.find({"corem_id": int(corem_id)}, {'_id': 0, 'cols': 1})[0]["cols"]
    cond_ids = [int(c["col_id"]) for c in conds]
    conds = [{ "id": c['col_id'], "name": c['egrin2_col_name']}
                 for c in db.col_info.find({"col_id": { "$in": cond_ids} }, {'_id': 0, 'col_id': 1, 'egrin2_col_name': 1})]
    return jsonify(conditions=conds)


@app.route('/api/v1.0.0/corem_expressions/<corem_id>')
def corem_expressions(corem_id):
    """retrieve the co-regulation expressions for the specified corem"""
    expressions = []
    corem = db.corem.find({"corem_id": int(corem_id)}, {'_id': 0, 'cols': 1, 'rows': 1})[0]
    cols = [int(c['col_id']) for c in corem['cols']]
    gene_ids = corem["rows"]
    expressions.append({"cols": cols, "rows": gene_ids})
    exps = db.gene_expression.find({'col_id': {"$in": cols}, 'row_id': {"$in": gene_ids}},
                                    {'_id': 0, 'row_id': 1, 'col_id': 1, 'standardized_expression': 1, 'raw_expression': 1})
    exps.sort('col_id', 1)
    sys_names = db.row_info.find({"row_id": {"$in": gene_ids}}, {"_id": 0, "row_id": 1, "sysName": 1})
    sys_name_map = {sn['row_id']: sn['sysName'] for sn in sys_names}
    col_names = db.col_info.find({"col_id": {"$in": cols}}, {"_id": 0, "col_id": 1, "egrin2_col_name": 1})
    col_name_map = {cn['col_id']: cn['egrin2_col_name'] for cn in col_names}
    series = defaultdict(list)
    for exp in exps:
        gene = sys_name_map[exp['row_id']]
        value = exp['standardized_expression']
        if np.isnan(value):
            value = 0.0
        series[gene].append(value)

    return jsonify(expressions=series, conditions=[col_name_map[cid] for cid in cols])


######################################################################
### API functions global to the model
######################################################################

@app.route('/api/v1.0.0/corems')
def corems():
    """return a list of entries (corem_id, num_genes, num_conditions, num_gres)"""
    corems = [{ "id": c['corem_id'], "num_conds": len(c['cols']), "num_genes": len(c['rows'])}
                  for c in db.corem.find({}, {'_id': 0, 'corem_id': 1, 'cols': 1, 'rows': 1})]
    return jsonify(corems=corems)


@app.route('/api/v1.0.0/conditions')
def conditions():
    conds = [{ "id": c['col_id'], "name": c['egrin2_col_name']}
                  for c in db.col_info.find({}, {'_id': 0, 'col_id': 1, 'egrin2_col_name': 1})]
    return jsonify(conditions=conds)


def _query_genes(start, num_entries):
    """reusable gene query function. Since this is a database access
    function, we could possibly move that into a database access module
    """
    end = start + num_entries
    chroms = db.genome.find({}, {"_id": 0, "scaffoldId": 1, "NCBI_RefSeq": 1 })
    chrom_map = { int(c['scaffoldId']): c['NCBI_RefSeq'] for c in chroms }
    cursor = db.row_info.find({}, {'_id': 0, 'row_id': 1, 'sysName': 1, 'name': 1,
                                       'accession': 1, 'desc': 1, 'start': 1, 'stop': 1,
                                       'strand': 1, 'scaffoldId': 1})[start:end]
    return [{ "id": r['row_id'],
                  "gene_name": r['sysName'],
                  "common_name": r['name'],
                  "accession":  str(r['accession']),
                  "description": r['desc'],
                  "start": r['start'], "stop": r['stop'], "strand": r['strand'],
                  "chromosome": chrom_map[r['scaffoldId']]} for r in cursor]


@app.route('/api/v1.0.0/genes')
def genes():
    """Returns a list of genes. The maximum number of entries returned is BATCH_SIZE.
    Optional request parameters:
    - start: start position
    - length: number of elements to return
    """
    try:
        start = int(request.args.get('start'))
    except Exception as e:
        start = 0
    try:
        num_entries = int(request.args.get('length'))
    except Exception as e:
        num_entries = BATCH_SIZE
    genes = _query_genes(start, num_entries)
    return jsonify(genes=genes)


@app.route('/api/v1.0.0/biclusters')
def biclusters():
    """Returns a list of biclusters. The maximum number of entries returned is BATCH_SIZE.
    Optional request parameters:
    - start: start position
    - length: number of elements to return
    """
    try:
        start = int(request.args.get('start'))
    except Exception as e:
        start = 0
    try:
        num_entries = int(request.args.get('length'))
    except Exception as e:
        num_entries = BATCH_SIZE
    end = start + num_entries

    cursor = db.bicluster_info.find({}, {'_id': 1, 'rows': 1, 'columns': 1, 'residual': 1})[start:end]
    clusters = [{ "id": str(c['_id']), "num_genes": len(c['rows']), "num_conditions": len(c["columns"]), "residual": c["residual"]}
                    for c in cursor]
    return jsonify(biclusters=clusters)


@app.route('/api/v1.0.0/api_info')
def api_info():
    db_info = list(db.ensemble_info.find({}))
    return jsonify(model_type='egrin2', api_version='1.0.0', num_runs=len(db_info),
                   organism=db_info[0]['organism'], species=db_info[0]['species'])


@app.route('/api/v1.0.0/summary')
def summary():
    num_genes = db.row_info.count()
    num_conditions = db.col_info.count()
    num_corems = db.corem.count()
    num_biclusters = db.bicluster_info.count()

    #num_gres = len(db.motif_info.distinct('gre_id', { '$and': [ {'gre_id': {'$ne': "NaN"} }, {'gre_id': {'$ne': None}}]} ));
    num_gres = 163  # TODO we need a link to show the selection of the relevant GREs

    return jsonify(num_genes=num_genes, num_conditions=num_conditions, num_corems=num_corems,
                   num_biclusters=num_biclusters, num_gres=num_gres)


@app.route('/')
def index():
    return render_template('index.html')


"""
TODO:

  - GREs are only found in db.motif_info as gre_id attribute and it can be null
  - cluster-motif assosciations are stored in motif_info and fimo
  - corem - cluster relationships ?
  - corem - GRE relationships ?
"""

if __name__ == '__main__':
    handler = logging.StreamHandler()
    handler.setLevel(logging.INFO)

    app.debug = True
    app.secret_key = 'supercalifragilistic'
    app.logger.addHandler(handler)
    app.run(host='0.0.0.0', debug=True)
