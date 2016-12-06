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
from bson.objectid import ObjectId
import pandas as pd

import egrin2_query as e2q

app = Flask(__name__)
CORS(app)

app.config.from_envvar('EGRIN2API_SETTINGS')
client = pymongo.MongoClient(host=app.config['MONGODB_HOST'], port=app.config['MONGODB_PORT'])
db = client[app.config['MONGODB_DB']]


"""
for i in range(len(corem_ids)):
    try:
        genes = find_corem_info(db, [corem_ids.values[i,0]], x_type="corem_id", y_type="genes")
        genes_list = list(genes['genes'])
        corem_GREs = agglom(db, genes_list, x_type="genes", y_type="gres", logic="or")
        corem_GREs_sig = corem_GREs[corem_GREs.loc[:,'qval_BH']<0.01]
        corem_GREs_sig['COREM'] = str(corem_ids.values[i,0])
        corem_GREs_sig['genes'] = str(genes_list)
        #corem_GREs_sig.columns = ['GRE_id', 'counts', 'all_counts', 'pval', 'qval_BH', 'qval_Bonf', 'COREM', 'genes']
        corem_GREs_sig_163 = corem_GREs_sig[corem_GREs_sig.index<=163]
        corem_GREs_sig_163.to_csv('gres/corem_gres2/'+str(corem_ids.values[i,0])+'.csv')
    except:
        # note that errors are silently ignored
        pass
"""


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
    print("gene id: ", gene_id)
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
    print("cond: %s, res: %s" % (condition_id, str(cond_docs)))

    df = pd.read_csv("cond_blocks.csv")
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


@app.route('/api/v1.0.0/genes')
def genes():
    chroms = db.genome.find({}, {"_id": 0, "scaffoldId": 1, "NCBI_RefSeq": 1 })
    chrom_map = { int(c['scaffoldId']): c['NCBI_RefSeq'] for c in chroms }
    genes = [{ "id": r['row_id'],
                   "gene_name": r['sysName'],
                   "common_name": r['name'],
                   "accession":  str(r['accession']),
                   "description": r['desc'],
                   "start": r['start'], "stop": r['stop'], "strand": r['strand'],
                   "chromosome": chrom_map[r['scaffoldId']]}
                  for r in db.row_info.find({},
                                            {'_id': 0, 'row_id': 1, 'sysName': 1, 'name': 1,
                                                 'accession': 1, 'desc': 1, 'start': 1, 'stop': 1,
                                            'strand': 1, 'scaffoldId': 1})]
    return jsonify(genes=genes)


@app.route('/api/v1.0.0/biclusters')
def biclusters():
    """TODO: we currently limit the output to 100 entries, because rendering everything exceeds
    the memory limit, this lsit should most likely be loaded in batches
    """
    clusters = [{ "id": str(c['_id']), "num_genes": len(c['rows']), "num_conditions": len(c["columns"]), "residual": c["residual"]}
                    for c in db.bicluster_info.find({}, {'_id': 1, 'rows': 1, 'columns': 1, 'residual': 1})]
    return jsonify(biclusters=clusters[:100])


@app.route('/api/v1.0.0/api_info')
def api_info():
    db_info = list(db.ensemble_info.find({}))
    return jsonify(model_type='egrin2', api_version='1.0.0', num_runs=len(db_info),
                   organism=db_info[0]['organism'], species=db_info[0]['species'])


@app.route('/api/v1.0.0/summary')
def summary():
    #db_info = list(db.ensemble_info.find({}))
    num_genes = db.row_info.count()
    num_conditions = db.col_info.count()
    num_corems = db.corem.count()
    num_biclusters = db.bicluster_info.count()

    num_gres = len(db.motif_info.distinct('gre_id', { '$and': [ {'gre_id': {'$ne': "NaN"} }, {'gre_id': {'$ne': None}}]} ));

    return jsonify(num_genes=num_genes, num_conditions=num_conditions, num_corems=num_corems,
                   num_biclusters=num_biclusters, num_gres=num_gres)


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
