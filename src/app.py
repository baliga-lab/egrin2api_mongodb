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
import mysql.connector

import egrin2_query as e2q

BATCH_SIZE = 100

app = Flask(__name__)
CORS(app)

app.config.from_envvar('EGRIN2API_SETTINGS')
client = pymongo.MongoClient(host=app.config['MONGODB_HOST'], port=app.config['MONGODB_PORT'])
db = client[app.config['MONGODB_DB']]

def __request_batch_params():
    try:
        start = int(request.args.get('start'))
    except Exception as e:
        start = 0
    try:
        num_entries = int(request.args.get('length'))
    except Exception as e:
        num_entries = BATCH_SIZE
    try:
        search = request.args.get('search')
    except Exception as e:
        search = None
    return start, num_entries, search


def make_sites(cluster_id, motif_num, start, stop):
    """Generate the site objects for the specified motif and number. Note the explicit
    int() cast to make sure our CentOS server pymongo can understand it
    """
    sites = db.fimo.find({"cluster_id": cluster_id, "motif_num": motif_num,
                              "start": {"$gte": int(start)}, "stop": {"$lte": int(stop)}},
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
    corem_infos = [{"id": r["corem_id"], "genes": r["rows"]}
                       for r in db.corem.find({"rows": gene_id}, {"_id": 0, "corem_id": 1, "rows": 1})]
    return jsonify(corem_infos=corem_infos)


def cachedb_conn():
    return mysql.connector.connect(host=app.config['MYSQL_HOST'],
                                       user=app.config['MYSQL_USER'],
                                       password=app.config['MYSQL_PASSWORD'],
                                       database=app.config['MYSQL_CACHEDB'])


def __gene_gre_counts_mysql(gene, start, stop):
    conn = cachedb_conn()
    cursor = conn.cursor()
    result = defaultdict(list)
    try:
        cursor.execute("""select gre_id,position,count from gene_gre_counts ggc
join genes g on ggc.gene_id=g.id where g.name=%s and position >= %s and position <= %s
order by gre_id,position""",
                           [gene, int(start), int(stop)])
        for row in cursor.fetchall():
            gre_id, pos, count = row
            result['GRE_' + str(gre_id)].append({"pos": pos, "count": count})
        return result
    finally:
        cursor.close()
        conn.close()

@app.route('/api/v1.0.0/gene_gres/<gene>')
def gene_gre_counts(gene):
    """Get the GRE counts (we only use the first 163) for a specific gene
    """

    # Read the output window from the TSS start stop file
    gene_start_stop_path = app.config["GENE_TSS_START_STOP_FILE"]
    df = pd.read_csv(gene_start_stop_path)
    df = df[df['gene'] == gene].reset_index()
    tss_start = df['start'][0]
    tss_stop = df['end'][0]
    tss_strand = df['strand'][0]

    # Read the chipseq peaks for the gene
    chipseq_path = app.config["CHIPSEQ_FILE"]
    chipseq_df = pd.read_csv(chipseq_path)

    gene_infos = [(r["row_id"], r["start"], r["stop"], r["strand"])
                    for r in db.row_info.find({"sysName": gene},
                                                  {"_id": 0, "row_id": 1, "start": 1, "stop": 1, "strand": 1})]
    gene_num, gene_start, gene_stop, gene_strand = gene_infos[0]
    biclusters = db.bicluster_info.find({"rows": gene_num}, {"_id": 1})
    cluster_ids = [c['_id'] for c in biclusters]
    motif_infos = db.motif_info.find({"$and": [{"cluster_id": {"$in": cluster_ids}},
                                                   {"gre_id": {"$ne": "NaN"}},
                                                   {"gre_id": {"$ne": None}},
                                                {"gre_id": {"$lte": 163 }}]},
                                         {"_id": 0, "gre_id": 1, "motif_num": 1, "cluster_id": 1})
    # set the window to the TSS, and if the TSS does not overlap
    # with the coding region, include 50 additional bases
    # TODO: because of the circular nature of the chromosome, we need to handle Rv0001
    # separately
    window_start = tss_start
    window_stop = tss_stop
    if gene_start > window_stop:
        window_stop = gene_start + 50
    app.logger.debug("window %d-%d", window_start, window_stop)
    chipseq_df = chipseq_df[(chipseq_df['position'] >= window_start) & (chipseq_df['position'] <= window_stop)].reset_index()
    chipseq_peaks = {e['tf']: int(e['position']) for index, e in chipseq_df.iterrows()}

    gres = __gene_gre_counts_mysql(gene, window_start, window_stop)
    return jsonify(gene={'name': gene, 'start': gene_start, 'stop': gene_stop, 'strand': gene_strand},
                   tss={'start': int(tss_start), 'stop': int(tss_stop)},
                       gres=gres, chipseq_peaks=chipseq_peaks)


@app.route('/api/v1.0.0/corem_info/<corem_num>')
def corem_info(corem_num):
    genes = e2q.corem_genes(db, corem_num)
    corem_gres = e2q.gene_gres(db, genes)
    gre_ids = sorted([int(i) for i in corem_gres.index if i <= 163])
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

def __condition_block_map():
    cond_blocks_path = app.config["COND_BLOCKS_FILE"]
    df = pd.read_csv(cond_blocks_path)
    block2cond = defaultdict(set)

    for i in range(df.shape[0]):
        egrin2_block = df['EGRIN2.block'][i]
        conds = df["conds"][i].split()
        for cond in conds:
            block2cond[egrin2_block].add(cond)
    result = {}
    for block, conds in block2cond.items():
        result[block] = list(conds)
    return result


@app.route('/api/v1.0.0/condition_blocks')
def condition_block_info():
    return jsonify(__condition_block_map())


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
    start, num_entries, search = __request_batch_params()
    end = start + num_entries

    gene_ids = db.corem.find({"corem_id": int(corem_id)}, {'_id': 0, 'rows': 1})[0]["rows"]
    chroms = db.genome.find({}, {"_id": 0, "scaffoldId": 1, "NCBI_RefSeq": 1 })
    chrom_map = { int(c['scaffoldId']): c['NCBI_RefSeq'] for c in chroms }

    search_params = {"row_id": { "$in": gene_ids }}
    if search is not None and len(search) > 0:
        match_name = {'$or': [{'name': {'$regex': search}},  {'sysName': {'$regex': search}}] }
        search_params = {"$and": [search_params, match_name]}

    cursor = db.row_info.find(search_params,
                                  {'_id': 0, 'row_id': 1, 'sysName': 1, 'name': 1,
                                       'accession': 1, 'desc': 1, 'start': 1, 'stop': 1,
                                  'strand': 1, 'scaffoldId': 1})
    total = cursor.count()
    cursor = cursor[start:end]
    genes = [{ "id": r['row_id'],
                   "gene_name": r['sysName'],
                   "common_name": r['name'],
                   "accession":  str(r['accession']),
                   "description": r['desc'],
                   "start": r['start'], "stop": r['stop'], "strand": r['strand'],
                   "chromosome": chrom_map[r['scaffoldId']]}
                  for r in cursor]
    return jsonify(genes=genes, total=total)


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
    blocks = request.args.get('blocks')
    blocks = map(int, blocks.split(','))

    if (len(blocks) > 1 or (len(blocks) == 1 and blocks[0] != 0)):
        # find the condition ids that are needed
        corem_cond_blocks = {b['id']: b['name'] for b in __corem_condition_blocks(corem_id)}
        block_names = [corem_cond_blocks[block_id] for block_id in blocks]
        condblock_conds = __condition_block_map()
        used_conds = set()
        for name in block_names:
            used_conds.update(condblock_conds[name])
        all_cols = db.col_info.find({}, {'_id': 0, 'col_id': 1, 'egrin2_col_name': 1})
        all_col_names = {cn['egrin2_col_name']: cn['col_id'] for cn in all_cols}
        used_cond_ids = sorted([all_col_names[c] for c in used_conds])
    else:
        used_cond_ids = None


    corem = db.corem.find({"corem_id": int(corem_id)}, {'_id': 0, 'cols': 1, 'rows': 1})[0]
    if used_cond_ids is None:
        cols = [int(c['col_id']) for c in corem['cols']]
    else:
        cols = used_cond_ids
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

def __corem_condition_blocks(corem_id):
    """reusable blocks for the specified corem"""
    cond_blocks_path = app.config["COREM_COND_BLOCKS_FILE"]
    df = pd.read_csv(cond_blocks_path)
    df = df[df['COREM'] == int(corem_id)].reset_index()
    # we give the blocks a fake id, so we can easily reload
    # expressions
    return [{'id': i + 1,
             'name': df['EGRIN2.block'][i], 'q_value': df['BH.adjusted.p.value'][i]}
            for i in range(df.shape[0])]


@app.route('/api/v1.0.0/corem_condition_enrichment/<corem_id>')
def corem_condition_enrichment(corem_id):
    return jsonify(condition_blocks=__corem_condition_blocks(corem_id))


@app.route('/api/v1.0.0/corem_categories/<corem_id>')
def corem_category_enrichment(corem_id):
    corem_category_path = app.config["COREM_CATEGORIES_FILE"]
    df = pd.read_csv(corem_category_path)
    df = df[df['corem_ID'] == int(corem_id)].reset_index()
    categories = [{'category': df['tuberculist.category'][i], 'p_adj': df['p.adj'][i]}
                   for i in range(df.shape[0])]
    return jsonify(categories=categories)


def __pwm2pssm_rows(pwm):
    """reduce the information to a 2D-array of numbers"""
    return [[r['a'], r['c'], r['g'], r['t']] for r in pwm]

def __gre_motifs(gre_id):
    """extract the motifs from the GRE"""
    motifs = db.motif_info.find({'gre_id': gre_id})
    return [{'pssm': __pwm2pssm_rows(m['pwm'])} for m in motifs]


def __gre_pssm(gre_id):
    gre_pssms_path = app.config["GRE_PSSMS_DIR"]
    path = os.path.join(gre_pssms_path, 'gre_pssm_%04d.json' % gre_id)
    if os.path.exists(path):
        with open(path) as infile:
            return json.load(infile)
    else:
        return {}

def __gre_motif_evalue(df, gre_id):
    df = df[df['GRE'] == gre_id].reset_index()
    result = df['motif e-value'][0]
    if np.isnan(result):
        result = 0.0
    return result

@app.route('/api/v1.0.0/corem_gres/<corem_id>')
def corem_gres(corem_id):
    corem_gres_path = app.config["COREM_GRES_FILE"]
    gre_summary_path = app.config["GRE_SUMMARY_FILE"]
    df = pd.read_csv(corem_gres_path)
    df = df[df['corem'] == int(corem_id)].reset_index()
    gre_df = pd.read_csv(gre_summary_path, sep='\t')

    gres = [{'gre': int(df['gre'][i]), 'q_value': df['qval_BH'][i],
                 'pssm': __gre_pssm(int(df['gre'][i])),
                 'motif_evalue': __gre_motif_evalue(gre_df, int(df['gre'][i]))}
                for i in range(df.shape[0])]
    return jsonify(gres=gres)


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
    """Returns a list of genes. The maximum number of entries returned is BATCH_SIZE.
    Optional request parameters:
    - start: start position
    - length: number of elements to return
    """
    start, num_entries, search = __request_batch_params()
    end = start + num_entries
    search_params = {}
    if search is not None and len(search) > 0:
        search_params = {'$or': [{'name': {'$regex': search}},  {'sysName': {'$regex': search}}] }

    chroms = db.genome.find({}, {"_id": 0, "scaffoldId": 1, "NCBI_RefSeq": 1 })
    chrom_map = { int(c['scaffoldId']): c['NCBI_RefSeq'] for c in chroms }
    cursor = db.row_info.find(search_params, {'_id': 0, 'row_id': 1, 'sysName': 1, 'name': 1,
                                        'accession': 1, 'desc': 1, 'start': 1, 'stop': 1,
                                        'strand': 1, 'scaffoldId': 1})[start:end]
    total = db.row_info.find(search_params).count()
    genes = [{ "id": r['row_id'],
                   "gene_name": r['sysName'],
                   "common_name": r['name'],
                   "accession":  str(r['accession']),
                   "description": r['desc'],
                   "start": r['start'], "stop": r['stop'], "strand": r['strand'],
                   "chromosome": chrom_map[r['scaffoldId']]} for r in cursor]

    return jsonify(genes=genes, total=total)


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


@app.route('/api/v1.0.0/gres')
def gres():
    try:
        start = int(request.args.get('start'))
    except Exception as e:
        start = 0
    try:
        num_entries = int(request.args.get('length'))
    except Exception as e:
        num_entries = BATCH_SIZE
    end = start + num_entries

    corem_gres_path = app.config["COREM_GRES_FILE"]
    df = pd.read_csv(corem_gres_path)
    gre_summary_path = app.config["GRE_SUMMARY_FILE"]
    gre_df = pd.read_csv(gre_summary_path, sep='\t')

    gres = [{
        'gre': i,
        'corems': list(df[df['gre'] == i]['corem']),
        'pssm': __gre_pssm(i),
        'motif_evalue': __gre_motif_evalue(gre_df, i)
        } for i in range(1, 164)]
    total = 163

    return jsonify(gres=gres[start:end], total=total)


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
