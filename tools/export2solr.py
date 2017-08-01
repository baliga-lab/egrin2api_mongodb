#!/usr/bin/env python3

"""
This script generates an import document for a Solr 6.5.0 server.
Our default location is garda.
It is setup to accept uploads through the post tool in Solr 6.5.0

The cores to hold the data are

  - mtb_corems
  - mtb_clusters

respectively
"""

import pymongo
from bson.objectid import ObjectId

def __gene_map(db):
    return { g['row_id']: g['sysName']
        for g in db.row_info.find({}, {'_id': 0, 'row_id': 1, 'sysName': 1, 'name': 1})}

def __cond_map(db):
    return { c['col_id']: c['egrin2_col_name']
        for c in db.col_info.find({}, {'_id': 0, 'col_id': 1, 'egrin2_col_name': 1})}

def export_biclusters(db):
    gene_map = __gene_map(db)
    cond_map = __cond_map(db)
    cursor = db.bicluster_info.find({}, {'_id': 1, 'rows': 1, 'columns': 1, 'residual': 1})  #[0:10]

    print("<add>")
    for c in cursor:
        genes = sorted([gene_map[rid] for rid in c["rows"] if rid in gene_map])
        conditions = sorted([cond_map[cid] for cid in c["columns"]])
        print("<doc>")
        print("  <field name=\"id\">%s</field>" % (c["_id"]))
        for gene in genes:
            print("  <field name=\"genes\">%s</field>" % gene)
        for cond in conditions:
            print("  <field name=\"conditions\">%s</field>" % cond)
        print("  <field name=\"residual\">%f</field>" % c['residual'])
        print("</doc>")
    print("</add>")


def export_corems(db):
    gene_map = __gene_map(db)
    cond_map = __cond_map(db)
    cursor = db.corem.find({}, {'_id': 0, 'corem_id': 1, 'rows': 1, 'cols': 1})
    print("<add>")
    for c in cursor:
        genes = sorted([gene_map[rid] for rid in c["rows"] if rid in gene_map])
        col_ids = [col['col_id'] for col in c['cols']]
        conditions = sorted([cond_map[cid] for cid in col_ids])
        print("<doc>")
        print("  <field name=\"id\">%s</field>" % (c["corem_id"]))
        for gene in genes:
            print("  <field name=\"genes\">%s</field>" % gene)
        for cond in conditions:
            print("  <field name=\"conditions\">%s</field>" % cond)
        print("</doc>")
    print("</add>")


if __name__ == '__main__':
    client = pymongo.MongoClient(host='como', port=27018)
    db = client['mtu_db']
    #export_biclusters(db)
    export_corems(db)
