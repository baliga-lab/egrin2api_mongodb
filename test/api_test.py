#!/usr/bin/env python3

"""png_util_test.py
"""
import unittest
import xmlrunner
import sys
import json
from app import app


BATCHSIZE = 100
BICLUSTER1_ID = '552ef64cb744335b89d3d6b6'

class APITest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for png_util"""
    def setUp(self):
        self.app = app.test_client()

    def test_summary(self):
        """test the summary() function"""
        summary = json.loads(self.app.get('/api/v1.0.0/summary').data.decode('utf-8'))
        self.assertTrue(summary['num_conditions'] > 0)
        self.assertTrue(summary['num_genes'] > 0)
        self.assertTrue(summary['num_corems'] > 0)
        self.assertTrue(summary['num_biclusters'] > 0)
        self.assertTrue(summary['num_gres'] > 0)

    def test_api_info(self):
        """test the api_info() function"""
        info = json.loads(self.app.get('/api/v1.0.0/api_info').data.decode('utf-8'))
        self.assertEquals(info['model_type'], 'egrin2')
        self.assertEquals(info['api_version'], '1.0.0')
        self.assertEquals(info['organism'], 'mtu')
        self.assertEquals(info['species'], 'Mycobacterium_tuberculosis_H37Rv')
        self.assertEquals(info['num_runs'], 278)


    def test_conditions_no_args(self):
        """test the conditions() function without arguments"""
        conditions = json.loads(self.app.get('/api/v1.0.0/conditions').data.decode('utf-8'))["conditions"]
        # without arguments we currently receive all conditions
        self.assertEquals(len(conditions), 1753)


    def test_corems_no_args(self):
        """test the corems() function without arguments"""
        corems = json.loads(self.app.get('/api/v1.0.0/corems').data.decode('utf-8'))["corems"]
        # without arguments we currently receive all corems
        self.assertEquals(len(corems), 560)

    def test_genes_no_args(self):
        """test the genes() function without arguments"""
        genes = json.loads(self.app.get('/api/v1.0.0/genes').data.decode('utf-8'))["genes"]
        # without arguments we currently receive the first BATCHSIZE genes
        self.assertEquals(len(genes), BATCHSIZE)

    def test_genes_with_args(self):
        """test the genes() function with start/length arguments"""
        genes = json.loads(self.app.get('/api/v1.0.0/genes?start=100&length=2').data.decode('utf-8'))["genes"]
        self.assertEquals(len(genes), 2)

    def test_genes_with_search(self):
        """test the genes() function with start/length arguments"""
        result = json.loads(self.app.get('/api/v1.0.0/genes?search=dna').data.decode('utf-8'))
        genes = result['genes']
        total = result['total']
        self.assertEquals(total, 12)
        self.assertTrue(len(genes) <= 12)

    def test_biclusters_no_args(self):
        """test the biclusters() function without arguments"""
        biclusters = json.loads(self.app.get('/api/v1.0.0/biclusters').data.decode('utf-8'))["biclusters"]
        # without arguments we currently receive the first BATCHSIZE biclusters
        self.assertEquals(len(biclusters), BATCHSIZE)

    def test_biclusters_with_args(self):
        """test the biclusters() function with start/length arguments"""
        biclusters = json.loads(self.app.get('/api/v1.0.0/biclusters?start=100&length=3').data.decode('utf-8'))["biclusters"]
        self.assertEquals(len(biclusters), 3)

    def test_corem_conditions_1(self):
        """test the corem_conditions() function with corem 1"""
        conds = json.loads(self.app.get('/api/v1.0.0/corem_conditions/1').data.decode('utf-8'))["conditions"]
        self.assertEquals(len(conds), 659)

    def test_corem_genes_1(self):
        """test the corem_genes() function with corem 1"""
        genes = json.loads(self.app.get('/api/v1.0.0/corem_genes/1').data.decode('utf-8'))["genes"]
        self.assertEquals(len(genes), 3)

    def test_corem_genes_2_search(self):
        """test the corem_genes() function with corem 1"""
        result = json.loads(self.app.get('/api/v1.0.0/corem_genes/2?search=Rv13').data.decode('utf-8'))
        genes = result['genes']
        total = result['total']
        self.assertEquals(len(genes), 1)
        self.assertEquals(total, 1)

    def test_bicluster_pssms_1(self):
        """test the bicluster_pssms() function with the first bicluster"""
        motifs = json.loads(self.app.get('/api/v1.0.0/bicluster_pssms/' + BICLUSTER1_ID).data.decode('utf-8'))["motifs"]
        self.assertEquals(len(motifs), 2)

    def test_bicluster_conditions_1(self):
        """test the bicluster_conditions() function with the first bicluster"""
        conds = json.loads(self.app.get('/api/v1.0.0/bicluster_conditions/' + BICLUSTER1_ID).data.decode('utf-8'))["conditions"]
        self.assertEquals(len(conds), 91)

    def test_bicluster_genes_1(self):
        """test the bicluster_genes() function with the first bicluster"""
        genes = json.loads(self.app.get('/api/v1.0.0/bicluster_genes/' + BICLUSTER1_ID).data.decode('utf-8'))["genes"]
        self.assertEquals(len(genes), 16)

    def test_bicluster_info_1(self):
        """test the bicluster_info() function with the first bicluster"""
        bicluster = json.loads(self.app.get('/api/v1.0.0/bicluster_info/' + BICLUSTER1_ID).data.decode('utf-8'))["bicluster"]
        self.assertEquals(bicluster['id'], BICLUSTER1_ID)
        self.assertEquals(bicluster['num_conditions'], 91)
        self.assertEquals(bicluster['num_genes'], 16)
        self.assertAlmostEquals(bicluster['residual'], 0.4244137794133839)

    def test_condition_biclusters_1(self):
        """test the condition_biclusters() function with the first condition"""
        biclusters = json.loads(self.app.get('/api/v1.0.0/condition_biclusters/1').data.decode('utf-8'))["biclusters"]
        self.assertEquals(len(biclusters), 6235)

    def test_gene_biclusters_Rv0001(self):
        """test the gene_biclusters() function with a specified gene"""
        biclusters = json.loads(self.app.get('/api/v1.0.0/gene_biclusters/Rv0001').data.decode('utf-8'))["biclusters"]
        self.assertEquals(len(biclusters), 556)

    def test_gene_info_Rv0001(self):
        """test the gene_info() function with a specified gene"""
        gene = json.loads(self.app.get('/api/v1.0.0/gene_info/Rv0001').data.decode('utf-8'))["gene"]
        self.assertEquals(gene['accession'], 'NP_214515.1')
        self.assertEquals(gene['chromosome'], 'NC_000962')
        self.assertEquals(gene['common_name'], 'dnaA')
        self.assertEquals(gene['gene_name'], 'Rv0001')
        self.assertEquals(gene['id'], 0)
        self.assertEquals(gene['start'], 1)
        self.assertEquals(gene['stop'], 1524)
        self.assertEquals(gene['strand'], '+')
        self.assertEquals(gene['description'], 'chromosomal replication initiation protein (NCBI)')

    def test_condition_info_1(self):
        """test the condition_info() function with the first condition"""
        info = json.loads(self.app.get('/api/v1.0.0/condition_info/1').data.decode('utf-8'))
        self.assertEquals(len(info['blocks']['carbon.source.fatty.acid']), 224)

    def test_corem_info_1(self):
        """test the corem_info() function with the first corem"""
        info = json.loads(self.app.get('/api/v1.0.0/corem_info/1').data.decode('utf-8'))
        self.assertEquals(info['corem'], "1")
        self.assertEquals(len(info['genes']), 3)
        self.assertEquals(len(info['gres']), 51)

    def test_gene_gres_Rv0116c(self):
        """test the gene_gres() function with the specified gene"""
        gres = json.loads(self.app.get('/api/v1.0.0/gene_gres/Rv0116c').data.decode('utf-8'))
        self.assertEquals(gres['gene'], 'Rv0116c')
        self.assertEquals(len(gres['gres']), 6)

    def test_corems_with_gene_Rv1100(self):
        """test the corems_with_gene() function with the specified gene"""
        corems = json.loads(self.app.get('/api/v1.0.0/corems_with_gene/Rv1100').data.decode('utf-8'))['corem_infos']
        self.assertEquals(len(corems), 7)

    def test_corem_expressions1(self):
        """test the corem_expressions() function with the specified corem"""
        exps = json.loads(self.app.get('/api/v1.0.0/corem_expressions/1').data.decode('utf-8'))
        num_conditions = len(exps['conditions'])
        self.assertTrue(num_conditions > 0)
        self.assertEquals(len(exps['expressions']['Rv0440']), num_conditions)


    def test_corem_condition_enrichment1(self):
        """test the corem_condition_enrichment() function with the specified corem"""
        exps = json.loads(self.app.get('/api/v1.0.0/corem_condition_enrichment/1').data.decode('utf-8'))
        num_blocks = len(exps['condition_blocks'])
        self.assertEquals(num_blocks, 3)

    def test_corem_category_enrichment1(self):
        """test the corem_condition_enrichment() function with the specified corem"""
        cats = json.loads(self.app.get('/api/v1.0.0/corem_categories/1').data.decode('utf-8'))
        num_cats = len(cats['categories'])
        self.assertEquals(num_cats, 1)

    def test_corem_gres2(self):
        """test the corem_gres() function with the specified corem"""
        gres = json.loads(self.app.get('/api/v1.0.0/corem_gres/2').data.decode('utf-8'))['gres']
        self.assertEquals(len(gres), 2)
        self.assertTrue(len(gres[0]['pssm']['values']) > 0)

    def test_gene_gre_counts(self):
        """test the gene_gre_counts() function with the specified gene"""
        gre_counts = json.loads(self.app.get('/api/v1.0.0/gene_gres/Rv0135c').data.decode('utf-8'))
        self.assertEquals(len(gre_counts['chipseq_peaks']), 3)
        self.assertTrue(len(gre_counts['gres']['GRE_1']) > 0)


if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(APITest))
    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(SUITE))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
