#!/usr/bin/env python3

"""png_util_test.py
"""
import unittest
import xmlrunner
import sys
import json
from app import app


BATCHSIZE = 100


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

    def test_biclusters_no_args(self):
        """test the biclusters() function without arguments"""
        biclusters = json.loads(self.app.get('/api/v1.0.0/biclusters').data.decode('utf-8'))["biclusters"]
        # without arguments we currently receive the first BATCHSIZE biclusters
        self.assertEquals(len(biclusters), BATCHSIZE)

    def test_biclusters_with_args(self):
        """test the biclusters() function with start/length arguments"""
        biclusters = json.loads(self.app.get('/api/v1.0.0/biclusters?start=100&length=3').data.decode('utf-8'))["biclusters"]
        self.assertEquals(len(biclusters), 3)


if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(APITest))
    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(SUITE))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))