#!/usr/bin/env python3

"""png_util_test.py
"""
import unittest
import xmlrunner
import sys
import json
from app import app


class APITest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for png_util"""
    def setUp(self):
        self.app = app.test_client()

    def test_summary(self):
        """test the chunks() function"""
        summary = json.loads(self.app.get('/api/v1.0.0/summary').data.decode('utf-8'))
        self.assertTrue(summary['num_conditions'] > 0)
        self.assertTrue(summary['num_genes'] > 0)
        self.assertTrue(summary['num_corems'] > 0)
        self.assertTrue(summary['num_biclusters'] > 0)
        self.assertTrue(summary['num_gres'] > 0)


if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(APITest))
    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(SUITE))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
