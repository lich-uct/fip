#!/usr/bin/env python3
import unittest

from fiprofiling import relationmatrices as rm

class TestPrecomputedRMs(unittest.TestCase):
    def test_corm_import(self):
        from fiprofiling.data.precomputed_rm.drugbank.maccs import corm, coprm, pmirm, zpmirm
        self.assertIsInstance(corm, rm.CORM)
        self.assertIsInstance(coprm, rm.COPRM)
        self.assertIsInstance(pmirm, rm.PMIRM)
        self.assertIsInstance(zpmirm, rm.ZPMIRM)

if __name__ == '__main__':
    unittest.main()
