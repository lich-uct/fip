#!/usr/bin/env python3
import unittest

from fiprofiling import util

class TestBitArrays(unittest.TestCase):
    def test_bool2binstring(self):
        boollist = [True, False, True, True, False]
        binstring = util.bool2binstring(boollist)
        boollist2 = util.binstring2bool(binstring)

    def test_bool2string_single(self):
        self.assertEqual(util.bool2binstring(True), '1')
        self.assertEqual(util.bool2binstring(False), '0')

if __name__ == '__main__':
    unittest.main()
